# ============================================================================
# main_figures.R
# ----------------------------------------------------------------------------
# Exact reproduction of the five main-text figures in
#   Claessen, Traber & Schoonvelde, "Moderation or Amplification? Mapping
#   Speech Distinctiveness Among Green Parties Using a Machine Learning
#   Approach", Legislative Studies Quarterly.
#
# The plotting code below is copied verbatim from the analysis pipeline
# (Salience_Graphs_Sentences.R, Se_Gr_disti.R, Se_Models.R), so the output is
# byte-for-byte the same as the figures in the article — not a look-alike
# reconstruction. Every figure is built from text-free, sentence-level data, so
# it is auditable down to individual sentences:
#
#   data/green_party_sentences_xl.rds     sentence-level Green-party data; Fig 2
#                                         aggregates it and Figs 3-5 refit the
#                                         logit models from it
#   data/salience_sentences_xl.csv        one row per Green env/welfare sentence
#   data/salience_monthly_totals_xl.csv   monthly Green sentence totals
#                                         (Fig 1 = numerator / denominator of these)
#
# Outputs (PDF + PNG) are written to output/:
#   figure1_salience.pdf
#   figure2_distinctiveness_bars.pdf
#   figure3_odds_ratios_environment.pdf / _welfare.pdf
#   figure4_predicted_gov_opp.pdf
#   figure5_electoral_cycle_environment.pdf / _welfare.pdf
#
# Run from the repository root:  Rscript main_figures.R
#
# NOTE ON R VERSION: tested on R 4.5.2 and R 4.6.0; either reproduces the
# figures in the article. The published figures were themselves made under two R
# versions (Figs 1-2 under 4.5.2, Figs 3-5 under 4.6.0), so which subset comes
# out byte-identical depends on the R you run -- under 4.5.2 Figs 1-2 are
# byte-identical, under 4.6.0 Figs 3-5 are. The rest differ only by sub-pixel
# text anti-aliasing from the ggplot2 version (mean diff < 0.4/255, at glyph
# edges); the plotted geometry, colours and legends are unchanged.
# ============================================================================

# ---- preflight ---------------------------------------------------------------
# The Figure 3 coefficient plots use sjPlot::plot_models, which needs
# insight >= 1.5.0.6. Two failure modes, both of which otherwise surface as a
# cryptic "Error in purrr::map()" from deep inside sjPlot:
#   (a) insight is too old on disk            -> install.packages("insight")
#   (b) an older insight is already LOADED in a long-running session (it was
#       loaded before the package was updated) -> RESTART R; a loaded namespace
#       cannot be upgraded in place.
# main_models.R does not use sjPlot and is unaffected by either.
local({
  need <- package_version("1.5.0.6")
  if ("insight" %in% loadedNamespaces() &&
      package_version(getNamespaceVersion("insight")) < need) {
    stop("insight ", getNamespaceVersion("insight"), " is already loaded in this ",
         "session, but Figure 3 (sjPlot) needs >= ", need, ".\n",
         "  Restart R and source this script again -- a namespace that is already ",
         "loaded cannot be upgraded in place.", call. = FALSE)
  }
  inst <- tryCatch(packageVersion("insight"), error = function(e) NULL)
  if (is.null(inst) || inst < need) {
    stop("Figure 3 (sjPlot) needs insight >= ", need, " (found: ",
         if (is.null(inst)) "not installed" else format(inst), ").\n",
         "  Run install.packages(\"insight\"), then restart R.", call. = FALSE)
  }
})

suppressWarnings(suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(sjPlot)
  library(sandwich)
  library(lmtest)
  library(marginaleffects)
  library(scales)
  library(grid)
  library(forcats)
}))

dir.create("output", showWarnings = FALSE)

green <- as.data.frame(readRDS("data/green_party_sentences_xl.rds"))
# The pipeline splits the data by a country_topic key; rebuild it (identical to
# the pipeline's "DE_environment" / "IR_welfare" / ... labels).
green$country_topic_model <- paste(green$country, green$topic, sep = "_")

# ============================================================================
# SHARED HELPERS  (verbatim from Modelling_Helpers.R / Se_Models.R)
# ============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

# One-way speaker cluster vector, NA-safe.
make_speaker_cluster <- function(df, speaker_col = "speaker") {
  if (speaker_col %in% names(df) && !all(is.na(df[[speaker_col]]))) return(df[[speaker_col]])
  seq_len(nrow(df))
}

# Pre-compute HC0 cluster-robust vcov (clustered by speaker) once per model.
build_vcov_cache <- function(models, df, type = "HC0") {
  cl <- make_speaker_cluster(df)
  cache <- vector("list", length(models)); names(cache) <- names(models)
  for (nm in names(models)) {
    m <- models[[nm]]
    cache[[nm]] <- list(
      key  = paste(deparse(stats::formula(m)), collapse = " "),
      vcov = tryCatch(sandwich::vcovCL(m, cluster = cl, type = type),
                      error = function(e) sandwich::vcovCL(m, cluster = cl, type = "HC0"))
    )
  }
  attr(cache, "df") <- df; attr(cache, "type") <- type
  cache
}

# Closure mapping a model -> cached vcov (for sjPlot's vcov.fun); cache misses
# (e.g. the refitted electoral-cycle model) fall through to a fresh HC0.
make_vcov_lookup <- function(cache) {
  cache_keys <- vapply(cache, function(e) e$key %||% "", character(1))
  cache_df   <- attr(cache, "df")
  cache_type <- attr(cache, "type") %||% "HC0"
  function(model, ...) {
    key <- paste(deparse(stats::formula(model)), collapse = " ")
    idx <- match(key, cache_keys)
    if (!is.na(idx)) return(cache[[idx]]$vcov)
    cl <- make_speaker_cluster(cache_df)
    tryCatch(sandwich::vcovCL(model, cluster = cl, type = cache_type),
             error = function(e) sandwich::vcovCL(model, cluster = cl, type = "HC0"))
  }
}

# marginaleffects-based robust predictions (honour the cluster-robust vcov).
predict_robust <- function(model, terms, vcov, n = 50) {
  data_used <- if (!is.null(model$model)) model$model else model$data
  dg_args <- list()
  for (v in terms) {
    col <- data_used[[v]]
    if (is.null(col)) stop(sprintf("predict_robust: variable '%s' not in model frame", v))
    if (is.factor(col)) dg_args[[v]] <- levels(col)
    else                dg_args[[v]] <- seq(min(col, na.rm = TRUE), max(col, na.rm = TRUE), length.out = n)
  }
  nd    <- do.call(marginaleffects::datagrid, c(list(model = model), dg_args))
  preds <- marginaleffects::predictions(model, newdata = nd, vcov = vcov)
  pdf   <- as.data.frame(preds)
  std_err <- if ("std.error" %in% names(pdf)) pdf$std.error
             else (pdf$conf.high - pdf$conf.low) / (2 * stats::qnorm(0.975))
  out <- data.frame(
    x = pdf[[terms[1]]], predicted = pdf$estimate, std.error = std_err,
    conf.low = pdf$conf.low, conf.high = pdf$conf.high, stringsAsFactors = FALSE
  )
  if (length(terms) >= 2) out$group <- pdf[[terms[2]]]
  out
}

# Fit the three nested logit models used throughout the paper.
fit_three <- function(df, y) {
  f_base <- as.formula(sprintf("%s ~ govopp", y))
  f_controls <- as.formula(sprintf("%s ~ govopp + non_green_issue_salience + non_green_distinctiveness_avg + seats_perc + log_sentence_nchar + election + year", y))
  f_time <- as.formula(sprintf("%s ~ govopp + time_until_election_neg + I(time_until_election_neg^2) +
                                    non_green_issue_salience + non_green_distinctiveness_avg +
                                    seats_perc + log_sentence_nchar + election + year", y))
  m1 <- glm(f_base, data = df, family = binomial(link = 'logit'), x = TRUE)
  m2 <- glm(f_controls, data = df, family = binomial(link = 'logit'), x = TRUE)
  m3 <- glm(f_time, data = df, family = binomial(link = 'logit'), x = TRUE)
  list(m1 = m1, m2 = m2, m3 = m3)
}

# ============================================================================
# FIGURE 1 — Issue salience over time
#   Plot code verbatim from Salience_Graphs_Sentences.R; the monthly salience is
#   computed here from the sentence-level coding so the figure is auditable down
#   to individual sentences (no raw text needed):
#
#     percentage(topic, month) = 100 * (# Green sentences of that topic that month)
#                                     / (total Green sentences that month, all topics)
#
#   Numerator: data/salience_sentences_xl.csv  — one row per Green environment /
#     welfare sentence (country, date, topic), full corpus date range.
#   Denominator: data/salience_monthly_totals_xl.csv — total Green sentences per
#     month across ALL topics (the salience denominator).
#   This aggregation reproduces the published monthly percentages exactly.
# ============================================================================
cat("\n[Figure 1] issue salience over time (from sentence-level coding) ...\n")

break.vec.de <- c(as.Date("1990-01-01"),
                  seq(from = as.Date("1990-01-01"), to = as.Date("2015-01-01"), by = "5 years"),
                  as.Date("2015-01-01"))
break.vec.ir <- c(as.Date("1990-01-01"),
                  seq(from = as.Date("1990-01-01"), to = as.Date("2010-01-01"), by = "5 years"),
                  as.Date("2010-01-01"))

sal_sent <- read.csv("data/salience_sentences_xl.csv", stringsAsFactors = FALSE)
sal_tot  <- read.csv("data/salience_monthly_totals_xl.csv", stringsAsFactors = FALSE)
sal_sent$date  <- as.Date(sal_sent$date)
sal_sent$year  <- as.integer(format(sal_sent$date, "%Y"))
sal_sent$month <- as.integer(format(sal_sent$date, "%m"))
salience <- sal_sent %>%
  dplyr::count(country, topic, year, month, name = "sum_topic") %>%
  dplyr::left_join(sal_tot, by = c("country", "year", "month")) %>%
  dplyr::mutate(percentage = (sum_topic / total_sentences) * 100,
                date  = as.Date(sprintf("%04d-%02d-01", year, month)),
                topic = tools::toTitleCase(topic))
data_de <- salience[salience$country == "DE", ]
data_ir <- salience[salience$country == "IR", ]

p_de <- ggplot(data_de, aes(x = date, y = percentage, color = topic)) +
  geom_smooth(aes(fill = topic), alpha = 0.3, span = 0.3, method = "loess", size = 1.15) +
  scale_color_manual(values = c("Environment" = "darkgreen", "Welfare" = "#762A83")) +
  scale_fill_manual(values = c("Environment" = "darkgreen", "Welfare" = "#762A83")) +
  coord_cartesian(ylim = c(0, 12.5)) +
  ggplot2::annotate("rect", xmin = as.Date("1998-10-27"), xmax = as.Date("2005-11-22"),
                    ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) +
  ggplot2::annotate(geom = "label", x = as.Date("2002-06-05"), y = 11,
                    label = "Greens in\nGovernment", fontface = 'italic', size = 4.5) +
  theme_classic() +
  scale_x_date(breaks = break.vec.de, date_labels = '%Y', expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  labs(x = "Time", y = "Percentage of Total Sentences per Party per Month",
       title = paste("German Green Party"), color = NULL, fill = NULL) +
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = 'bold')) +
  guides(colour = guide_legend(nrow = 1))

p_ir <- ggplot(data_ir, aes(x = date, y = percentage, color = topic)) +
  geom_smooth(aes(fill = topic), alpha = 0.3, span = 0.3, method = "loess", size = 1.15) +
  scale_color_manual(values = c("Environment" = "darkgreen", "Welfare" = "#762A83")) +
  scale_fill_manual(values = c("Environment" = "darkgreen", "Welfare" = "#762A83")) +
  coord_cartesian(ylim = c(0, 12.5)) +
  ggplot2::annotate("rect", xmin = as.Date("2007-06-14"), xmax = as.Date("2011-03-09"),
                    ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.3) +
  ggplot2::annotate(geom = "label", x = as.Date("2008-11-01"), y = 11,
                    label = "Greens in\nGovernment", fontface = 'italic', size = 4.5) +
  theme_classic() +
  scale_x_date(breaks = break.vec.ir, date_labels = '%Y', expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  labs(x = "Time", y = "Percentage of Total Sentences per Party per Month",
       title = paste("Irish Green Party"), color = NULL, fill = NULL) +
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = 'bold')) +
  guides(colour = guide_legend(nrow = 1))

fig1 <- ggarrange(p_de, p_ir, ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave("output/figure1_salience.pdf", fig1, width = 11.2, height = 5.3)
ggsave("output/figure1_salience.png", fig1, width = 11.2, height = 5.3, dpi = 300)

# ============================================================================
# FIGURE 2 — Distinctiveness by government / opposition
#   Plot code verbatim from Se_Gr_disti.R; the bars are computed here from the
#   sentence-level data, so the figure is auditable down to individual sentences.
#   Each bar is the mean of correct_SVC (1 = the SVC classifier attributed the
#   sentence to the Green Party) over the Green sentences of that issue while the
#   party was in government vs. opposition. Government/opposition is assigned from
#   the sitting date exactly as in the pipeline (Germany: in government 1998-10-27
#   to 2005-11-01, and post-2005 opposition is excluded from the figure; Ireland:
#   in government 2007-06-14 to 2011-03-09).
# ============================================================================
cat("[Figure 2] distinctiveness by government/opposition (from sentence-level data) ...\n")

.bars_de <- green %>%
  dplyr::filter(country == "DE") %>%
  dplyr::mutate(date = as.Date(date), issue = topic,
                govopp2 = dplyr::case_when(
                  date < as.Date("2005-11-01") & date > as.Date("1998-10-27") ~ "Government",
                  date > as.Date("2005-11-01")                                ~ "opposition after",
                  TRUE                                                        ~ "Opposition")) %>%
  dplyr::filter(govopp2 != "opposition after") %>%
  dplyr::group_by(issue, govopp2) %>%
  dplyr::summarise(mean = mean(correct_SVC, na.rm = TRUE), sd.value = sd(correct_SVC, na.rm = TRUE),
                   count = dplyr::n(), se.mean = sd.value / sqrt(count), .groups = "drop")
.bars_de$govopp2 <- forcats::fct_reorder(.bars_de$govopp2, .bars_de$count)
error_bars_de.long2 <- as.data.frame(.bars_de)

.bars_ir <- green %>%
  dplyr::filter(country == "IR") %>%
  dplyr::mutate(date = as.Date(date), issue = topic,
                govopp = ifelse(date < as.Date("2007-06-14") | date > as.Date("2011-03-09"),
                                "Opposition", "Government")) %>%
  dplyr::group_by(issue, govopp) %>%
  dplyr::summarise(mean = mean(correct_SVC, na.rm = TRUE), sd.value = sd(correct_SVC, na.rm = TRUE),
                   count = dplyr::n(), se.mean = sd.value / sqrt(count), .groups = "drop")
.bars_ir$govopp <- factor(.bars_ir$govopp, levels = c("Opposition", "Government"))
error_bars_IR.long2 <- as.data.frame(.bars_ir)

gs1 <- ggplot(error_bars_de.long2, aes(x = govopp2, y = mean, colour = issue, fill = issue)) +
  geom_col(aes(color = issue, fill = issue), position = 'dodge', width = .5) +
  ylab("Distinctiveness") + theme_classic() + xlab('') +
  scale_fill_manual(values = c("environment" = "darkgreen", "welfare" = "#762A83"),
                    labels = c("environment" = "Environment", "welfare" = "Welfare")) +
  scale_colour_manual(values = c("environment" = "darkgreen", "welfare" = "#762A83"),
                      labels = c("environment" = "Environment", "welfare" = "Welfare")) +
  coord_cartesian(ylim = c(0, 0.71)) +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle('German Green Party') +
  guides(fill = guide_legend(title = NULL), colour = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0.000, 0.000))

gs1.1 <- gs1 +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10), legend.title = element_text(size = 11),
        plot.title = element_text(hjust = .5, size = 14, face = 'bold'),
        strip.text.x = element_text(face = 'bold'),
        strip.background = element_rect(fill = '#EEEEEE'), legend.position = 'bottom')

gs2 <- ggplot(error_bars_IR.long2, aes(x = govopp, y = mean, colour = issue, fill = issue)) +
  geom_col(aes(color = issue, fill = issue), position = 'dodge', width = .5) +
  ylab("Distinctiveness") + theme_classic() + xlab('') +
  scale_fill_manual(values = c("environment" = "darkgreen", "welfare" = "#762A83"),
                    labels = c("environment" = "Environment", "welfare" = "Welfare")) +
  scale_colour_manual(values = c("environment" = "darkgreen", "welfare" = "#762A83"),
                      labels = c("environment" = "Environment", "welfare" = "Welfare")) +
  coord_cartesian(ylim = c(0, 0.71)) +
  scale_x_discrete(guide = guide_axis(angle = 0)) +
  ggtitle('Irish Green Party') +
  guides(fill = guide_legend(title = NULL), colour = guide_legend(title = NULL)) +
  scale_y_continuous(expand = c(0.000, 0.000))

gs3 <- gs2 +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10), legend.title = element_text(size = 11),
        plot.title = element_text(hjust = .5, size = 14, face = 'bold'),
        strip.text.x = element_text(face = 'bold'),
        strip.background = element_rect(fill = '#EEEEEE'), legend.position = 'bottom')

fig2 <- ggarrange(gs1.1, gs3, nrow = 1, common.legend = TRUE, legend = 'bottom')
ggsave("output/figure2_distinctiveness_bars.pdf", fig2, width = 12, height = 6)
ggsave("output/figure2_distinctiveness_bars.png", fig2, width = 12, height = 6, dpi = 300)

# ============================================================================
# FIGURES 3-5 — refit the logit models and draw the coefficient, predicted-
# probability and electoral-cycle plots (verbatim from Se_Models.R)
# ============================================================================

# --- Coefficient plot (Fig 3) --------------------------------------------------
make_coefplot_exact <- function(models, vcov_lookup, title_, topic, country) {
  outcome_var <- as.character(models$m1$formula[2])
  rm_terms <- c("seats_perc", "election", "terms", "log_terms", "year",
                "sentence_nchar", "log_sentence_nchar",
                "time_until_election_neg", "I(time_until_election_neg^2)",
                "non_green_issue_salience", "non_green_distinctiveness_avg")
  all_term_names <- unique(unlist(lapply(models, function(m) names(stats::coef(m)))))
  rm_terms <- unique(c(rm_terms, grep("^election", all_term_names, value = TRUE)))

  plot_95 <- sjPlot::plot_models(
    models$m1, models$m2, models$m3,
    vcov.fun = vcov_lookup, vcov.args = list(),
    auto.label = TRUE, show.values = TRUE, vline.color = "white",
    spacing = 0.65, dot.size = 2.75, line.size = 0.55,
    p.threshold = c(0.05, 0.01, 0.001), rm.terms = rm_terms
  )
  patch_plus <- function(d) {
    if ("p.stars" %in% names(d)) d$p.stars <- as.character(d$p.stars)
    if ("p.label" %in% names(d)) d$p.label <- as.character(d$p.label)
    idx <- !is.na(d$p.value) & d$p.value < 0.1 &
           "p.stars" %in% names(d) &
           (is.na(d$p.stars) | trimws(d$p.stars) %in% c("", "n.s."))
    if (any(idx)) {
      d$p.stars[idx] <- "+"
      if ("p.label" %in% names(d)) d$p.label[idx] <- paste0(trimws(d$p.label[idx]), " +")
    }
    d
  }
  plot_95$data <- patch_plus(plot_95$data)

  plot_90 <- sjPlot::plot_models(
    models$m1, models$m2, models$m3,
    vcov.fun = vcov_lookup, vcov.args = list(),
    auto.label = TRUE, show.values = TRUE, ci.lvl = .90,
    vline.color = "white", spacing = 0.65, dot.size = 2.75, line.size = 0.75,
    p.threshold = c(0.05, 0.01, 0.001), rm.terms = rm_terms
  )
  plot_90$data <- patch_plus(plot_90$data)

  plot_95$data$group <- dplyr::recode(plot_95$data$group,
    !!paste0(outcome_var, '.3') := 'with controls (incl. time)',
    !!paste0(outcome_var, '.2') := 'with controls',
    !!paste0(outcome_var, '.1') := 'base')
  plot_90$data$group <- dplyr::recode(plot_90$data$group,
    !!paste0(outcome_var, '.3') := 'with controls (incl. time)',
    !!paste0(outcome_var, '.2') := 'with controls',
    !!paste0(outcome_var, '.1') := 'base')

  plot_95$data$term <- dplyr::recode(plot_95$data$term,
    "govoppGovernment" = if (topic == "welfare") "Government" else "In Government")

  conf_90 <- dplyr::select(plot_90$data, conf.low, conf.high)
  colnames(conf_90) <- c('conf.low90', 'conf.high90')
  plot_95$data <- cbind(plot_95$data, conf_90)

  all_ci_values <- c(plot_95$data$conf.low, plot_95$data$conf.high,
                     plot_95$data$conf.low90, plot_95$data$conf.high90)
  min_val <- min(all_ci_values, na.rm = TRUE); max_val <- max(all_ci_values, na.rm = TRUE)
  padding <- (max_val - min_val) * 0.1
  ylim <- c(max(0.1, min_val - padding), max_val + padding)

  plot_95 +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = 'dashed', col = '#5A5A5A') +
    ylab("Odds Ratios") + labs(color = "Models", title = title_) +
    theme(axis.text = element_text(size = 13), axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13), legend.text = element_text(size = 13),
          legend.title = element_text(size = 13),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0.5),
          legend.position = "right") +
    coord_flip(ylim = ylim) +
    geom_errorbar(aes_string(ymin = "conf.low", ymax = "conf.high"),
                  width = 0.25, position = position_dodge(0.65)) +
    geom_linerange(aes_string(ymin = "conf.low90", ymax = "conf.high90"),
                   size = 1.5, position = position_dodge(0.65))
}

# --- Predicted probability plot (individual; Fig 4 is the combined version) ----
make_pred_plot_exact <- function(models, df, vcov_lookup, title_, topic) {
  color_ <- if (topic == "welfare") "#762A83" else "darkgreen"
  V_m3 <- vcov_lookup(models$m3)
  dat3 <- predict_robust(models$m3, terms = "govopp", vcov = V_m3)
  ggplot(dat3, aes(x = x, y = predicted)) +
    geom_point(color = color_, size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, linewidth = 0.75, color = color_) +
    theme_classic() + ylab("Predicted Values of Distinctiveness") +
    labs(color = "Topic", title = title_) +
    theme(axis.text = element_text(size = 14), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0), axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 14), legend.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0.5), legend.position = "right") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = .2), labels = seq(0, 1, by = .2))
}

# --- Electoral-cycle plot (Fig 5) ---------------------------------------------
make_electoral_cycle_plot <- function(df, outcome_var, vcov_lookup, title_, topic) {
  color_ <- if (topic == "welfare") "#762A83" else "darkgreen"
  formula_time <- as.formula(sprintf("%s ~ time_until_election_neg + I(time_until_election_neg^2) +
                                      govopp + non_green_issue_salience + non_green_distinctiveness_avg +
                                      seats_perc + log_sentence_nchar + election + year", outcome_var))
  model_time <- glm(formula_time, data = df, family = binomial(link = 'logit'), x = TRUE)
  model_time$call$data <- df
  V_time <- vcov_lookup(model_time)
  dat <- predict_robust(model_time, terms = "time_until_election_neg", vcov = V_time, n = 100)
  ylim <- c(0, 1); breaks <- seq(0, 1, by = 0.2)
  electoral_plot <- ggplot(dat, aes(x = x, y = predicted)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = color_, alpha = 0.2) +
    geom_line(color = color_, linewidth = 1) +
    geom_rug(data = df, aes(x = time_until_election_neg), inherit.aes = FALSE,
             alpha = 0.3, length = unit(0.02, "npc"), sides = "b", color = color_) +
    theme_classic() + ylab("Predicted Distinctiveness Values") +
    xlab("Time until Election in Years") +
    theme(axis.text = element_text(size = 14), axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14), legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          plot.title = element_text(size = 17, face = 'bold', hjust = 0.5),
          legend.position = "right", panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), strip.text.x = element_text(face = 'bold', size = 14),
          strip.background = element_rect(fill = '#EEEEEE')) +
    ggtitle(title_) + coord_cartesian(ylim = ylim) +
    scale_y_continuous(breaks = breaks, labels = breaks, expand = c(0.001, 0.001)) +
    scale_x_continuous(expand = c(0, 0))
  electoral_plot$data$group <- title_
  electoral_plot
}

# --- Main loop: fit models and build the per-country/topic plots --------------
plots <- list()
ctm_list <- unique(green$country_topic_model)
data_by_ctm <- split(green, green$country_topic_model)

for (ctm in ctm_list) {
  df_ctm  <- data_by_ctm[[ctm]]
  country <- substr(ctm, 1, 2)
  topic   <- sub(".*_", "", ctm)

  # German welfare speech before the Greens entered government (identical to pipeline).
  if (country == "DE" && topic == "welfare") {
    df_ctm <- df_ctm %>% dplyr::filter(date < '2005-11-22')
  }

  y <- "correct_SVC"
  work <- df_ctm[!is.na(df_ctm[[y]]), ]
  work$govopp <- droplevels(work$govopp)

  cat(sprintf("  [%s] fitting logit models on %d sentences ...\n", ctm, nrow(work)))
  mods <- fit_three(work, y)

  vcov_cache  <- build_vcov_cache(mods, work, type = "HC0")
  vcov_lookup <- make_vcov_lookup(vcov_cache)

  country_name <- if (country == "DE") "German" else "Irish"
  plot_title <- paste(country_name, "Green Party")

  mods$m3$call$data <- work
  dat3 <- predict_robust(mods$m3, terms = "govopp", vcov = vcov_cache$m3$vcov)
  dat3 <- dat3 %>% dplyr::filter(!is.na(x)) %>% dplyr::mutate(x = droplevels(as.factor(x)))
  dat3$topic <- topic; dat3$country <- country

  coef_plot      <- make_coefplot_exact(mods, vcov_lookup, plot_title, topic, country)
  pred_plot      <- make_pred_plot_exact(mods, work, vcov_lookup, plot_title, topic)
  electoral_plot <- make_electoral_cycle_plot(work, y, vcov_lookup, plot_title, topic)

  plots[[ctm]] <- list(coef = coef_plot, pred = pred_plot, pred_data = dat3, electoral = electoral_plot)
}

# --- FIGURE 3: combined coefficient plots -------------------------------------
cat("[Figure 3] odds-ratio coefficient plots ...\n")
env_coef <- list(); wel_coef <- list()
for (ctm in names(plots)) {
  side <- if (grepl("DE", ctm)) "DE" else "IR"
  if (grepl("environment", ctm)) env_coef[[side]] <- plots[[ctm]]$coef
  if (grepl("welfare", ctm))     wel_coef[[side]] <- plots[[ctm]]$coef
}
combined_env_coef <- ggarrange(env_coef$DE, env_coef$IR, nrow = 1, widths = c(1, 1),
                               common.legend = TRUE, legend = 'bottom')
combined_env_coef <- annotate_figure(combined_env_coef,
  top = text_grob("Issue: Environment", color = "black", face = "bold", size = 22))
ggsave("output/figure3_odds_ratios_environment.pdf", combined_env_coef, width = 10, height = 5.3)
ggsave("output/figure3_odds_ratios_environment.png", combined_env_coef, width = 10, height = 5.3, dpi = 300)

combined_wel_coef <- ggarrange(wel_coef$DE, wel_coef$IR, nrow = 1, widths = c(1, 1),
                               common.legend = TRUE, legend = 'bottom')
combined_wel_coef <- annotate_figure(combined_wel_coef,
  top = text_grob("Issue: Welfare", color = "black", face = "bold", size = 22))
ggsave("output/figure3_odds_ratios_welfare.pdf", combined_wel_coef, width = 10, height = 5.3)
ggsave("output/figure3_odds_ratios_welfare.png", combined_wel_coef, width = 10, height = 5.3, dpi = 300)

# --- FIGURE 4: combined predicted-probability plots ---------------------------
cat("[Figure 4] predicted distinctiveness by gov/opposition ...\n")
pred_de <- list(); pred_ir <- list()
for (ctm in names(plots)) {
  if (grepl("DE", ctm)) pred_de[[ctm]] <- plots[[ctm]]$pred_data
  else                  pred_ir[[ctm]] <- plots[[ctm]]$pred_data
}
make_combined_pred <- function(pred_list, title_) {
  d <- dplyr::bind_rows(
    pred_list[[1]] %>% dplyr::mutate(Issue = tools::toTitleCase(topic)),
    pred_list[[2]] %>% dplyr::mutate(Issue = tools::toTitleCase(topic)))
  d <- d %>% dplyr::filter(!is.na(x), !is.na(predicted)) %>% dplyr::mutate(x = droplevels(as.factor(x)))
  d$Issue <- factor(d$Issue, levels = c("Environment", "Welfare"))
  dodge_position <- position_dodge(width = 0.65)
  ggplot(d, aes(x = x, y = predicted, color = Issue, group = Issue)) +
    geom_point(position = dodge_position, size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.5, position = dodge_position, linewidth = 0.75) +
    scale_color_manual(values = c("Environment" = "darkgreen", "Welfare" = "#762A83")) +
    theme_classic() +
    labs(title = title_, y = "Predicted Values of Distinctiveness", x = "", color = NULL) +
    theme(axis.text = element_text(size = 14), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0), axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 14), legend.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = 'bold', hjust = 0.5), legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))
}
combined_germany <- make_combined_pred(pred_de, "German Green Party")
combined_ireland <- make_combined_pred(pred_ir, "Irish Green Party")
final_combined_pred <- ggarrange(combined_germany, combined_ireland, nrow = 1, widths = c(1, 1),
                                 common.legend = TRUE, legend = 'bottom')
ggsave("output/figure4_predicted_gov_opp.pdf", final_combined_pred, width = 9, height = 4.7)
ggsave("output/figure4_predicted_gov_opp.png", final_combined_pred, width = 9, height = 4.7, dpi = 300)

# --- FIGURE 5: electoral-cycle plots ------------------------------------------
cat("[Figure 5] electoral-cycle predicted distinctiveness ...\n")
env_elec <- list(); wel_elec <- list()
for (ctm in names(plots)) {
  side <- if (grepl("DE", ctm)) "DE" else "IR"
  if (grepl("environment", ctm)) env_elec[[side]] <- plots[[ctm]]$electoral
  if (grepl("welfare", ctm))     wel_elec[[side]] <- plots[[ctm]]$electoral
}
combined_env_elec <- ggarrange(env_elec$DE, env_elec$IR, nrow = 1, widths = c(1, 1))
combined_env_elec <- annotate_figure(combined_env_elec,
  top = text_grob("Issue: Environment", color = "black", face = "bold", size = 18))
ggsave("output/figure5_electoral_cycle_environment.pdf", combined_env_elec, width = 12, height = 6)
ggsave("output/figure5_electoral_cycle_environment.png", combined_env_elec, width = 12, height = 6, dpi = 300)

combined_wel_elec <- ggarrange(wel_elec$DE, wel_elec$IR, nrow = 1, widths = c(1, 1))
combined_wel_elec <- annotate_figure(combined_wel_elec,
  top = text_grob("Issue: Welfare", face = "bold", size = 18))
ggsave("output/figure5_electoral_cycle_welfare.pdf", combined_wel_elec, width = 12, height = 6)
ggsave("output/figure5_electoral_cycle_welfare.png", combined_wel_elec, width = 12, height = 6, dpi = 300)

cat("\nDone. Figures written to output/.\n")
