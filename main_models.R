# =============================================================================
# Main regression models
# "Moderation or Amplification? Mapping Speech Distinctiveness Among Green
#  Parties Using a Machine Learning Approach" (Legislative Studies Quarterly)
# Claessen, Traber & Schoonvelde
#
# Reproduces the main-text regression tables for all four country x issue
# subsets, exactly as they appear in the paper (Appendix: Tables for the Main
# Results) -- same coefficients, same standard errors, the same variable order,
# and the average marginal effects (AME) column:
#
#   Section 1  Moderation or Amplification (government/opposition)
#              4 columns:  Base | w. Controls | Complete | AME
#   Section 2  Electoral cycle
#              3 columns:  Base | Complete | AME
#
# The outcome is the SVC classifier (correct_SVC), the specification reported
# in the paper. All models are sentence-level logistic regressions with
# speaker-clustered (HC0) standard errors.
#
# Data: data/green_party_sentences_xl.rds  (plain-text export: .csv)
#   One row = one Green Party sentence classified as environmental or welfare
#   speech (Germany 1991-2005, Ireland 1992-2011). The raw sentence text is
#   not redistributed here; see README.
# =============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(sandwich, lmtest, texreg, margins)

# texreg has no built-in extract method for `margins` objects -- supply one so
# the average marginal effects print as a column of the table (as in the paper).
setMethod("extract", signature = "margins",
  definition = function(model, ...) {
    s <- summary(model)
    texreg::createTexreg(
      coef.names  = as.character(s$factor),
      coef        = s$AME,
      se          = s$SE,
      pvalues     = s$p,
      gof.names   = character(0),
      gof         = numeric(0),
      gof.decimal = logical(0))
  })

# ---- data ------------------------------------------------------------------
# The .rds and the .csv hold the same 25,531 rows and give identical results;
# read whichever you prefer. `speaker_code` is an encoding-independent integer
# id for `speaker` (one code per speaker) -- clustering on either gives the same
# numbers. It is kept because a CSV cannot store per-string encoding flags, so
# clustering on the code is the safest choice when working from the .csv:
#
#   d <- read.csv("data/green_party_sentences_xl.csv", stringsAsFactors = FALSE)
#   CLUSTER <- ~ speaker_code
#
# (Upstream note: the speaker column previously carried mixed string encodings,
# which silently split one speaker's cluster and slightly distorted the German
# cluster-robust SEs. That is fixed at source in Harmonise_Speaker_Names.R;
# see the README section "The speaker-encoding fix".)
d <- as.data.frame(readRDS("data/green_party_sentences_xl.rds"))
CLUSTER <- ~ speaker
# factor with opposition as reference category, as in the paper
d$govopp <- relevel(factor(d$govopp), ref = "Opposition")

SUBSETS <- list(
  "Germany - Environment" = d[d$country == "DE" & d$topic == "environment", ],
  "Germany - Welfare"     = d[d$country == "DE" & d$topic == "welfare", ],
  "Ireland - Environment" = d[d$country == "IR" & d$topic == "environment", ],
  "Ireland - Welfare"     = d[d$country == "IR" & d$topic == "welfare", ]
)

# --- specifications ----------------------------------------------------------
# Section 1: Moderation or Amplification (government/opposition)
F1_M1 <- correct_SVC ~ govopp
F1_M2 <- correct_SVC ~ govopp + non_green_issue_salience + non_green_distinctiveness_avg +
                       seats_perc + election + log_sentence_nchar + year
F1_M3 <- correct_SVC ~ govopp + time_until_election_neg + I(time_until_election_neg^2) +
                       non_green_issue_salience + non_green_distinctiveness_avg +
                       seats_perc + election + log_sentence_nchar + year

# Section 2: Electoral cycle
F2_M1 <- correct_SVC ~ time_until_election_neg + I(time_until_election_neg^2)
F2_M2 <- correct_SVC ~ time_until_election_neg + I(time_until_election_neg^2) + govopp +
                       non_green_issue_salience + non_green_distinctiveness_avg +
                       seats_perc + election + log_sentence_nchar + year

# --- helpers -----------------------------------------------------------------
# Logit + speaker-clustered (HC0) standard errors, as in the paper. The data are
# embedded in m$call$data so vcovCL / margins can resolve them afterwards.
fit_ct <- function(formula, data) {
  m  <- glm(formula, data = data, family = binomial(link = "logit"), x = TRUE)
  m$call$data <- data
  ct <- lmtest::coeftest(m, vcov. = sandwich::vcovCL(m, cluster = CLUSTER, type = "HC0"))
  list(m = m, ct = ct)
}
fit_ame <- function(m) {
  margins::margins(m, vcov = sandwich::vcovCL(m, cluster = CLUSTER, type = "HC0"))
}

# Coefficient labels and display order, matching the tables in the paper
# (the intercept is printed last).
S1_COEF_NAMES <- c("(Intercept)", "In Government",
                   "Issue Salience (non-Green)", "Issue Distinctiveness (non-Green)",
                   "Perc. of Seats", "Election",
                   "No. of Char. in Sentence (log)", "Year",
                   "Time Until Election", "Time Until Election^2")
S1_ORDER <- c(2, 9, 10, 3, 4, 5, 6, 7, 8, 1)

S2_COEF_NAMES <- c("(Intercept)", "Time Until Election", "Time Until Election^2",
                   "In Government",
                   "Issue Salience (non-Green)", "Issue Distinctiveness (non-Green)",
                   "Perc. of Seats", "Election",
                   "No. of Char. in Sentence (log)", "Year")
S2_ORDER <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 1)

for (nm in names(SUBSETS)) {
  df <- SUBSETS[[nm]]
  cat("\n\n=====================================================================\n")
  cat(nm, " (N =", nrow(df), ")\n")
  cat("=====================================================================\n")

  s1m1 <- fit_ct(F1_M1, df)
  s1m2 <- fit_ct(F1_M2, df)
  s1m3 <- fit_ct(F1_M3, df)
  s1ame <- fit_ame(s1m3$m)

  s2m1 <- fit_ct(F2_M1, df)
  s2m2 <- fit_ct(F2_M2, df)
  s2ame <- fit_ame(s2m2$m)

  cat("\n----", nm, ": Moderation or Amplification (Government/Opposition)\n")
  cat(screenreg(list(s1m1$ct, s1m2$ct, s1m3$ct, s1ame),
        custom.model.names = c("Base", "w. Controls", "Complete", "AME"),
        custom.coef.names  = S1_COEF_NAMES,
        reorder.coef       = S1_ORDER,
        stars = c(.001, .01, .05, .1), symbol = "+", digits = 2,
        custom.gof.rows = list(
          "N"   = c(nobs(s1m1$m), nobs(s1m2$m), nobs(s1m3$m), nobs(s1m3$m)),
          "AIC" = c(round(AIC(s1m1$m), 1), round(AIC(s1m2$m), 1),
                    round(AIC(s1m3$m), 1), round(AIC(s1m3$m), 1))),
        include.deviance = FALSE, include.loglik = FALSE, include.bic = FALSE,
        include.nobs = FALSE, include.aic = FALSE), "\n")

  cat("\n----", nm, ": Electoral Cycle\n")
  cat(screenreg(list(s2m1$ct, s2m2$ct, s2ame),
        custom.model.names = c("Base", "Complete", "AME"),
        custom.coef.names  = S2_COEF_NAMES,
        reorder.coef       = S2_ORDER,
        stars = c(.001, .01, .05, .1), symbol = "+", digits = 2,
        custom.gof.rows = list(
          "N"   = c(nobs(s2m1$m), nobs(s2m2$m), nobs(s2m2$m)),
          "AIC" = c(round(AIC(s2m1$m), 1), round(AIC(s2m2$m), 1), round(AIC(s2m2$m), 1))),
        include.deviance = FALSE, include.loglik = FALSE, include.bic = FALSE,
        include.nobs = FALSE, include.aic = FALSE), "\n")
}

cat("\nDone. Coefficients, standard errors, significance stars, variable order\n")
cat("and AMEs correspond to the tables in the paper (Appendix: Tables for the\n")
cat("Main Results).\n")
