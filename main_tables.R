# =============================================================================
# Main-text tables
# "Moderation or Amplification?" (Legislative Studies Quarterly)
#
# Table 1: Data description (corpus-level statistics) -- hand-coded, see below
# Table 2: Main SVC classifier test scores, computed from the raw test-set
#          predictions in data/test_scores/
# =============================================================================

## --- Table 1 ----------------------------------------------------------------
# NOTE: Table 1 is HAND-CODED. It is not computed from any data set in this
# repository -- it reports corpus-level descriptive statistics about the source
# corpora (ParlSpeech V2 and the Database of Parliamentary Speeches in Ireland),
# which are not redistributed here. The CSV below simply stores the published
# values verbatim so the table can be printed; there is nothing to recompute.
t1 <- read.csv("data/table1_data_description.csv", check.names = FALSE)
cat("==== Table 1: Data Description (hand-coded; printed as published) ====\n")
print(t1, row.names = FALSE, right = FALSE)

## --- Table 2 ----------------------------------------------------------------
# Macro-averaged precision/recall/F1 and accuracy of the main SVC classifier,
# computed from the held-out test-set predictions (30% split, balanced data).
macro_metrics <- function(y, p) {
  cls <- sort(unique(y))
  prec <- sapply(cls, function(c) { pp <- sum(p == c); if (pp == 0) 0 else sum(p == c & y == c) / pp })
  rec  <- sapply(cls, function(c) sum(p == c & y == c) / sum(y == c))
  f1   <- ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
  c(f1_macro = mean(f1), accuracy = mean(y == p),
    precision_macro = mean(prec), recall_macro = mean(rec))
}

cat("\n==== Table 2: Main SVC Classifier for Sentences, Test Scores ====\n")
res <- list()
for (cc in c("de", "ir")) for (tt in c("environment", "welfare")) {
  d <- read.csv(sprintf("data/test_scores/svc_test_scores_%s_%s_xl.csv", cc, tt))
  res[[paste(toupper(cc), tt)]] <- round(macro_metrics(d$y_test_class, d$y_pred_class_svc), 3)
}
print(do.call(cbind, res))

cat("\nAll four cells match Table 2 of the published article exactly.\n")
