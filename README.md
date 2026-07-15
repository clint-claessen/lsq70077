# Replication materials — Moderation or Amplification?

**Moderation or Amplification? Mapping Speech Distinctiveness Among Green Parties Using a Machine Learning Approach**
Clint Claessen, Denise Traber, Martijn Schoonvelde. *Legislative Studies Quarterly*.

This repository reproduces **everything in the main paper**: Tables 1–2,
Figures 1–5, and the main regression models (government/opposition and
electoral-cycle specifications) for Green Party speech distinctiveness in the
German Bundestag and the Irish Dáil.

## Contents

| File | Reproduces |
|---|---|
| `main_models.R` | The main-text regression tables for all 4 country × issue subsets — both sections, speaker-clustered SEs, AME columns, printed in the paper's variable order |
| `main_figures.R` | Figures 1–5, written to `output/` — exact reproductions of the published figures (see *Figures* below) |
| `main_tables.R` | Table 1 (hand-coded corpus statistics, printed as published — nothing to recompute) and Table 2 (SVC test scores, computed from the raw test-set predictions) |
| `data/green_party_sentences_xl.rds` | Sentence-level input: 25,531 Green Party sentences with all model variables. Figure 2 aggregates it; Figures 3–5 and all regression models refit from it. Read by default — reproduces the published tables exactly. |
| `data/green_party_sentences_xl.csv` | The same 25,531 rows as plain text, for inspection or use outside R. For CSV-based reproduction cluster on `speaker_code` (see *The speaker-encoding fix* below). |
| `data/salience_sentences_xl.csv` | **Sentence-level salience coding** for Figure 1: one row per Green environment/welfare sentence (`country`, `date`, `topic`), full corpus date range, no text. The Figure 1 numerator is counted from this. |
| `data/salience_monthly_totals_xl.csv` | Total Green sentences per month across all topics — the Figure 1 denominator (`country`, `year`, `month`, `total_sentences`) |
| `data/figure1_salience_monthly.csv` | Monthly Green Party issue salience: the aggregation of the two files above, i.e. the exact values plotted in Figure 1 |
| `data/figure2_distinctiveness_bars.csv` | Mean distinctiveness by government/opposition: the aggregation of `green_party_sentences_xl` computed by `main_figures.R`, i.e. the exact values plotted in Figure 2 |
| `data/table1_data_description.csv` | The hand-coded corpus statistics printed as Table 1 (stored verbatim; they describe the source corpora, which are not redistributed here) |
| `data/test_scores/svc_test_scores_*.csv` | Held-out test-set predictions of the main SVC classifier (Table 2) |
| `output/` | Generated figures (PNG + PDF) |

## How to run

```r
# from the repository root
source("main_tables.R")    # Tables 1 and 2
source("main_figures.R")   # Figures 1-5 -> output/
source("main_models.R")    # main regression models (console tables)
```

Dependencies: `ggplot2`, `dplyr`, `ggpubr`, `sandwich`, `lmtest`, `sjPlot`,
`marginaleffects`, `forcats`, `scales`, `texreg`, `margins`. Tested on R 4.5.2
and R 4.6.0.

`main_models.R` and `main_tables.R` have no unusual requirements. The Figure 3
coefficient plots in `main_figures.R` use `sjPlot::plot_models`, which needs
**`insight` ≥ 1.5.0.6**; the script checks this up front. If you hit

```
! Namespace 'insight' 1.3.1 is already loaded, but >= 1.5.0.6 is required
```

then an older `insight` was loaded earlier in a long-running session and cannot
be upgraded in place — **restart R** and re-run (and `install.packages("insight")`
first if the installed version is itself older than 1.5.0.6).

## Figures

`main_figures.R` is not a look-alike reconstruction: the plotting code is the
verbatim code from the analysis pipeline, run against the bundled data, so the
output matches the article. Both a PDF and a 300-dpi PNG are written to
`output/` for each figure.

Verified pixel-by-pixel against the published PDFs. The published figures were
themselves produced under two R versions — Figures 1–2 under R 4.5.2 and
Figures 3–5 under R 4.6.0 — so which subset comes out *byte*-identical depends on
the R you run:

| Run under | Figures 1–2 | Figures 3–5 |
|---|---|---|
| R 4.5.2 (ggplot2 3.5.2) | **pixel-identical** | visually identical |
| R 4.6.0 (ggplot2 4.0.3) | visually identical | **pixel-identical** |

"Visually identical" means the only difference is sub-pixel text anti-aliasing
from the ggplot2 version: the curves, bars, ribbons, points, intervals, colours
and legends are unchanged (mean difference < 0.4/255, confined to glyph edges).
Either R reproduces the figures in the article.

Every figure is built from text-free, sentence-level data, so each is auditable
down to individual sentences:

- **Figure 1** — `main_figures.R` counts, per month, the Green environment /
  welfare sentences in `data/salience_sentences_xl.csv` (numerator) and divides
  by the total Green sentences that month in `data/salience_monthly_totals_xl.csv`
  (denominator): `percentage = 100 × sum_topic / sum_total`.
- **Figure 2** — each bar is the mean of `correct_SVC` over the Green sentences of
  that issue while the party was in government vs. opposition (government/opposition
  assigned from the sitting date), aggregated from `data/green_party_sentences_xl.rds`.
- **Figures 3–5** — the logit models are refit from `data/green_party_sentences_xl.rds`
  and plotted with the pipeline's `sjPlot` / `marginaleffects` code.

Each aggregation reproduces the plotted values exactly
(`data/figure1_salience_monthly.csv`, `data/figure2_distinctiveness_bars.csv`).

## Model specification

The outcome throughout is `correct_SVC`, the SVC classifier reported in the
paper. All regressions are logistic regressions at the sentence level with
speaker-clustered (HC0) standard errors:

```
# Section 1 — Moderation or Amplification (government/opposition)
correct_SVC ~ govopp                                                   # Base
correct_SVC ~ govopp + <controls>                                      # w. Controls
correct_SVC ~ govopp + time_until_election_neg +
              I(time_until_election_neg^2) + <controls>                # Complete

# Section 2 — Electoral cycle
correct_SVC ~ time_until_election_neg + I(time_until_election_neg^2)   # Base
correct_SVC ~ time_until_election_neg + I(time_until_election_neg^2) +
              govopp + <controls>                                      # Complete

# <controls> = non_green_issue_salience + non_green_distinctiveness_avg +
#              seats_perc + election + log_sentence_nchar + year
```

`main_models.R` prints each table with the coefficients, speaker-clustered SEs,
significance stars, variable order and average marginal effects (AME) column of
the corresponding table in the paper. Figures 4–5 show predicted probabilities at
means/modes computed by the delta method on the clustered covariance matrix.

### The speaker-encoding fix

The `.rds` and the `.csv` hold the same 25,531 rows and give identical results;
read whichever you prefer, and cluster on `speaker` or on `speaker_code` — all
four combinations agree.

That is worth a note, because it was not always true. `speaker` is the clustering
variable for every standard error in the paper, and in the source data the
speaker names used to carry **mixed string encodings**: 28 of the 303 sentences by
"Katrin Göring-Eckardt" were flagged UTF-8 and the other 275 native (she is the
one Green affected). R's `identical()`, `match()` and `table()` all translate
encodings and reported them as equal, so the split was invisible — but `rowsum()`,
which `sandwich::vcovCL()` uses to build the clustered meat matrix, matches the
raw string and did not, silently leaving those 28 rows out of her cluster.

This is now **fixed at source**: the pipeline's speaker-harmonisation step
normalises the whole column with `enc2utf8()`, so every row of a speaker clusters
together. The tables in the published article were regenerated accordingly.
Coefficients, N and AIC were never affected (`speaker` enters only the variance),
Ireland was never affected, and no figure value or significance star changed.

**`speaker_code`** is an encoding-independent integer id (one code per speaker).
It is kept because a CSV cannot store per-string encoding flags, so clustering on
the code is the most robust choice when working from the `.csv`:

```r
d <- read.csv("data/green_party_sentences_xl.csv", stringsAsFactors = FALSE)
# cluster = ~ speaker_code   # identical results to ~ speaker
```

## Data description (`green_party_sentences_xl.csv`)

One row is one sentence spoken by the German Green Party (Bundestag, 1991–2005)
or the Irish Green Party (Dáil, 1992–2011) classified as environmental or
welfare speech by the guided-BERTopic pipeline described in the paper.

| Variable | Description |
|---|---|
| `sentence_id` | Sentence identifier within the source speech corpus |
| `country` | `DE` (Germany) / `IR` (Ireland) |
| `topic` | `environment` / `welfare` |
| `date`, `year`, `year_month` | Sitting date of the parent speech |
| `speaker` | Speaker of the parent speech (the clustering variable for the SEs, as published) |
| `speaker_code` | Encoding-independent integer id for `speaker` (one code per speaker). Cluster on this when working from the `.csv` — see *The speaker-encoding fix* above. |
| `govopp` | Government vs. opposition status of the Green Party |
| `election` | Electoral period indicator |
| `seats_perc` | Green Party seat share (%) |
| `time_until_election`, `time_until_election_neg` | Years until the next election (positive / negative coding; the models use the negative coding, −5 to 0) |
| `non_green_issue_salience` | Monthly issue salience of all non-Green parties (%) |
| `non_green_distinctiveness_avg` | Monthly average speech distinctiveness of all non-Green parties |
| `sentence_nchar`, `log_sentence_nchar` | (Log) number of characters in the sentence |
| `terms` | Number of terms in the parent speech |
| `correct_SVC` | Outcome: 1 if the SVC classifier attributed the sentence to the Green Party, 0 otherwise |

**Note on text.** The raw sentence text is not redistributed. The underlying
speeches are available from the ParlSpeech V2 dataset (Rauh & Schwalbach 2020)
and the Database of Parliamentary Speeches in Ireland (Herzog & Mikhaylov
2017); the full processing pipeline (topic modeling, classification, data
preparation) is available from the authors upon reasonable request.

## Citation

> Claessen, C., Traber, D., & Schoonvelde, M. Moderation or Amplification?
> Mapping Speech Distinctiveness Among Green Parties Using a Machine Learning
> Approach. *Legislative Studies Quarterly*.

## Contact

Clint Claessen — clint.claessen@plus.ac.at
