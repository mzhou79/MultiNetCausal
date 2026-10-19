# uConnect missing-outcome comparison

## What this script does

`missingness_check.R` independently builds a descriptive table comparing participants with observed versus missing Wave 3 HIV-risk-behavior outcomes within Mo's baseline/treatment-complete cohort. It does not run or modify `uconnect_prep_mo.R`, fit models, impute data, or test/prove MCAR. The local artifact contains code only; use the real inputs and retain its outputs inside Stronghold.

## Run inside Stronghold

1. Put this reviewed script in your own Stronghold working folder through the permitted code-import route. Do not overwrite Mo's source.
2. Review the input paths at the top. They refer to the node/edge and Wave 1/Wave 3 files seen in the current prep source.
3. In RStudio, source this script and explicitly run:

```r
run_missingness_check()
```

Sourcing alone loads definitions and does not read data or generate output. No additional R packages are required. Each successful run writes a new timestamped folder below `getwd()/uconnect_missingness_check/`:

- `table.csv`: treatment, seven baseline covariates, Wave 3 outcome availability/observed prevalence, sexual and augmented degree, and isolate summaries.
- `cohort_flow.csv`: actual sample counts, exclusions, Wave 3 linkage/missing-field counts, network counts, and comparison with manuscript-reported reference counts.
- `source_reconciliation.txt`: exact input paths, selected Wave 3 field, definitions, denominator/SMD conventions, and limitations.

No participant/ID/edge-level files are saved. All aggregate files remain in Stronghold pending its normal disclosure process.

## Important reconciliation with the prep source and manuscript

The reference source was visually read through line 426 in Stronghold, without executing, editing, copying through the clipboard, or transferring the file. The formal current `main_manuscript.tex` was read independently by the literature specialist, and its application section was selectively checked by main.

- The manuscript studies **baseline PrEP awareness** and **Wave 3 HIV risk behavior**, not HIV seroconversion or PrEP use. Its stated cohort flow is 618 enrolled -> 586 baseline/treatment-complete -> 468 observed outcomes and 118 missing outcomes. These are reference values, not independently reproduced counts or enforced targets.
- The executable prep excludes missing healthcare, notwithstanding old comments saying healthcare missingness is retained as a category. This script follows that executable baseline mask.
- Employment uses original codes 1/2 -> employed and 3 -> not employed. All seven baseline covariates and treatment must be nonmissing. There is no baseline-HIV-negative restriction.
- The displayed prep has a Wave 3 lookup inconsistency: its Wave 3 mapping uses the `_w1` field name, and its actual lookup reads `risk_w1` rather than the loaded `risk_w3`. The comparison script instead joins the actual Wave 3 file by `su_id`; it never trusts the original `obs_y` or sources the original prep.
- The Wave 3 file's full schema was not opened. `W3_OUTCOME_COLUMN <- NULL` accepts only an unambiguous supplied v4 risk field (`hivriskbehav_v4_w3`, or a legacy `_w1`-suffixed field inside the actual Wave 3 file). The exact binding is printed and recorded. If ambiguous or absent, execution stops for explicit column selection. A legacy suffix does not itself verify the field's wave/content; confirm the recorded binding with Mo's variable construction before interpreting the real-data table.
- Mo's derived risk fields are used as supplied. Component construction, partial-response rules, HIV-status eligibility and structural skips are not independently verified. The manuscript itself uses both “missing” and “not well defined or missing.” A missing derived outcome must not automatically be called dropout/nonresponse.
- Degrees use the same full-cohort undirected simple network for both missingness groups. Sexual ties are type 2, augmented ties union types 1/2. Between-group edges count; isolates remain. This is not degree recalculated in separate group-induced graphs.
- Absolute SMDs are descriptive summaries, not network-adjusted tests. Similarity cannot establish MCAR, and this comparison does not investigate the initial baseline/treatment exclusions.

## Evidence boundary

Code preparation and synthetic checks are distinct from execution with protected data. No real-data missingness table or manuscript count reproduction has been completed by preparing this artifact. Read the real run's `cohort_flow.csv` and selected Wave 3 binding before treating its comparison as application evidence. Deviations are reported rather than hidden or forced away.
