# CLAUDE.md — g2f-maize-quantitative-genetics

## Project Overview
This is a learning and portfolio project by Nikolaos, a plant biology PhD candidate at Ghent University transitioning to computational biology roles in ag-biotech. The project performs GWAS and ML-based yield prediction on the G2F (Genomes to Fields) maize dataset.

The codebase will eventually be containerised as a Nextflow pipeline.

---

## User Preferences
- **Coding style:** Exercise-driven learning. Prefer guided tasks over complete solutions. Give hints, not answers, unless explicitly asked.
- **Language:** R (strong), Python (learning)
- **Style:** Tidyverse conventions throughout. Pipes (`|>`) preferred over nested calls.
- **Output:** Always print informative messages with `cat()` or `sprintf()`. Save figures to the appropriate `figures/` subdirectory.

---

## Repository Structure

```
g2f-maize-quantitative-genetics/
├── 01-gwas-flowering-time/
│   ├── scripts/                  # Scripts 00-03, complete
│   ├── figures/                  # Committed figures
│   └── results/
│       └── gwas-output/          # Local only — gitignored (large GAPIT files)
├── 02-yield-prediction/
│   ├── scripts/                  # In progress
│   ├── figures/
│   └── results/
├── data/
│   ├── raw/                      # Local only — gitignored
│   └── processed/                # Committed outputs
└── CLAUDE.md
```

---

## Data Files (all in `data/raw/`)

| File | Contents | Notes |
|---|---|---|
| `1_Training_Trait_Data_2014_2023.csv` | All phenotypes, 173,960 rows × 26 cols | Main trait file |
| `2_Training_Meta_Data_2014_2023.csv` | Environment metadata | 272 envs, 100% overlap with trait |
| `3_Training_Soil_Data_2015_2023.csv` | Soil properties, 36 cols | 186 envs, 68.4% of envs / 74.1% of yield plots — weakest file |
| `4_Training_Weather_Data_2014_2023_seasons_only.csv` | **Daily** NASA POWER weather over the growing season (`Env`, `Date`, 16 vars) — up to 241 rows/env, must be aggregated before joining | 269 envs, 98.9% overlap |
| `4_Training_Weather_Data_2014_2023_full_year.csv` | Same, full calendar year | Not yet used |
| `6_Training_EC_Data_2014_2023.csv` | Derived environmental covariates, 655 cols | 241 envs, 88.6% overlap |
| `key_inbreds_G2F_2014-2025.txt` | Inbred parent key | Not yet explored |
| `5_Genotype_Data_All_2014_2025_Hybrids.vcf` | Raw VCF | Not re-parsed — use numerical file |
| `5_Genotype_Data_All_2014_2025_Hybrids_numerical.txt` | SNP dosage matrix | Preferred for ML |

**Key columns in trait file:**
- Environment ID: `Env` (e.g. `DEH1_2014`)
- Hybrid ID: `Hybrid`
- Target variable: `Yield_Mg_ha` (94.8% non-missing)
- Flowering traits: `Pollen_DAP_days`, `Silk_DAP_days`
- All files join on `Env`; genotype file joins on `Hybrid`

---

## GWAS Project — Re-run in progress ⚠️

**Scripts 00–01** are complete — do not modify. **Script 02 (and re-run of 03)** is being fixed, see below.

**Genotype coding bug (found 2026-10):** `02-process-genotypes.R` used `extract.gt(..., as.numeric = TRUE)`.
vcfR converts GT strings with C `atof()`, which stops at the `/`, so the value is the *first allele only*
(`"0/1"` → 0, `"1/0"` → 1, `"1/1"` → 1, `"2/2"` → 2). Consequences in `genotypes_for_GWAS.csv`:
hybrids sharing the first-named parent are identical across all SNPs (964/964 groups; only 1,139 unique
profiles among 2,420 hybrids), tester alleles are invisible, and the "only 0/1" pattern is an artefact.
MAF was also filtered on all 5,899 hybrids before subsetting to the 2,420 analysed.
Findings below are **provisional** until the fix exercises are done.

**Fix exercises (user runs locally — needs raw VCF + R):**
- [ ] G1: Inspect raw GT strings (`table()`), check phasing / allele order
- [ ] G2: Drop non-biallelic sites; convert GT strings to ALT-allele dosage 0/1/2 (handle `./.`)
- [ ] G3: Subset to analysis hybrids *before* MAF/missingness filtering
- [ ] G4: Validate — genotype 2 present; hybrids sharing a parent no longer identical
- [ ] G5: Re-run `03-run-gwas.R`; compare old vs new top hits; re-check S5_11998173 by genotype-stratified BLUE means
- [ ] G6: Update `GWAS_summary.md` and the key findings / technical notes here

**Key findings (provisional — from first-allele-only genotypes):**
- S5_11998173 — plantacyanin gene, 3' UTR variant (REF=T, ALT=C), −4.1 days flowering, MAF 1.07%. Flagship finding.
- S8_113222492 — Vgt1 region (Chr8), known flowering locus. Validates the pipeline.
- Chr3 @ 156 Mb — strongest statistical signal but near-zero effect; likely LD proxy, needs fine-mapping.
- 2,060 SNPs after QC (filters to be re-applied after the coding fix)
- FarmCPU preferred over GLM — corrects for population structure and kinship

**Processed data available:**
- `data/processed/flowering_BLUEs_for_GWAS.csv` — 2,420 hybrids × 2 cols
- `data/processed/genotypes_for_GWAS.csv` — 2,420 hybrids × 2,061 cols — **wrong coding, do not use for ML**
- `data/processed/SNP_map_for_GWAS.csv` — 2,060 SNPs × 3 cols

---

## ML Yield Prediction Project — In Progress 🔄

**Goal:** Predict `Yield_Mg_ha` using progressively richer feature sets.

**Three phases:**

| Phase | Features | Question |
|---|---|---|
| A | Environment + soil + weather | What does growing context alone explain? |
| B | + Genotype (SNPs) | How much does genetics add? |
| C | + Flowering time | Does phenology improve prediction? |

**Todo (current position: mid-Exercise 2):**
- [x] Exercise 1: Orient in trait data (173,960 rows, 94.8% yield coverage, 2014–2023) — `scripts/01-explore-trait-data.R`
- [x] Exercise 2: Check auxiliary file coverage (table above) — `scripts/02-check-auxiliary-coverage.R`
  - Rebuilt 2026-10; all `[OK]` checks pass on real data. 0 envs lack yield entirely; 45 locations, 5,205 hybrids
- [ ] Exercise 2 remaining: Meta feature audit → `results/meta_feature_audit.csv` (column, decision, reason) — apply the legitimacy filter below
- [ ] Exercise 3: Yield distribution and coverage deep-dive
- [ ] Exercise 3 also: variance decomposition `Yield ~ (1|Env) + (1|Hybrid) + (1|Env:Hybrid)` — sets the ceiling per phase; confirm 2-year hybrid/tester blocks across all years
- [ ] Exercise 4: Load and validate numerical genotype file (must contain 0/1/2; hybrids sharing a parent must differ)
- [ ] Exercise 5: Join all files, check row survival, quantify missingness
- [ ] Exercise 6: Final cleaned dataset → `data/processed/`

**Modelling plan:**
- Baselines: training mean → environment means → ridge/lasso
- Tree-based: Random Forest, XGBoost
- Genomic benchmark: GBLUP (rrBLUP) — industry standard to beat
- GxE: SNP PCs × env PCs (trees) or reaction-norm model
- Language split: R for data prep (tidyverse), Python for ML (sklearn), GBLUP in R (rrBLUP/sommer)
- Metrics: Pearson r, R², RMSE — reported **overall and within-environment** (mean of per-env r), since overall r is dominated by environment means

**CV design (decided):**
- *Story:* serve both a seed company forecasting next season and a breeder ranking new lines.
- *Data:* all years 2014–2023; unit = Env × Hybrid mean yield (replicate plots averaged).
- *Loop 1 — leave-one-year-out (10 folds):* headline for Phase A. Tag each test prediction as
  known/new hybrid × known/new location (vs. training years) and report metrics per cell.
  Note: hybrids run in 2-year tester blocks (e.g. 2018–19 LH195/PHT69, 2020–21 PHZ51/PHP02/PHK76),
  so most LOYO test hybrids were seen in the paired year.
- *Loop 2 — new lines in known environments:* GroupKFold grouped by **line** = first-named parent in
  `Hybrid` (testers are second-named). Headline for Phase B. Phase C re-runs both loops.
- *Weather:* hindcast — actual season weather/EC used as if known. State this in every write-up.
- *Leakage rules:* hyperparameter tuning uses the same grouping inside training folds; anything fitted
  (PCA, scaling, imputation, flowering BLUEs) is fitted on training folds only.
- *Fold assignments* saved once to `data/processed/` and shared by R and Python.
- *Excluded features (post-hoc measurements):* `Date_Harvested`, `Grain_Moisture`, `Twt_kg_m3`, lodging,
  plant/ear height. `Stand_Count_plants` borderline — justify if used.

**Feature selection (decided) — two separate questions:**
1. *Is a feature allowed?* — answered by domain rules in the audits (Exercise 2, and later soil/weather/EC).
   Drop **only** for legitimacy: leakage (recorded after the fact, e.g. damage comments, harvest dates),
   IDs in disguise (farm/field/trial codes — the model memorises sites instead of learning site properties),
   empty/constant columns, exact duplicates. **Default is keep** — "maybe useless" is for the model to decide.
   Rationale: leakage and ID columns *win* performance contests, so score-based selection can't police them.
2. *Is a feature useful?* — answered empirically, **always inside the CV loop** (selecting on all data, then
   CV-ing, gives optimistic scores):
   - Phase A selection runs on **environment-mean yield (272 rows)**, since env features only explain env-level
     variance — fast enough for thousands of fits. The chosen env set is then carried into row-level Phase B/C models.
   - Methods: lasso/elastic-net path (cheap first pass); block ablation (drop meta / soil / weather / EC as
     whole groups); recursive feature elimination by permutation importance (drop weakest 10–20% per round).
     Not one-at-a-time backward elimination — ~p² fits and it overfits the CV.
   - Pick with the **1-SE rule**: simplest model within one standard error of the best.
   - Report **selection stability**: how many of the 10 folds selected each feature. Expect correlated features
     (e.g. adjacent-month temperatures) to trade places — interpret groups, not single winners.
3. *Portfolio leakage demo:* deliberately add one leakage column and one farm ID; show the score rises under
   random CV and in the LOYO known-location cell, then collapses in the new-location cell.

---

## Coding Conventions

- Use `tidyverse` throughout; load at top of every script
- Always set paths relative to project root
- Raw data lives in `data/raw/` — never modify these files
- Processed outputs go to `data/processed/`
- Figures go to the relevant `figures/` subdirectory
- Use `ggsave()` for all plots with explicit `width`, `height`, `dpi`
- Comment the *why*, not the *what*

---

## Git Rules
- Raw data and GAPIT output are gitignored — never commit them
- Large files (>50MB) must never be staged — check `git status` before committing
- If a large file is accidentally committed: `git filter-branch` to scrub history, then force push
- Kinship matrix (`GAPIT.Genotype.Kin_Zhang.csv`, ~100MB) has already caused one incident — it is gitignored

---

## Key Technical Notes

- **Hybrid genotypes:** Correct dosage is 0/1/2 — genotype 2 occurs whenever both parents carry ALT. The earlier "only 0/1" observation was the vcfR `as.numeric` bug (see GWAS section)
- **GAPIT effect sizes:** Reported assuming 0/1/2 dosage — always validate against observed genotype-stratified means
- **BLUEs vs BLUPs:** BLUEs used for GWAS (fixed effects, no shrinkage); BLUPs would attenuate allele effects
- **Soil gaps:** 86 envs lack soil — all of 2014, but also 2018 (only 38% covered) and 2015 (56%). Filling from the same `Field_Location` in another year recovers 76 of 86; 10 stay missing → also add a missing-indicator flag
- **Meta column names are messy:** Long strings with special characters; will need renaming before use
- **EC file has 655 columns:** Derived crop model outputs by growth stage — rich but needs selection/reduction
