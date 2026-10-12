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
| `3_Training_Soil_Data_2015_2023.csv` | Soil properties | 186 envs, 68.4% overlap — weakest file |
| `4_Training_Weather_Data_2014_2023_seasons_only.csv` | Weather by growth stage | 269 envs, 98.9% overlap |
| `6_Training_EC_Data_2014_2023.csv` | Derived environmental covariates, 655 cols | 241 envs, 88.6% overlap |
| `5_Genotype_Data_All_2014_2025_Hybrids.vcf` | Raw VCF | Not re-parsed — use numerical file |
| `5_Genotype_Data_All_2014_2025_Hybrids_numerical.txt` | SNP dosage matrix | Preferred for ML |

**Key columns in trait file:**
- Environment ID: `Env` (e.g. `DEH1_2014`)
- Hybrid ID: `Hybrid`
- Target variable: `Yield_Mg_ha` (94.8% non-missing)
- Flowering traits: `Pollen_DAP_days`, `Silk_DAP_days`
- All files join on `Env`; genotype file joins on `Hybrid`

---

## GWAS Project — Complete ✅

**Scripts 00–03** are complete and committed. Do not modify them.

**Key findings:**
- S5_11998173 — plantacyanin gene, 3' UTR variant (REF=T, ALT=C), −4.1 days flowering, MAF 1.07%. Flagship finding.
- S8_113222492 — Vgt1 region (Chr8), known flowering locus. Validates the pipeline.
- Chr3 @ 156 Mb — strongest statistical signal but near-zero effect; likely LD proxy, needs fine-mapping.
- 2,060 SNPs after QC; all hybrids show only genotypes 0/1 (F1 hybrid biology, expected)
- FarmCPU preferred over GLM — corrects for population structure and kinship

**Processed data available:**
- `data/processed/flowering_BLUEs_for_GWAS.csv` — 2,420 hybrids × 2 cols
- `data/processed/genotypes_for_GWAS.csv` — 2,420 hybrids × 2,061 cols
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
- [x] Exercise 1: Orient in trait data (173,960 rows, 94.8% yield coverage, 2014–2023)
- [x] Exercise 2: Check auxiliary file coverage (table above)
- [ ] Exercise 2 remaining: Identify useful ML features from meta file
- [ ] Exercise 3: Yield distribution and coverage deep-dive
- [ ] Exercise 4: Load and validate numerical genotype file
- [ ] Exercise 5: Join all files, check row survival, quantify missingness
- [ ] Exercise 6: Final cleaned dataset → `data/processed/`

**Modelling plan:**
- Baseline: linear regression
- Tree-based: Random Forest, XGBoost
- Genomic benchmark: GBLUP (rrBLUP) — industry standard to beat
- CV strategy: 5-fold stratified by environment (not random — avoids data leakage)
- Metrics: Pearson r, R², RMSE

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

- **F1 hybrid genotypes:** Only 0 (REF/REF) and 1 (REF/ALT) observed — genotype 2 absent by biology, not error
- **GAPIT effect sizes:** Reported assuming 0/1/2 dosage — always validate against observed genotype-stratified means
- **BLUEs vs BLUPs:** BLUEs used for GWAS (fixed effects, no shrinkage); BLUPs would attenuate allele effects
- **Soil data starts 2015:** 2014 rows will have missing soil features — needs imputation strategy
- **Meta column names are messy:** Long strings with special characters; will need renaming before use
- **EC file has 655 columns:** Derived crop model outputs by growth stage — rich but needs selection/reduction
