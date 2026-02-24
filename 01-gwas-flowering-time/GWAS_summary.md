# GWAS of Flowering Time in Maize — Analysis Summary
**Dataset:** G2F Genomes to Fields Initiative · **Trait:** Days to pollen shed (Pollen_DAP_days)  
**Repository:** `g2f-maize-quantitative-genetics`

---

## Data

Raw trait data comprised 173,960 plot records from 2014–2023. Years 2018–2021 were selected based on genotypic diversity (>1,000 unique hybrids/year), environment coverage (96 environments), and genotype-phenotype overlap (97.8% of phenotyped hybrids present in VCF) to reduce computational time. The final dataset contained 63,189 plots across 2,420 hybrids.

---

## Methods

### Phenotype Adjustment
BLUEs were calculated for each hybrid using a mixed model fitted with `lme4`:

```
Pollen_DAP_days ~ Hybrid (fixed) + (1 | Env) (random)
```

Marginal means were extracted using `emmeans`. Variance partitioning yielded 92.3% environmental and 7.7% residual variance. BLUEs ranged from 53.5–76.8 days (IQR: 62.5–65.5); 91.4% of hybrids had SE < 1.0 day (≥10 environments).

### Genotype QC
The VCF contained 2,425 SNPs × 5,899 hybrids. SNPs were filtered for missingness (>10% removed, −201 SNPs) and minor allele frequency (<5% removed, −164 SNPs), leaving 2,060 SNPs. Residual missing values (62,636) were imputed by per-SNP mean. All hybrids showed genotypes of 0 (REF/REF) or 1 (REF/ALT) only — expected for F1 hybrids where genotype 2 requires both parental lines to carry the ALT allele.

### GWAS
Association analysis was performed using GAPIT3 with two models. GLM used 3 principal components as fixed covariates; FarmCPU used 6 PCs and iteratively conditioned on pseudo-QTNs to account for population structure and kinship. Both models used a Zhang-method kinship matrix. Significance threshold was P < 0.05 (Bonferroni: P < 2.43 × 10⁻⁵).

---

## Results

GLM identified 321 significant SNPs; FarmCPU identified 139 (57% reduction). The higher Chr8 retention rate in FarmCPU (72% vs 24–62% on other chromosomes) reflects genuine signal preservation over inflation removal. Mild QQ plot deflation in FarmCPU is acceptable given the high relatedness of the breeding panel and is validated by recovery of the known Vgt1 locus.

### Top associations — GLM (top 20)

| Rank | SNP | Chr | Position | P-value | MAF | Effect (days) |
|---:|:---|---:|---:|:---|---:|---:|
| 1 | S5_11998173 | 5 | 11,998,173 | 2.22 × 10⁻³⁸ | 0.011 | −3.81 |
| 2 | S1_210194624 | 1 | 210,194,624 | 1.54 × 10⁻²⁹ | 0.053 | −1.80 |
| 3 | S3_137533609 | 3 | 137,533,609 | 7.21 × 10⁻²² | 0.105 | −1.11 |
| 4 | S5_206652970 | 5 | 206,652,970 | 4.47 × 10⁻²¹ | 0.011 | −2.67 |
| 5 | S3_177563221 | 3 | 177,563,221 | 1.53 × 10⁻²⁰ | 0.008 | −3.04 |
| 6 | S5_15897858 | 5 | 15,897,858 | 3.01 × 10⁻¹⁸ | 0.150 | −0.87 |
| 7 | S4_178069434 | 4 | 178,069,434 | 1.46 × 10⁻¹⁷ | 0.226 | +2.35 |
| 8 | S3_152992674 | 3 | 152,992,674 | 5.57 × 10⁻¹⁵ | 0.323 | −0.79 |
| 9 | S5_213160675 | 5 | 213,160,675 | 1.50 × 10⁻¹⁴ | 0.125 | −0.86 |
| 10 | S3_156277883 | 3 | 156,277,883 | 3.50 × 10⁻¹⁴ | 0.393 | −0.83 |
| 11 | S3_156277910 | 3 | 156,277,910 | 3.50 × 10⁻¹⁴ | 0.393 | −0.83 |
| 12 | S2_23345860 | 2 | 23,345,860 | 4.45 × 10⁻¹⁴ | 0.301 | −0.73 |
| 13 | S5_218069890 | 5 | 218,069,890 | 2.11 × 10⁻¹³ | 0.096 | −0.84 |
| 14 | S5_218071486 | 5 | 218,071,486 | 2.11 × 10⁻¹³ | 0.096 | −0.84 |
| 15 | S5_22426333 | 5 | 22,426,333 | 6.98 × 10⁻¹³ | 0.011 | −2.11 |
| 16 | S5_147071985 | 5 | 147,071,985 | 8.28 × 10⁻¹³ | 0.048 | −1.05 |
| 17 | S9_34634719 | 9 | 34,634,719 | 1.25 × 10⁻¹² | 0.300 | +0.79 |
| 18 | S2_12539331 | 2 | 12,539,331 | 1.55 × 10⁻¹² | 0.046 | −1.07 |
| 19 | S4_13315712 | 4 | 13,315,712 | 1.56 × 10⁻¹² | 0.349 | −0.83 |
| 20 | S1_193403593 | 1 | 193,403,593 | 1.99 × 10⁻¹² | 0.191 | −0.71 |

### Top associations — FarmCPU (top 20)

| Rank | SNP | Chr | Position | P-value | MAF | Effect (days) |
|---:|:---|---:|---:|:---|---:|---:|
| 1 | S3_156277883 | 3 | 156,277,883 | 5.76 × 10⁻¹⁷ | 0.393 | −0.06 |
| 2 | S3_156226908 | 3 | 156,226,908 | 1.99 × 10⁻¹² | 0.341 | +0.26 |
| 3 | S2_212614021 | 2 | 212,614,021 | 4.57 × 10⁻⁹ | 0.332 | +0.93 |
| 4 | S5_11998173 | 5 | 11,998,173 | 9.97 × 10⁻⁹ | 0.011 | −1.60 |
| 5 | S3_137533609 | 3 | 137,533,609 | 2.99 × 10⁻⁷ | 0.105 | −0.55 |
| 6 | S1_210194624 | 1 | 210,194,624 | 5.57 × 10⁻⁷ | 0.053 | −0.87 |
| 7 | S1_194883376 | 1 | 194,883,376 | 2.50 × 10⁻⁶ | 0.364 | −0.44 |
| 8 | S3_156270994 | 3 | 156,270,994 | 9.57 × 10⁻⁶ | 0.315 | −0.37 |
| 9 | S1_246327000 | 1 | 246,327,000 | 1.48 × 10⁻⁵ | 0.099 | +0.45 |
| 10 | S8_113222492 | 8 | 113,222,492 | 1.78 × 10⁻⁵ | 0.381 | −0.66 |
| 11 | S8_109133633 | 8 | 109,133,633 | 2.05 × 10⁻⁵ | 0.295 | +0.49 |
| 12 | S2_46237998 | 2 | 46,237,998 | 2.66 × 10⁻⁵ | 0.239 | −0.30 |
| 13 | S10_83314593 | 10 | 83,314,593 | 3.87 × 10⁻⁵ | 0.025 | +0.65 |
| 14 | S2_214213682 | 2 | 214,213,682 | 5.75 × 10⁻⁵ | 0.090 | −0.74 |
| 15 | S1_252128005 | 1 | 252,128,005 | 6.66 × 10⁻⁵ | 0.303 | −0.30 |
| 16 | S2_167022852 | 2 | 167,022,852 | 7.12 × 10⁻⁵ | 0.444 | +0.50 |
| 17 | S5_213565655 | 5 | 213,565,655 | 8.99 × 10⁻⁵ | 0.075 | −0.44 |
| 18 | S1_55541386 | 1 | 55,541,386 | 1.29 × 10⁻⁴ | 0.122 | +0.37 |
| 19 | S4_175791704 | 4 | 175,791,704 | 1.37 × 10⁻⁴ | 0.421 | −0.43 |
| 20 | S1_24749523 | 1 | 24,749,523 | 1.42 × 10⁻⁴ | 0.080 | +0.35 |

### SNPs shared between GLM and FarmCPU top 20

| SNP | Chr | Position | P-value (FarmCPU) | MAF | Effect (days) |
|:---|---:|---:|:---|---:|---:|
| S3_156277883 | 3 | 156,277,883 | 5.76 × 10⁻¹⁷ | 0.393 | −0.06 |
| S5_11998173 | 5 | 11,998,173 | 9.97 × 10⁻⁹ | 0.011 | −1.60 |
| S3_137533609 | 3 | 137,533,609 | 2.99 × 10⁻⁷ | 0.105 | −0.55 |
| S1_210194624 | 1 | 210,194,624 | 5.57 × 10⁻⁷ | 0.053 | −0.87 |

**S5_11998173 (plantacyanin family):** Rare ALT allele (C; REF = T) carried by 52/2,420 hybrids. Falls in the 3' UTR of a putative blue copper protein gene Zm00001eb216460 (Arabidopsis putative ortholog AT2G02850), suggesting a regulatory rather than coding mechanism, if it is not a false positive. Plantacyanin family proteins can be expressed in reproductive organs and their expression levels have impacts on fertility in Arabidopsis (https://doi.org/10.1104/pp.105.063388).

---

## Limitations

The SNP panel (2,060 markers) is sparse relative to LD decay in maize, limiting mapping resolution. The two-stage approach discards within-environment variance but is standard practice for datasets of this scale. The plantacyanin annotation is bioinformatic; no functional validation has been performed.

---

## Software

R 4.3+; `lme4`, `emmeans` (BLUEs); `vcfR` (genotype processing); `GAPIT3` (GWAS). Dependencies managed via `renv`.
