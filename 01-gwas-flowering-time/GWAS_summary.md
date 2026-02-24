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

**Top associations:**
GLM

   SNP            Chr       Pos  P.value     MAF  nobs Effect `H&B.P.Value`
   <chr>        <dbl>     <dbl>    <dbl>   <dbl> <dbl>  <dbl>         <dbl>
 1 S5_11998173      5  11998173 2.22e-38 0.0107   2420 -3.81       4.58e-35
 2 S1_210194624     1 210194624 1.54e-29 0.0525   2420 -1.80       1.59e-26
 3 S3_137533609     3 137533609 7.21e-22 0.105    2420 -1.11       4.95e-19
 4 S5_206652970     5 206652970 4.47e-21 0.0114   2420 -2.67       2.30e-18
 5 S3_177563221     3 177563221 1.53e-20 0.00847  2420 -3.04       6.30e-18
 6 S5_15897858      5  15897858 3.01e-18 0.150    2420 -0.873      1.03e-15
 7 S4_178069434     4 178069434 1.46e-17 0.226    2420  2.35       4.28e-15
 8 S3_152992674     3 152992674 5.57e-15 0.323    2420 -0.786      1.44e-12
 9 S5_213160675     5 213160675 1.50e-14 0.125    2420 -0.859      3.44e-12
10 S3_156277883     3 156277883 3.50e-14 0.393    2420 -0.831      6.55e-12
11 S3_156277910     3 156277910 3.50e-14 0.393    2420 -0.831      6.55e-12
12 S2_23345860      2  23345860 4.45e-14 0.301    2420 -0.730      7.64e-12
13 S5_218069890     5 218069890 2.11e-13 0.0961   2420 -0.839      3.10e-11
14 S5_218071486     5 218071486 2.11e-13 0.0961   2420 -0.839      3.10e-11
15 S5_22426333      5  22426333 6.98e-13 0.0107   2420 -2.11       9.58e-11
16 S5_147071985     5 147071985 8.28e-13 0.0483   2420 -1.05       1.07e-10
17 S9_34634719      9  34634719 1.25e-12 0.300    2420  0.791      1.52e-10
18 S2_12539331      2  12539331 1.55e-12 0.0459   2420 -1.07       1.69e-10
19 S4_13315712      4  13315712 1.56e-12 0.349    2420 -0.825      1.69e-10
20 S1_193403593     1 193403593 1.99e-12 0.191    2420 -0.706      2.05e-10

FarmCPU
   SNP            Chr       Pos  P.value    MAF  nobs  Effect `H&B.P.Value`
   <chr>        <dbl>     <dbl>    <dbl>  <dbl> <dbl>   <dbl>         <dbl>
 1 S3_156277883     3 156277883 5.76e-17 0.393   2420 -0.0556      1.19e-13
 2 S3_156226908     3 156226908 1.99e-12 0.341   2420  0.263       2.05e- 9
 3 S2_212614021     2 212614021 4.57e- 9 0.332   2420  0.925       3.14e- 6
 4 S5_11998173      5  11998173 9.97e- 9 0.0107  2420 -1.60        5.13e- 6
 5 S3_137533609     3 137533609 2.99e- 7 0.105   2420 -0.547       1.23e- 4
 6 S1_210194624     1 210194624 5.57e- 7 0.0525  2420 -0.868       1.91e- 4
 7 S1_194883376     1 194883376 2.50e- 6 0.364   2420 -0.444       7.37e- 4
 8 S3_156270994     3 156270994 9.57e- 6 0.315   2420 -0.374       2.46e- 3
 9 S1_246327000     1 246327000 1.48e- 5 0.0994  2420  0.445       3.38e- 3
10 S8_113222492     8 113222492 1.78e- 5 0.381   2420 -0.656       3.66e- 3
11 S8_109133633     8 109133633 2.05e- 5 0.295   2420  0.494       3.83e- 3
12 S2_46237998      2  46237998 2.66e- 5 0.239   2420 -0.301       4.57e- 3
13 S10_83314593    10  83314593 3.87e- 5 0.0252  2420  0.649       6.14e- 3
14 S2_214213682     2 214213682 5.75e- 5 0.0895  2420 -0.737       8.46e- 3
15 S1_252128005     1 252128005 6.66e- 5 0.303   2420 -0.301       9.15e- 3
16 S2_167022852     2 167022852 7.12e- 5 0.444   2420  0.502       9.17e- 3
17 S5_213565655     5 213565655 8.99e- 5 0.0754  2420 -0.442       1.09e- 2
18 S1_55541386      1  55541386 1.29e- 4 0.122   2420  0.372       1.46e- 2
19 S4_175791704     4 175791704 1.37e- 4 0.421   2420 -0.425       1.46e- 2
20 S1_24749523      1  24749523 1.42e- 4 0.0795  2420  0.351       1.46e- 2

Shared SNPs
  SNP            Chr       Pos  P.value    MAF  Effect
  <chr>        <dbl>     <dbl>    <dbl>  <dbl>   <dbl>
1 S3_156277883     3 156277883 5.76e-17 0.393  -0.0556
2 S5_11998173      5  11998173 9.97e- 9 0.0107 -1.60  
3 S3_137533609     3 137533609 2.99e- 7 0.105  -0.547 
4 S1_210194624     1 210194624 5.57e- 7 0.0525 -0.868

**S5_11998173 (plantacyanin family):** Rare ALT allele (C; REF = T) carried by 52/2,420 hybrids. Falls in the 3' UTR of a putative blue copper protein gene Zm00001eb216460 (Arabidopsis putative ortholog AT2G02850), suggesting a regulatory rather than coding mechanism, if it is not a false positive. Plantacyanin family proteins can be expressed in reproductive organs and their expression levels have impacts on fertility in Arabidopsis (https://doi.org/10.1104/pp.105.063388).

---

## Limitations

The SNP panel (2,060 markers) is sparse relative to LD decay in maize, limiting mapping resolution. The two-stage approach discards within-environment variance but is standard practice for datasets of this scale. The plantacyanin annotation is bioinformatic; no functional validation has been performed.

---

## Software

R 4.3+; `lme4`, `emmeans` (BLUEs); `vcfR` (genotype processing); `GAPIT3` (GWAS). Dependencies managed via `renv`.
