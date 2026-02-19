# Load packages

library(GAPIT)
library(tidyverse)
library(ggplot2)

# Load phenotype file
flowering_blues_for_gwas <- read_csv("~/g2f-maize-quantitative-genetics/data/processed/flowering_BLUEs_for_GWAS.csv") %>%
  as.data.frame()
dim(flowering_blues_for_gwas)
class(flowering_blues_for_gwas)

# Load genotype file 
genotype_for_gwas <- read_csv("~/g2f-maize-quantitative-genetics/data/processed/genotypes_for_GWAS.csv") %>%
  as.data.frame()

# Load genotype map file

genotype_map <- read_csv("~/g2f-maize-quantitative-genetics/data/processed/SNP_map_for_GWAS.csv") %>%
  as.data.frame()
dim(genotype_for_gwas)

# Set working directory to output folder
setwd("01-gwas-flowering-time/results/gwas-output/")

# do colnames match?
all(genotype_for_gwas$Taxa == flowering_blues_for_gwas$Taxa)


GAPIT_glm <- GAPIT(
  Y = flowering_blues_for_gwas,
  GD = genotype_for_gwas,
  GM = genotype_map,
  PCA.total = 4,
  model = "GLM"
)

sum(is.na(genotype_for_gwas))

# Load GWAS results
gwas_results <- read_csv("~/g2f-maize-quantitative-genetics/01-gwas-flowering-time/results/gwas-output/GAPIT.Association.GWAS_Results.GLM.Flowering_BLUE(NYC).csv")

# Check structure
head(gwas_results)
dim(gwas_results)

# Summary of P-values
summary(gwas_results$P.value)

# How many significant (P < 0.05)?
sum(gwas_results$P.value < 0.05, na.rm = TRUE)

# Top 10 most significant SNPs
gwas_results %>%
  arrange(P.value) %>%
  head(10)

GAPIT_farmCPU <- GAPIT(
  Y = flowering_blues_for_gwas,
  GD = genotype_for_gwas,
  GM = genotype_map,
  PCA.total = 4,
  model = "farmCPU",
)

# Load GWAS results
gwas_results_farmcpu <- read_csv("~/g2f-maize-quantitative-genetics/01-gwas-flowering-time/results/gwas-output/GAPIT.Association.GWAS_Results.FarmCPU.Flowering_BLUE(NYC).csv")

# Check structure
head(gwas_results_farmcpu)
dim(gwas_results_farmcpu)

# Summary of P-values
summary(gwas_results_farmcpu$P.value)

# How many significant (P < 0.05)?
sum(gwas_results_farmcpu$P.value < 0.05, na.rm = TRUE)

# Top 10 most significant SNPs
gwas_results_farmcpu %>%
  arrange(P.value) %>%
  head(10)

# Create summary of key QTL
key_qtl <- gwas_results_farmcpu %>%
  filter(P.value < 1e-5) %>%
  arrange(P.value) %>%
  select(SNP, Chr, Pos, P.value, MAF, Effect)

# Group by chromosome
key_qtl %>%
  group_by(Chr) %>%
  summarise(
    n_SNPs = n(),
    min_pval = min(P.value),
    mean_effect = mean(abs(Effect))
  ) %>%
  arrange(min_pval)

# Which SNPs appear in BOTH top 20?
top_glm <- gwas_results %>% arrange(P.value) %>% head(20)
top_farmcpu <- gwas_results_farmcpu %>% arrange(P.value) %>% head(20)

shared_snps <- intersect(top_glm$SNP, top_farmcpu$SNP)
length(shared_snps)

# Look at these shared SNPs
gwas_results_farmcpu %>%
  filter(SNP %in% shared_snps) %>%
  select(SNP, Chr, Pos, P.value, MAF, Effect) %>%
  arrange(P.value)

# Distribution of significant SNPs by chromosome
farmcpu_sig_by_chr <- table(gwas_results_farmcpu$Chr[gwas_results_farmcpu$P.value < 0.05])
glm_sig_by_chr <- table(gwas_results$Chr[gwas_results$P.value < 0.05])

# Compare
data.frame(
  Chromosome = 1:10,
  GLM_hits = as.numeric(glm_sig_by_chr),
  FarmCPU_hits = as.numeric(farmcpu_sig_by_chr)
)

# Get genotypes for your top SNP
top_snp <- "S5_11998173"

# Extract genotype for this SNP
geno_data <- data.frame(
  Taxa = genotype_for_gwas$Taxa,
  Genotype = genotype_for_gwas[[top_snp]]
)

# Merge with phenotypes
snp_pheno <- flowering_blues_for_gwas %>%
  left_join(geno_data, by = "Taxa")

# Look at distribution
table(snp_pheno$Genotype)

# Calculate mean flowering time per genotype
snp_pheno %>%
  group_by(Genotype) %>%
  summarise(
    n = n(),
    mean_flowering = mean(Flowering_BLUE, na.rm = TRUE),
    sd_flowering = sd(Flowering_BLUE, na.rm = TRUE)
  )


# Figure for S5_11998173)
s5_plot <- ggplot(snp_pheno, aes(x = factor(Genotype), y = Flowering_BLUE)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.3, fill = "darkblue", alpha = 0.7) +
  geom_jitter(alpha = 0.2, width = 0.1) +
  stat_summary(fun = mean, geom = "point", 
               shape = 23, size = 4, fill = "red") +
  labs(title = "Rare Allele with Major Effect: S5_11998173",
       subtitle = "Plantacyanin gene - 4.1 day flowering advancement (n=52 carriers)",
       x = "Genotype (0=reference, 1=rare allele)",
       y = "Flowering Time BLUE (days)",
       caption = "P = 9.97e-9 | Effect = -4.1 days | MAF = 1.07%") +
  annotate("text", x = 1.5, y = 66, 
           label = "Effect = -4.1 days", size = 6, color = "red", fontface = "bold") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("~/g2f-maize-quantitative-genetics/01-gwas-flowering-time/figures/plantacyanin_rare_allele_effect.png",
       s5_plot, width = 8, height = 6)
