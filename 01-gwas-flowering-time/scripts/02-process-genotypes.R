# Load vcfR and GAPIT
library(vcfR)
library(GAPIT)
library(tidyverse)
library(ggplot2)

# Read vcf file
vcf <- read.vcfR("data/raw/5_Genotype_Data_All_2014_2025_Hybrids.vcf")

head(vcf@fix)
# What chromosomes are represented
chromosomes <- vcf@fix %>%
  as.data.frame() %>%
  select(CHROM) %>%
  unique()
chromosomes

# What does the genotype matrix look like?
genotype_matrix <- vcf@gt %>%
  as.data.frame()

head(genotype_matrix,5)

head(vcf@gt,c(5,5))

#### Extract genotypes and convert to numeric dosage matrix ####

gt_matrix <- extract.gt(vcf, element = "GT")

# Check what you got
dim(gt_matrix)  # Should be 2425 rows (SNPs) × 5899 cols (hybrids)
head(gt_matrix[, 1:5])  # Look at first few

# Convert to numeric dosage (0, 1, 2)
# Hint: vcfR has a parameter for this!
gt_numeric <- extract.gt(vcf, 
                         element = "GT",
                         as.numeric = TRUE)  # This does the conversion!

# Check the conversion
head(gt_numeric[, 1:5])

#Percentage of NAs in dataset
missing_data <- sum(is.na(gt_numeric)) / length(gt_numeric) * 100
paste(missing_data, "% NAs in genotype matrix.")

# Distribution of genotypes 

table(gt_numeric)

gt_numeric[gt_numeric == 3] <- NA


#### Quality filters ####
# Remove SNPs with >10% missing data

# Calculate % missing per SNP (row-wise) and identify SNPs above 10% threshold
missing_per_snp <-rowSums(is.na(gt_numeric)) / ncol(gt_numeric) * 100 %>%
  as.vector()
gt_numeric <- gt_numeric %>%
  as.data.frame() %>%
  mutate(missing_percent = missing_per_snp,
         missing_above10 = missing_percent > 10)

# How many are above 10%?
table(gt_numeric$missing_above10)

# Remove them from the matrix

gt_numeric_above10 <- gt_numeric %>%
  filter(missing_above10 == F)

# cleaner approach to keep as matrix

# Start over, keep as matrix
gt_numeric <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
gt_numeric[gt_numeric == 3] <- NA

# Calculate missingness per SNP
missing_per_snp <- rowSums(is.na(gt_numeric)) / ncol(gt_numeric)

# Logical vector: which SNPs to KEEP
keep_snps <- missing_per_snp <= 0.10

# Filter the matrix
gt_numeric_filtered <- gt_numeric[keep_snps, ]

# Check
dim(gt_numeric_filtered)  # Should be 2224 × 5899
sum(keep_snps)  # Should be 2224

### Calculate minor allele frequency (MAF) and remove SNPs below 0.05 MAF ###

# allele frequency
allele_freq <- rowMeans(gt_numeric_filtered, na.rm = TRUE) / 2

# calculate MAF (pmin = element-wise minimum)
MAF <- pmin(allele_freq, 1 - allele_freq)

MAF_keep <- MAF >= 0.05

gt_numeric_filtered_maf <- gt_numeric_filtered[MAF_keep, ]
# Check
dim(gt_numeric_filtered_maf) 
sum(MAF_keep) 

sum(keep_snps) - sum(MAF_keep)

#### Match Genotypes to Phenotypes ####
# Goal: Subset genotype matrix to only the 2,420 hybrids that have phenotype data

# Load phenotype file
flowering_blues_for_gwas <- read_csv("data/processed/flowering_BLUEs_for_GWAS.csv")

# extract hybrids from pheno and genotype
pheno_hybrids <- flowering_blues_for_gwas$Taxa
geno_hybrids <- colnames(gt_numeric_filtered_maf)

# which hybrids from phenotype file are in genotype file?
hybrid_overlap <- geno_hybrids %in% pheno_hybrids
table(hybrid_overlap)

# Subset and keep only overlapping hybrids

gt_numeric_filtered_maf_subset <- gt_numeric_filtered_maf[, hybrid_overlap]
sum(is.na(gt_numeric_filtered_maf_subset))
# check dimensions
dim(gt_numeric_filtered_maf_subset)

# checks
length(colnames(gt_numeric_filtered_maf_subset))

# do colnames match?
all(colnames(gt_numeric_filtered_maf_subset) == flowering_blues_for_gwas$Taxa)

# Middle/mean impute missing SNP data

for(i in 1:nrow(gt_numeric_filtered_maf_subset)) {
  row_mean <- mean(gt_numeric_filtered_maf_subset[i, ], na.rm = TRUE)
  # Replace NAs with row mean
  gt_numeric_filtered_maf_subset[i, is.na(gt_numeric_filtered_maf_subset[i, ])] <- round(row_mean)
}

sum(is.na(gt_numeric_filtered_maf_subset))

# transpose and prepare for GAPIT - rownames are Taxa, first column)
gt_numeric_filtered_maf_subset_transpose <- t(gt_numeric_filtered_maf_subset) %>%
  as.data.frame() %>%
  mutate(Taxa = rownames(gt_numeric_filtered_maf_subset_transpose)) %>%
  select(Taxa, everything())

# Check how it looks
dim(gt_numeric_filtered_maf_subset_transpose)
head(gt_numeric_filtered_maf_subset_transpose, c(5,5))
sum(is.na(gt_numeric_filtered_maf_subset_transpose))
# Save file
write_csv(gt_numeric_filtered_maf_subset_transpose, file = "~/g2f-maize-quantitative-genetics/data/processed/genotypes_for_GWAS.csv")

# Create genotype map and filter for snps and maf
genotype_map <- vcf@fix %>%
  as.data.frame() %>%
  select(ID, CHROM, POS) %>%
  rename(SNP = ID,
         Chromosome = CHROM,
         Position = POS)

gneotype_map_keep_snps <- genotype_map[keep_snps, ]

genotype_map_keep_snps_maf <- gneotype_map_keep_snps[MAF_keep, ]

sum(is.na(genotype_map_keep_snps_maf))

# Save file 
write_csv(genotype_map_keep_snps_maf, file = "~/g2f-maize-quantitative-genetics/data/processed/SNP_map_for_GWAS.csv")
