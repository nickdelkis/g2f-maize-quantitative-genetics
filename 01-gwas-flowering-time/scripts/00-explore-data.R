# Load packages
library(tidyverse)

# Load data
trait_data <- read_csv("data/raw/1_Training_Trait_Data_2014_2023.csv")

# Get structure of dataset
glimpse(trait_data)

# summary of dataset
summary(trait_data)

# store column names
colnames <- colnames(trait_data)

# years in dataset
years <- unique(trait_data$Year)

#environments in dataset
envs <- unique(trait_data$Env) %>%
  length()

# unique hybrids
hybrids <- unique(trait_data$Hybrid) %>% 
  length()

# count by year
count(trait_data, Year)

#### Flowering Time Quality ####

# How many records have Pollen_DAP_days?
pollen_dap_days <- trait_data %>%
  filter(!is.na(Pollen_DAP_days)) %>%
  count() %>%
  as.numeric()

# How many have `Silk_DAP_days`?
silk_dap_days <- trait_data %>%
  filter(!is.na(Silk_DAP_days)) %>%
  count() %>%
  as.numeric()

# Distribution of values?
#Pollen DAP days
summary(trait_data$Pollen_DAP_days)
hist(trait_data$Pollen_DAP_days)

# Silk DAP days
summary(trait_data$Silk_DAP_days)
hist(trait_data$Silk_DAP_days)

# See outliers visually
boxplot(trait_data$Pollen_DAP_days, main = "Pollen DAP Distribution")
boxplot(trait_data$Silk_DAP_days, main = "Silk DAP Distribution")

# Calculate outlier thresholds using IQR method
pollen_q1 <- quantile(trait_data$Pollen_DAP_days, 0.25, na.rm = TRUE)
pollen_q3 <- quantile(trait_data$Pollen_DAP_days, 0.75, na.rm = TRUE)
pollen_iqr <- pollen_q3 - pollen_q1

# Outliers are beyond 1.5 * IQR from quartiles
pollen_lower <- pollen_q1 - 1.5 * pollen_iqr
pollen_upper <- pollen_q3 + 1.5 * pollen_iqr

# Count outliers
trait_data %>%
  filter(Pollen_DAP_days < pollen_lower | Pollen_DAP_days > pollen_upper) %>%
  count()

#### Completeness Summary by year ####

# Create completeness summary by year
year_summary <- trait_data %>%
  group_by(Year) %>%
  summarise(
    total_plots = n(),
    has_pollen = sum(!is.na(Pollen_DAP_days)),
    has_silk = sum(!is.na(Silk_DAP_days)),
    pollen_pct = round(100 * has_pollen / total_plots, 1),
    silk_pct = round(100 * has_silk / total_plots, 1),
    unique_hybrids = n_distinct(Hybrid)
  ) %>%
  arrange(desc(total_plots))

print(year_summary)

# See which environments have extreme values
outliers <- trait_data %>%
  filter(Pollen_DAP_days < pollen_lower | Pollen_DAP_days > pollen_upper)

# Which environments/years are these from?
outliers %>%
  count(Year, Env) %>%
  arrange(desc(n))

# Are they clustered in specific environments?
outliers %>%
  count(Env) %>%
  arrange(desc(n))

# Are these environments outliers in OTHER traits too?
trait_data %>%
  filter(Env %in% c("GEH1_2020", "NEH2_2015", "KSH1_2016")) %>%
  summary()

#### Check for genotype availability ####
# Load genotyped hybrid names from VCF
vcf_lines <- read_lines("data/raw/5_Genotype_Data_All_2014_2025_Hybrids.vcf", 
                        n_max = 500)

header_line <- vcf_lines[str_detect(vcf_lines, "^#CHROM")]

genotyped_hybrids <- header_line %>%
  str_split("\t") %>%
  pluck(1) %>%
  tail(-9)

length(genotyped_hybrids)

# Check overlap with trait data
trait_data_with_geno <- trait_data %>%
  filter(Hybrid %in% genotyped_hybrids)

# How much data do we keep?
nrow(trait_data_with_geno) / nrow(trait_data) * 100

# By year
trait_data_with_geno %>%
  filter(!is.na(Pollen_DAP_days)) %>%  # Only those with flowering data
  group_by(Year) %>%
  summarise(
    plots = n(),
    hybrids = n_distinct(Hybrid),
    environments = n_distinct(Env)
  ) %>%
  arrange(desc(plots))

# Filter to your selected years with complete data
gwas_data <- trait_data_with_geno %>%
  filter(Year %in% c(2018, 2019, 2020, 2021)) %>%
  filter(!is.na(Pollen_DAP_days))  # Must have flowering time

# Check what you have
nrow(gwas_data)
n_distinct(gwas_data$Hybrid)
n_distinct(gwas_data$Env)

# Summary by year
gwas_data %>%
  group_by(Year) %>%
  summarise(
    plots = n(),
    hybrids = n_distinct(Hybrid),
    environments = n_distinct(Env),
    mean_pollen = mean(Pollen_DAP_days, na.rm = TRUE),
    sd_pollen = sd(Pollen_DAP_days, na.rm = TRUE)
  )

# Visualize why you chose these years
year_summary <- trait_data_with_geno %>%
  filter(!is.na(Pollen_DAP_days)) %>%
  group_by(Year) %>%
  summarise(
    plots = n(),
    hybrids = n_distinct(Hybrid)
  )

# Plot hybrids per year
ggplot(year_summary, aes(x = Year, y = hybrids)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Unique Hybrids per Year (with flowering data & genotypes)",
       subtitle = "2018-2021 have >1000 hybrids; 2016-2017,2022-2023 too few",
       y = "Number of Unique Hybrids",
       x = "Year") +
  scale_x_continuous("Year", labels = as.character(year_summary$Year), breaks = year_summary$Year)
  theme_minimal()

# Save the plot
ggsave("01-gwas-flowering-time/figures/year_selection.png", width = 8, height = 5)

# Save for next analysis step
write_csv(gwas_data, "data/processed/gwas_phenotypes_2018-2021.csv")

# === DECISION: Use 2018-2021 for GWAS ===
# 
# Rationale:
# - 63,189 plots with flowering time data
# - 1,039-1,179 unique hybrids per year (>1000 threshold)
# - 96 total environments across 4 years
# - 97.8% of data has matching genotypes
# - 2022-2023 excluded: only ~550 hybrids (insufficient genetic diversity)
# 
# Next steps:
# - Calculate environment means (BLUEs/BLUPs)
# - Match with genotype data
# - Run GWAS analysis