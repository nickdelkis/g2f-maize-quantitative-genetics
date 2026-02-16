# Purpose: Calculate environment-adjusted phenotypes for GWAS using two-stage approach

# Load packages
library(tidyverse)
library(lme4)
library(emmeans)  # For extracting marginal means
library(ggplot2)

# Load GWAS data

gwas_data <- read_csv("data/processed/gwas_phenotypes_2018-2021.csv")

glimpse(gwas_data)

# How many observations per hybrid?
gwas_data %>%
  count(Hybrid) %>%
  summary()

# How many observations per environment?
gwas_data %>%
  count(Env) %>%
  summary()

# Is data balanced?
# No

model_blues <- lmer(data = gwas_data,
                    REML = TRUE,
                    formula = Pollen_DAP_days ~ Hybrid + (1|Env))

summary(model_blues)


### Step 1: `emmeans(model_blues, specs = ~ Hybrid)`

# **"specs = ~ Hybrid"** means: *"Give me the estimated marginal mean for each level of Hybrid"*
  
#  **What it does internally:**
#  1. Takes your fitted model
#  2. For each hybrid, calculates: **"What is the predicted flowering time for this hybrid, averaged across all environments?"**
#  3. Accounts for the random environment effects
#  4. Returns the **absolute values** (not contrasts)

# **Think of it as:**
#  ```
# For each Hybrid:
#   BLUE = Intercept + Hybrid_effect + average(Environment_effects)

# Extract hybrid BLUEs using emmeans
library(emmeans)

hybrid_blues <- emmeans(model_blues, specs = ~ Hybrid, pbkrtest.limit = 63189)

# Convert to data frame
blues_df <- as.data.frame(hybrid_blues)

# Examine
head(blues_df, 10)  # Look at first 10
nrow(blues_df)      # How many hybrids?

# Summary statistics
summary(blues_df$emmean)  # Range of BLUEs
summary(blues_df$SE)      # Range of standard errors

# Create a histogram showing the distribution of your BLUE values

hist(blues_df$emmean, bin)