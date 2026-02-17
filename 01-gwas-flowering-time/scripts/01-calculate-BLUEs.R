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

hybrid_blues <- emmeans(model_blues, specs = ~ Hybrid, pbkrtest.limit = 63189)

# Convert to data frame
blues_df <- as.data.frame(hybrid_blues) %>%
  select(Hybrid, emmean, SE)

# Examine
head(blues_df, 10)  # Look at first 10
nrow(blues_df)      # How many hybrids?

# Summary statistics
summary(blues_df$emmean)  # Range of BLUEs
summary(blues_df$SE)      # Range of standard errors

# Create a histogram showing the distribution of your BLUE values

hist(blues_df$emmean)

blues_hist <- ggplot(blues_df, aes(x = emmean)) +
  theme_bw() +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  geom_vline(xintercept = median(blues_df$emmean), 
                     color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 63, y = 210, label = "Median") +
  labs(title = "Histogram of BLUEs calculated means", x = "Flowering Time BLUE (Days After Planting)", y = "Frequency") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
ggsave("01-gwas-flowering-time/figures/BLUE_distribution.png", 
       blues_hist, width = 8, height = 5)
  


# Is the distribution approximately normal?
# Any signs of bimodality (two peaks)?
# Where do most hybrids cluster?

# Show that hybrids tested in more environments have lower uncertainty

# Group by hybrid, distinct environments and total observations
env_by_hybrid <- gwas_data %>%
  group_by(Hybrid) %>%
  summarize(n_environments = n_distinct(Env),
            n_observations = n())

# Join with blues DF by Hybrid
blues_df <- blues_df %>%
  left_join(env_by_hybrid, by = "Hybrid")

# Find which n_environments and after are SE<1

se_by_env <- blues_df %>%
  group_by(n_environments) %>%
  summarise(median_SE_per_n_env = median(SE)) %>%
  filter(median_SE_per_n_env<1.1) %>%
  arrange(desc(median_SE_per_n_env))

head(se_by_env)
# After 10 distinct environments, median SE drops below 1

# Make the scatter plot
scatter_plot <- ggplot(blues_df, aes(x = n_environments, y = SE)) +
  geom_point() +
  geom_smooth() +
  labs(title = "Scatter plot Distinct Environments by Standard Error",
       x = "Distinct Environments per Hybrid",
       y = "Flowering time prediction Standard Error") +
  theme_bw() +
theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  annotate("text",x = 31.5, y = 1.04, label = "Median SE<1 after 10 Distinct Environments")
ggsave("01-gwas-flowering-time/figures/SE_vs_environments.png", 
       scatter_plot, width = 8, height = 5)

# Histogram of n_environments
# Calculate median per Hybrid
n_env_median <- median(blues_df$n_environments)

n_env_hist <- ggplot(blues_df, aes(x = n_environments)) +
  theme_bw() +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  geom_vline(xintercept = median(n_env_median), 
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 27, y = 210, label = paste("Median = ", n_env_median)) +
  labs(title = "Histogram of distinct Environments per Hybrid", x = "Distinct Environments per Hybrid", y = "Frequency") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
ggsave("01-gwas-flowering-time/figures/n_env_distribution.png", 
       n_env_hist, width = 8, height = 5)

# How many hybrids were tested in fewer than 10 environments?

blues_df %>%
  filter(n_environments < 10) %>%
  summarise(hybrids_below_10_env = n_distinct(Hybrid))
 
# Select Hybrid and emmean and rename to Taxa and Flowering_BLUE for GWAS

flowering_blues_for_gwas <- blues_df %>%
  select(Hybrid, emmean) %>%
  rename(Taxa = Hybrid,
         Flowering_BLUE = emmean)

# How many hybrids for GWAS?
nrow(flowering_blues_for_gwas)

# Any NAs?
sum(is.na(flowering_blues_for_gwas))
