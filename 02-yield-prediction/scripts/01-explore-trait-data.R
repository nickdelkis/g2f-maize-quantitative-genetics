# Exercise 1: Orient in the trait data before any yield modelling
# Purpose: confirm size, time span and yield coverage of the main trait file,
#          and see how plots / environments / hybrids are spread across years.

library(tidyverse)

dir.create("02-yield-prediction/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("02-yield-prediction/results", recursive = TRUE, showWarnings = FALSE)

# Compare against numbers already recorded in CLAUDE.md, so a re-run flags drift
check <- function(label, observed, expected) {
  status <- if (isTRUE(all.equal(observed, expected, tolerance = 1e-3))) "OK" else "MISMATCH"
  cat(sprintf("[%s] %s: observed %s, expected %s\n", status, label, observed, expected))
}

#### Load ####

# guess_max raised because several columns (Range, Pass, Twt_kg_m3) are NA for
# the first thousands of rows, which makes readr guess the wrong type
trait_data <- read_csv("data/raw/1_Training_Trait_Data_2014_2023.csv",
                       guess_max = 50000, show_col_types = FALSE)

glimpse(trait_data)

#### Overall size ####

cat(sprintf("Trait file: %s rows x %d columns\n",
            format(nrow(trait_data), big.mark = ","), ncol(trait_data)))
check("Rows", nrow(trait_data), 173960)
check("Columns", ncol(trait_data), 26)

overview <- trait_data |>
  summarise(
    years        = n_distinct(Year),
    first_year   = min(Year),
    last_year    = max(Year),
    environments = n_distinct(Env),
    locations    = n_distinct(Field_Location),
    hybrids      = n_distinct(Hybrid)
  )
print(overview)
check("Environments", overview$environments, 272)

#### Yield coverage ####

yield_pct <- mean(!is.na(trait_data$Yield_Mg_ha)) * 100
cat(sprintf("Yield_Mg_ha non-missing: %.1f%%\n", yield_pct))
check("Yield coverage (%)", round(yield_pct, 1), 94.8)

summary(trait_data$Yield_Mg_ha)

# Environments with no yield at all can't contribute to training and would
# silently vanish in later joins, so name them now
envs_without_yield <- trait_data |>
  group_by(Env) |>
  summarise(n_yield = sum(!is.na(Yield_Mg_ha))) |>
  filter(n_yield == 0)
cat(sprintf("Environments with zero yield records: %d\n", nrow(envs_without_yield)))
print(envs_without_yield$Env)

#### By year ####

year_overview <- trait_data |>
  group_by(Year) |>
  summarise(
    plots          = n(),
    environments   = n_distinct(Env),
    hybrids        = n_distinct(Hybrid),
    yield_plots    = sum(!is.na(Yield_Mg_ha)),
    yield_pct      = round(100 * yield_plots / plots, 1),
    mean_yield     = round(mean(Yield_Mg_ha, na.rm = TRUE), 2)
  )
print(year_overview)

write_csv(year_overview, "02-yield-prediction/results/trait_year_overview.csv")

#### Figure: plots per year, split by whether yield was recorded ####

yield_by_year_plot <- trait_data |>
  mutate(yield_recorded = if_else(is.na(Yield_Mg_ha), "Missing", "Recorded")) |>
  count(Year, yield_recorded) |>
  ggplot(aes(x = factor(Year), y = n, fill = yield_recorded)) +
  geom_col() +
  scale_fill_manual(values = c(Recorded = "steelblue", Missing = "grey70")) +
  labs(title = "Plots per year and yield availability",
       x = "Year", y = "Number of plots", fill = "Yield_Mg_ha") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

ggsave("02-yield-prediction/figures/plots_per_year_yield_coverage.png",
       yield_by_year_plot, width = 8, height = 5, dpi = 300)

cat("Exercise 1 done: results/trait_year_overview.csv and figures/plots_per_year_yield_coverage.png saved\n")
