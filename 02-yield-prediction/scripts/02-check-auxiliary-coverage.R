# Exercise 2: How well do the auxiliary files cover the trait environments?
# Purpose: every auxiliary file joins on Env. Before choosing features, find out
#          which environments (and how many yield plots) each file would cost us,
#          and whether gaps cluster in particular years.

library(tidyverse)

dir.create("02-yield-prediction/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("02-yield-prediction/results", recursive = TRUE, showWarnings = FALSE)

check <- function(label, observed, expected) {
  status <- if (isTRUE(all.equal(observed, expected, tolerance = 1e-3))) "OK" else "MISMATCH"
  cat(sprintf("[%s] %s: observed %s, expected %s\n", status, label, observed, expected))
}

#### Trait environments (the reference set) ####

trait_data <- read_csv("data/raw/1_Training_Trait_Data_2014_2023.csv",
                       col_select = c(Env, Year, Yield_Mg_ha),
                       show_col_types = FALSE)

trait_envs <- trait_data |>
  group_by(Env, Year) |>
  summarise(yield_plots = sum(!is.na(Yield_Mg_ha)), .groups = "drop")

cat(sprintf("Trait file: %d environments, %s plots with yield\n",
            nrow(trait_envs), format(sum(trait_envs$yield_plots), big.mark = ",")))

#### Auxiliary files ####

aux_files <- c(
  meta    = "data/raw/2_Training_Meta_Data_2014_2023.csv",
  soil    = "data/raw/3_Training_Soil_Data_2015_2023.csv",
  weather = "data/raw/4_Training_Weather_Data_2014_2023_seasons_only.csv",
  ec      = "data/raw/6_Training_EC_Data_2014_2023.csv"
)

# Only the Env column is needed here; col_select avoids parsing all 655 EC columns
read_envs <- function(path) {
  read_csv(path, col_select = Env, show_col_types = FALSE)
}

aux_envs <- map(aux_files, read_envs)

# Duplicated Env rows would multiply trait rows in a join, so flag them early.
# (A weather file with one row per day or per stage is expected to repeat Env.)
iwalk(aux_envs, \(df, file) {
  cat(sprintf("%-8s %6d rows, %4d distinct Env, max rows per Env = %d\n",
              file, nrow(df), n_distinct(df$Env), max(count(df, Env)$n)))
})

#### Coverage summary ####

# Environment coverage treats every env equally; plot coverage weights by how
# much yield data would actually be lost, which is what matters for training
coverage <- imap(aux_envs, \(df, file) {
  covered <- trait_envs$Env %in% df$Env
  tibble(
    file              = file,
    envs_in_file      = n_distinct(df$Env),
    envs_matched      = sum(covered),
    pct_trait_envs    = round(100 * mean(covered), 1),
    pct_yield_plots   = round(100 * sum(trait_envs$yield_plots[covered]) /
                                sum(trait_envs$yield_plots), 1),
    envs_not_in_trait = sum(!unique(df$Env) %in% trait_envs$Env)
  )
}) |>
  list_rbind()

print(coverage)
write_csv(coverage, "02-yield-prediction/results/aux_file_coverage.csv")

expected <- tribble(
  ~file,     ~envs_in_file, ~pct_trait_envs,
  "meta",    272,           100,
  "soil",    186,           68.4,
  "weather", 269,           98.9,
  "ec",      241,           88.6
)
pwalk(expected, \(file, envs_in_file, pct_trait_envs) {
  row <- filter(coverage, .data$file == .env$file)
  check(sprintf("%s envs in file", file), row$envs_in_file, envs_in_file)
  check(sprintf("%s %% of trait envs", file), row$pct_trait_envs, pct_trait_envs)
})

#### Coverage by year ####

# Gaps that cluster in one year (e.g. soil before 2015) need a different fix
# than gaps scattered across years
coverage_by_year <- imap(aux_envs, \(df, file) {
  trait_envs |>
    mutate(covered = Env %in% df$Env) |>
    group_by(Year) |>
    summarise(pct_envs = round(100 * mean(covered), 1), .groups = "drop") |>
    mutate(file = file)
}) |>
  list_rbind()

coverage_by_year |>
  pivot_wider(names_from = file, values_from = pct_envs) |>
  print()

coverage_heatmap <- coverage_by_year |>
  mutate(file = factor(file, levels = names(aux_files))) |>
  ggplot(aes(x = factor(Year), y = file, fill = pct_envs)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(pct_envs)), size = 3.5) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100)) +
  labs(title = "Auxiliary file coverage of trait environments by year",
       x = "Year", y = "File", fill = "% envs\ncovered") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

ggsave("02-yield-prediction/figures/aux_coverage_by_year.png",
       coverage_heatmap, width = 9, height = 4, dpi = 300)

#### Which environments are missing? ####

# Named lists make the imputation decisions in Exercise 5 concrete
missing_envs <- imap(aux_envs, \(df, file) {
  trait_envs |>
    filter(!Env %in% df$Env) |>
    mutate(file = file)
}) |>
  list_rbind() |>
  select(file, Year, Env, yield_plots) |>
  arrange(file, Year, Env)

missing_envs |> count(file, Year) |> print(n = Inf)
write_csv(missing_envs, "02-yield-prediction/results/aux_missing_envs.csv")

cat("Exercise 2 (coverage) done: results/aux_file_coverage.csv, results/aux_missing_envs.csv,",
    "figures/aux_coverage_by_year.png saved\n")

#### Exercise 2 remaining: meta feature audit — your turn ####
# Steps (details in chat / CLAUDE.md):
#   2a. Clean the meta column names
#   2b. Profile every column: type, % missing, n distinct, example values
#   2c. Classify keep / derive / drop — known before planting? varies across envs?
#       an ID in disguise? redundant?
#   2d. Derive candidates: field lat/long, irrigated flag, previous crop (grouped),
#       planting DOY (check Date_Planted is constant within Env in the trait file)
# Output: 02-yield-prediction/results/meta_feature_audit.csv
