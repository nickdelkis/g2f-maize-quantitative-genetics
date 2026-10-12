# Shared helpers for the yield-prediction exercise scripts.
# Source from the project root: source("02-yield-prediction/scripts/utils.R")

dir.create("02-yield-prediction/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("02-yield-prediction/results", recursive = TRUE, showWarnings = FALSE)

# Compare against numbers already recorded in CLAUDE.md, so a re-run flags drift.
# Counts must match exactly; percentages are recorded to 1 decimal, so allow
# half a unit of rounding (absolute, not relative — a relative tolerance would
# let a large row count drift by hundreds unnoticed)
check <- function(label, observed, expected, abs_tol = 0) {
  status <- if (abs(observed - expected) <= abs_tol) "OK" else "MISMATCH"
  cat(sprintf("[%s] %s: observed %s, expected %s\n", status, label, observed, expected))
}
