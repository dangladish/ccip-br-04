########################################################################
## 03_uncertainty_report.R
##
## Compares uncertainty in coral cover and COTS counts across monitoring
## methods (manta tow, RS Surface, RS Deep, RHIS): a GBR-wide pooled
## summary plus per-reef detail. Reef-level inputs only; fits no models.
## Inputs:  data/derived/*.RDS (from 00_read_and_build_data.R) and the raw
##          RS Deep point CSV (needed for sub-reef grid-cell spread).
## Outputs: PNGs in output/figures/.
##
## Two uncertainty bands are drawn on every plot:
##   - SD across sub-reef sampling units (thin, pale): heterogeneity of the
##     units; does not shrink as units are added.
##   - SE of the mean = SD / sqrt(n) (thick, nested in the SD band): shrinks
##     as units are pooled. Computed here (add_se()), so it differs from the
##     analytic `_se` columns in reef_level.RDS.
##
## Sub-reef sampling unit: tows (manta, RS Surface, RS Deep COTS), ~100 m
## grid cells (RS Deep coral cover), surveys (RHIS).
## COTS is plotted as raw counts per tow/survey (sparse and right-skewed, so
## SD bands can fall below zero) and on a log1p scale.
## The pooled summary pools all sub-reef units directly rather than
## averaging reef means. Per-reef/pooled mean and SD logic is in
## R/uncertainty_shared.R (shared with the Shiny data build).
## Per-reef plots are facetted by method and sized for one A4 page.
########################################################################

library(tidyverse)
library(arrow)
library(patchwork)

source("R/uncertainty_shared.R")

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Load inputs ----

tow_level  <- readRDS("data/derived/tow_level.RDS")
reef_level <- readRDS("data/derived/reef_level.RDS")
reef_xwalk <- readRDS("data/derived/reef_xwalk.RDS")
rhis       <- readRDS("data/derived/rhis.RDS")

rs_deep_raw <- read_csv(
  "data/raw/reefscan_deep/reefcloud-point-classification-GRB-ReefScan-Deep-2026-08-13-filtered.csv",
  show_col_types = FALSE
)

# Reload RS Deep with the voyage exclusion applied, matching the canonical
# tables (the read above is unfiltered).
rs_deep_raw <- load_rs_deep_raw_3voyage_excluded()

deep_cell_summ_all <- compute_deep_cell_summary(rs_deep_raw, reef_xwalk)

# RHIS restricted to the project's reefs (required by the pooled functions)
rhis_scope <- rhis |> filter(reef_label %in% reef_xwalk$reef_label)

# ---- Helpers: SE column and dual SD/SE band layer ----

# SE of the mean from the SD and n returned by the shared compute_*() functions
add_se <- function(df) dplyr::mutate(df, se = sd / sqrt(n))

# Thin pale band = +/-1 SD, thick band = +/-1 SE. dodge_width = NULL turns
# dodging off (one row per method).
dual_band <- function(dodge_width = 0.7) {
  pos <- if (is.null(dodge_width) || dodge_width == 0) "identity" else position_dodge(width = dodge_width)
  list(
    geom_linerange(aes(xmin = mean - sd, xmax = mean + sd), position = pos, linewidth = 0.6, alpha = 0.35),
    geom_linerange(aes(xmin = mean - se, xmax = mean + se), position = pos, linewidth = 2.2),
    geom_point(position = pos, size = 0.9, color = "grey15")
  )
}

uncertainty_caption <- "Thin pale band = \u00b11 SD across sub-reef units (spread); thick band = \u00b11 SE of the mean (precision)"

count_breaks <- c(0, 1, 2, 5, 10, 20, 50)

# ---- GBR-wide pooled summary: coral + COTS ----

pooled_coral    <- compute_coral_uncertainty_pooled(tow_level, deep_cell_summ_all, rhis_scope) |> add_se()
pooled_cots     <- compute_cots_uncertainty_pooled(tow_level, rhis_scope) |> add_se()
pooled_cots_log <- compute_cots_uncertainty_pooled(tow_level, rhis_scope, log1p = TRUE) |> add_se()

cat("GBR-wide pooled coral cover (n = sub-reef units per method):\n")
pooled_coral |> select(method, mean, sd, se, n) |> print()
cat("\nGBR-wide pooled COTS counts (n = sub-reef units per method):\n")
pooled_cots |> select(method, mean, sd, se, n) |> print()

p_coral_gbr <- ggplot(pooled_coral, aes(y = method, x = mean, color = method)) +
  dual_band(dodge_width = NULL) +
  guides(color = "none") +
  labs(
    y = NULL,
    x = "Hard coral cover (%)",
    title = "GBR-wide coral cover: method comparison, all reefs pooled",
    subtitle = "Point = mean across every sub-reef unit, all reefs combined",
    caption = uncertainty_caption
  )

ggsave("output/figures/gbr_pooled_coral_cover_comparison.png", p_coral_gbr, width = 8, height = 3.5, dpi = 150)
cat("\nSaved output/figures/gbr_pooled_coral_cover_comparison.png\n")

p_cots_gbr <- ggplot(pooled_cots, aes(y = method, x = mean, color = method)) +
  dual_band(dodge_width = NULL) +
  guides(color = "none") +
  labs(
    y = NULL,
    x = "COTS observed (count per tow or survey)",
    title = "GBR-wide COTS counts: method comparison, all reefs pooled",
    subtitle = "Point = mean across every sub-reef unit, all reefs combined (SD band may go below 0 -- see caption)",
    caption = uncertainty_caption
  )

ggsave("output/figures/gbr_pooled_cots_detection_comparison.png", p_cots_gbr, width = 8, height = 3.5, dpi = 150)
cat("Saved output/figures/gbr_pooled_cots_detection_comparison.png\n")

p_cots_log_gbr <- ggplot(pooled_cots_log, aes(y = method, x = mean, color = method)) +
  dual_band(dodge_width = NULL) +
  guides(color = "none") +
  scale_x_continuous(breaks = log1p(count_breaks), labels = count_breaks) +
  labs(
    y = NULL,
    x = "COTS observed (count per tow or survey, log1p scale)",
    title = "GBR-wide COTS counts (log1p scale): method comparison, all reefs pooled",
    subtitle = "Mean \u00b11 SD/SE computed on log1p(count) per sub-reef unit, then plotted -- axis labelled in original count units",
    caption = uncertainty_caption
  )

ggsave("output/figures/gbr_pooled_cots_count_log1p.png", p_cots_log_gbr, width = 8, height = 3.5, dpi = 150)
cat("Saved output/figures/gbr_pooled_cots_count_log1p.png\n")

# ---- Combined pooled figure (coral over COTS log1p), minimal labels ----
p_coral_gbr_simple <- p_coral_gbr +
  labs(title = "GBR-wide Coral cover", subtitle = NULL, caption = NULL)

p_cots_log_gbr_simple <- p_cots_log_gbr +
  labs(
    title = "GBR-wide COTS counts (log1p scale)",
    subtitle = NULL,
    caption = NULL,
    x = "COTS observed (count per tow or survey, log1p scale)"
  )

p_gbr_pooled_combined <- p_coral_gbr_simple / p_cots_log_gbr_simple

ggsave(
  "output/figures/gbr_pooled_combined.png",
  p_gbr_pooled_combined, width = 8, height = 7, dpi = 150
)
cat("Saved output/figures/gbr_pooled_combined.png\n")

# ---- Per-reef coral cover: mean +/- SD/SE ----

combined_coral <- compute_coral_uncertainty_reef(tow_level, deep_cell_summ_all, rhis, reef_xwalk) |> add_se()

cat("\nCoral cover uncertainty summary -- reefs per method:\n")
combined_coral |> count(method) |> print()

# Order reefs by ID for a stable layout
reef_order <- sort(unique(combined_coral$reef_label), decreasing = TRUE)
combined_coral <- combined_coral |> mutate(reef_label = factor(reef_label, levels = reef_order))

p_coral <- ggplot(combined_coral, aes(y = reef_label, x = mean, color = method)) +
  dual_band(dodge_width = NULL) +
  facet_grid(method ~ ., scales = "free_y", space = "free_y") +
  guides(color = "none") +
  labs(
    y = "Reef (ordered by reef ID)",
    x = "Hard coral cover (%)",
    title = "Reef-level Coral cover"
  ) +
  theme(
    axis.text.y = element_text(size = 5),
    strip.text = element_text(size = 7),
    panel.spacing = unit(3, "pt")
  )

ggsave(
  "output/figures/reef_coral_cover_method_comparison_all_reefs_with_rhis.png",
  p_coral, width = 8.27, height = 11.69, dpi = 300, limitsize = FALSE
)

cat("Saved output/figures/reef_coral_cover_method_comparison_all_reefs_with_rhis.png\n")

# ---- Per-reef COTS counts ----

# tow_level COTS counts are already zero-filled where imagery/tow was reviewed.
# RS Deep covers few tows per reef; shown as-is.
combined_cots <- compute_cots_uncertainty_reef(tow_level, rhis, reef_xwalk) |> add_se()

cat("\nCOTS count uncertainty summary -- reefs per method:\n")
combined_cots |> count(method) |> print()
cat("(rs_deep n per reef is expected to be small and uneven)\n")

reef_order_cots <- sort(unique(combined_cots$reef_label), decreasing = TRUE)
combined_cots <- combined_cots |> mutate(reef_label = factor(reef_label, levels = reef_order_cots))

p_cots <- ggplot(combined_cots, aes(y = reef_label, x = mean, color = method)) +
  dual_band(dodge_width = NULL) +
  facet_grid(method ~ ., scales = "free_y", space = "free_y") +
  guides(color = "none") +
  labs(
    y = "Reef (ordered by reef ID)",
    x = "COTS observed (count per tow or survey)",
    title = "Reef-level COTS counts"
  ) +
  theme(
    axis.text.y = element_text(size = 5),
    strip.text = element_text(size = 7),
    panel.spacing = unit(3, "pt")
  )

ggsave(
  "output/figures/reef_cots_detection_method_comparison_all_reefs.png",
  p_cots, width = 8.27, height = 11.69, dpi = 300, limitsize = FALSE
)

cat("Saved output/figures/reef_cots_detection_method_comparison_all_reefs.png\n")

# ---- Per-reef COTS counts, log1p scale ----

# Mean/SD/SE are computed on log1p(count) per sub-reef unit, then aggregated
# (log1p of already-computed raw means/SDs would be undefined below -1).
# Axis is labelled in original count units.
combined_cots_log <- compute_cots_uncertainty_reef(tow_level, rhis, reef_xwalk, log1p = TRUE) |> add_se()

reef_order_cots_log <- sort(unique(combined_cots_log$reef_label), decreasing = TRUE)
combined_cots_log <- combined_cots_log |> mutate(reef_label = factor(reef_label, levels = reef_order_cots_log))

p_cots_log <- ggplot(combined_cots_log, aes(y = reef_label, x = mean, color = method)) +
  dual_band(dodge_width = NULL) +
  facet_grid(method ~ ., scales = "free_y", space = "free_y") +
  guides(color = "none") +
  scale_x_continuous(breaks = log1p(count_breaks), labels = count_breaks) +
  labs(
    y = "Reef (ordered by reef ID)",
    x = "COTS observed (count per tow or survey, log1p scale)",
    title = "Reef-level COTS counts (log1p scale)"
  ) +
  theme(
    axis.text.y = element_text(size = 5),
    strip.text = element_text(size = 7),
    panel.spacing = unit(3, "pt")
  )

ggsave(
  "output/figures/reef_cots_count_log1p_all_reefs.png",
  p_cots_log, width = 8.27, height = 11.69, dpi = 300, limitsize = FALSE
)

cat("Saved output/figures/reef_cots_count_log1p_all_reefs.png\n")
