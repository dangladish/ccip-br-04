# uncertainty_shared.R
# Helpers for within-reef spread (mean and SD) of coral cover and COTS counts across each
# method's sub-reef sampling units, per reef, pooled across reefs, and per reef-ID band.
# Used by the uncertainty report script and the Shiny data-bundling job.
# Sampling unit: tow (manta, ReefScan Surface; ReefScan Deep COTS), ~100 m grid cell
# (ReefScan Deep coral cover) or survey (RHIS). This describes heterogeneity within reefs,
# not the precision (`_se` columns) of a reef mean.

library(tidyverse)

uncertainty_method_levels <- c("manta", "rs_surface", "rs_deep", "rhis")
uncertainty_method_labels <- c("Manta tow", "ReefScan Surface", "ReefScan Deep", "RHIS")

# ReefScan Deep coral cover has no reliable per-tow grain; a ~0.001 deg (~100 m) grid cell
# is used as the sub-reef unit. Cells with < 5 classified points are dropped.
compute_deep_cell_summary <- function(rs_deep_raw, reef_xwalk) {
  deep_site_to_reef <- reef_xwalk |>
    filter(has_rs_deep) |>
    separate_rows(rs_deep_sites, sep = "; ") |>
    select(rs_deep_site = rs_deep_sites, reef_label)

  rs_deep_raw |>
    filter(!is.na(site_reef_name)) |>
    inner_join(deep_site_to_reef, by = c("site_reef_name" = "rs_deep_site")) |>
    mutate(cell_lat = round(image_latitude, 3), cell_lon = round(image_longitude, 3)) |>
    group_by(reef_label, cell_lat, cell_lon) |>
    summarise(
      n_classified = sum(point_machine_group_code != "IN", na.rm = TRUE),
      pct_HC       = sum(point_machine_group_code == "HC", na.rm = TRUE) / n_classified,
      .groups = "drop"
    ) |>
    filter(n_classified >= 5)  # drop near-empty cells
}

# Read ReefScan Deep raw points and drop imagery whose survey date falls in the date
# window of an excluded manta voyage (the raw file has no tow_id), matching the exclusion
# applied when building the canonical tables.
load_rs_deep_raw_3voyage_excluded <- function(
  path = "data/raw/reefscan_deep/reefcloud-point-classification-GRB-ReefScan-Deep-2026-08-13-filtered.csv",
  manta_path = "data/raw/manta_tows/manta_tows.parquet",
  exclude_voyages = c("FMPCRS 2024-05", "FMPCRS 2025-01", "FMPCRS 2026-03")
) {
  rs_deep_raw <- readr::read_csv(path, show_col_types = FALSE)

  voyage_windows <- arrow::read_parquet(manta_path) |>
    dplyr::filter(!is.na(tow_distance_m), tow_distance_m > 0) |>
    dplyr::summarise(start = min(tow_date), end = max(tow_date), .by = voyage)

  in_excluded_window <- function(d) {
    Reduce(`|`, lapply(exclude_voyages, function(v) {
      w <- voyage_windows[voyage_windows$voyage == v, ]
      !is.na(d) & d >= w$start & d <= w$end
    }))
  }

  n_pre <- nrow(rs_deep_raw)
  rs_deep_raw <- rs_deep_raw[!in_excluded_window(as.Date(rs_deep_raw$`survey_start_date (UTC)`)), ]
  cat("RS Deep raw points (uncertainty layer): ", n_pre, " -> ", nrow(rs_deep_raw),
      " after 3-voyage exclusion (", round(100 * (n_pre - nrow(rs_deep_raw)) / n_pre, 1), "% excluded)\n", sep = "")
  rs_deep_raw
}

# Per-reef mean and SD of hard coral cover across sub-reef units.
compute_coral_uncertainty_reef <- function(tow_level, deep_cell_summ_all, rhis, reef_xwalk) {
  tow_summ_all <- tow_level |>
    group_by(reef_label) |>
    summarise(
      manta_mean      = mean(perc_hard, na.rm = TRUE),
      manta_sd        = sd(perc_hard, na.rm = TRUE),
      manta_n         = sum(!is.na(perc_hard)),
      rs_surface_mean = mean(pct_HC_surf, na.rm = TRUE),
      rs_surface_sd   = sd(pct_HC_surf, na.rm = TRUE),
      rs_surface_n    = sum(!is.na(pct_HC_surf)),
      .groups = "drop"
    ) |>
    mutate(across(c(rs_surface_mean, rs_surface_sd), ~ ifelse(rs_surface_n == 0, NA_real_, .x)))

  deep_summ_all <- deep_cell_summ_all |>
    group_by(reef_label) |>
    summarise(
      rs_deep_mean = mean(pct_HC, na.rm = TRUE),
      rs_deep_sd   = sd(pct_HC, na.rm = TRUE),
      rs_deep_n    = n(),
      .groups = "drop"
    )

  rhis_summ_all <- rhis |>
    filter(reef_label %in% reef_xwalk$reef_label, !is.na(hc_prop_survey)) |>
    group_by(reef_label) |>
    summarise(
      rhis_mean = mean(hc_prop_survey, na.rm = TRUE),
      rhis_sd   = sd(hc_prop_survey, na.rm = TRUE),
      rhis_n    = n(),
      .groups = "drop"
    )

  tow_summ_all |>
    full_join(deep_summ_all, by = "reef_label") |>
    pivot_longer(
      cols = -reef_label,
      names_to = c("method", ".value"),
      names_pattern = "(manta|rs_surface|rs_deep)_(mean|sd|n)"
    ) |>
    filter(!is.na(mean)) |>
    mutate(reef_label = as.character(reef_label)) |>
    bind_rows(
      rhis_summ_all |> transmute(reef_label, method = "rhis", mean = rhis_mean, sd = rhis_sd, n = rhis_n)
    ) |>
    mutate(
      method = factor(method, levels = uncertainty_method_levels, labels = uncertainty_method_labels),
      quantity = "Hard coral cover (proportion)"
    )
}

# Per-reef mean and SD of COTS count across sub-reef units. RHIS uses only surveys with the
# predator module (NA there means not surveyed). log1p = TRUE transforms each unit's count
# before taking mean/SD.
compute_cots_uncertainty_reef <- function(tow_level, rhis, reef_xwalk, log1p = FALSE) {
  tf <- if (log1p) base::log1p else identity

  cots_tow_summ <- tow_level |>
    mutate(
      manta_v      = tf(cots_observed),
      rs_surface_v = tf(cots_reefscan_surface),
      rs_deep_v    = tf(cots_reefscan_deep)
    ) |>
    group_by(reef_label) |>
    summarise(
      manta_mean      = mean(manta_v, na.rm = TRUE),
      manta_sd        = sd(manta_v, na.rm = TRUE),
      manta_n         = sum(!is.na(manta_v)),
      rs_surface_mean = mean(rs_surface_v, na.rm = TRUE),
      rs_surface_sd   = sd(rs_surface_v, na.rm = TRUE),
      rs_surface_n    = sum(!is.na(rs_surface_v)),
      rs_deep_mean    = mean(rs_deep_v, na.rm = TRUE),
      rs_deep_sd      = sd(rs_deep_v, na.rm = TRUE),
      rs_deep_n       = sum(!is.na(rs_deep_v)),
      .groups = "drop"
    )

  rhis_cots_summ <- rhis |>
    filter(reef_label %in% reef_xwalk$reef_label, has_predator_survey) |>
    mutate(rhis_v = tf(cots_count_survey)) |>
    group_by(reef_label) |>
    summarise(
      rhis_mean = mean(rhis_v, na.rm = TRUE),
      rhis_sd   = sd(rhis_v, na.rm = TRUE),
      rhis_n    = n(),
      .groups = "drop"
    )

  cots_tow_summ |>
    pivot_longer(
      cols = -reef_label,
      names_to = c("method", ".value"),
      names_pattern = "(manta|rs_surface|rs_deep)_(mean|sd|n)"
    ) |>
    filter(n > 0) |>
    mutate(reef_label = as.character(reef_label)) |>
    bind_rows(
      rhis_cots_summ |> transmute(reef_label, method = "rhis", mean = rhis_mean, sd = rhis_sd, n = rhis_n)
    ) |>
    mutate(
      method = factor(method, levels = uncertainty_method_levels, labels = uncertainty_method_labels),
      quantity = "COTS count (per tow or survey)"
    )
}

# Pooled across all reefs, one row per method (units pooled directly, not a mean of reef
# means). `rhis_scope` is RHIS restricted to the reefs in reef_xwalk.
compute_coral_uncertainty_pooled <- function(tow_level, deep_cell_summ_all, rhis_scope) {
  bind_rows(
    tow_level |>
      summarise(mean = mean(perc_hard, na.rm = TRUE), sd = sd(perc_hard, na.rm = TRUE), n = sum(!is.na(perc_hard))) |>
      mutate(method = "manta"),
    tow_level |>
      summarise(mean = mean(pct_HC_surf, na.rm = TRUE), sd = sd(pct_HC_surf, na.rm = TRUE), n = sum(!is.na(pct_HC_surf))) |>
      mutate(method = "rs_surface"),
    deep_cell_summ_all |>
      summarise(mean = mean(pct_HC, na.rm = TRUE), sd = sd(pct_HC, na.rm = TRUE), n = n()) |>
      mutate(method = "rs_deep"),
    rhis_scope |>
      filter(!is.na(hc_prop_survey)) |>
      summarise(mean = mean(hc_prop_survey, na.rm = TRUE), sd = sd(hc_prop_survey, na.rm = TRUE), n = n()) |>
      mutate(method = "rhis")
  ) |>
    mutate(
      method = factor(method, levels = uncertainty_method_levels, labels = uncertainty_method_labels),
      quantity = "Hard coral cover (proportion)"
    )
}

compute_cots_uncertainty_pooled <- function(tow_level, rhis_scope, log1p = FALSE) {
  tf <- if (log1p) base::log1p else identity

  bind_rows(
    tow_level |> mutate(v = tf(cots_observed)) |>
      summarise(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), n = sum(!is.na(v))) |>
      mutate(method = "manta"),
    tow_level |> mutate(v = tf(cots_reefscan_surface)) |>
      summarise(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), n = sum(!is.na(v))) |>
      mutate(method = "rs_surface"),
    tow_level |> mutate(v = tf(cots_reefscan_deep)) |>
      summarise(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), n = sum(!is.na(v))) |>
      mutate(method = "rs_deep"),
    rhis_scope |> filter(has_predator_survey) |> mutate(v = tf(cots_count_survey)) |>
      summarise(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), n = n()) |>
      mutate(method = "rhis")
  ) |>
    filter(n > 0) |>
    mutate(
      method = factor(method, levels = uncertainty_method_levels, labels = uncertainty_method_labels),
      quantity = "COTS count (per tow or survey)"
    )
}

# Pooled within reef-ID band: the numeric prefix of reef_label ("21-477" -> "21"), a coarse
# geographic grouping. `band_label` is the display string; `band` is kept for sorting.
reef_band <- function(reef_label) sub("-.*$", "", as.character(reef_label))

compute_coral_uncertainty_band <- function(tow_level, deep_cell_summ_all, rhis_scope) {
  bind_rows(
    tow_level |> mutate(band = reef_band(reef_label)) |> group_by(band) |>
      summarise(mean = mean(perc_hard, na.rm = TRUE), sd = sd(perc_hard, na.rm = TRUE), n = sum(!is.na(perc_hard)), .groups = "drop") |>
      mutate(method = "manta"),
    tow_level |> mutate(band = reef_band(reef_label)) |> group_by(band) |>
      summarise(mean = mean(pct_HC_surf, na.rm = TRUE), sd = sd(pct_HC_surf, na.rm = TRUE), n = sum(!is.na(pct_HC_surf)), .groups = "drop") |>
      mutate(method = "rs_surface"),
    deep_cell_summ_all |> mutate(band = reef_band(reef_label)) |> group_by(band) |>
      summarise(mean = mean(pct_HC, na.rm = TRUE), sd = sd(pct_HC, na.rm = TRUE), n = n(), .groups = "drop") |>
      mutate(method = "rs_deep"),
    rhis_scope |> filter(!is.na(hc_prop_survey)) |> mutate(band = reef_band(reef_label)) |> group_by(band) |>
      summarise(mean = mean(hc_prop_survey, na.rm = TRUE), sd = sd(hc_prop_survey, na.rm = TRUE), n = n(), .groups = "drop") |>
      mutate(method = "rhis")
  ) |>
    filter(n > 0) |>
    mutate(
      method = factor(method, levels = uncertainty_method_levels, labels = uncertainty_method_labels),
      band_label = paste0(band, "-XXX"),
      quantity = "Hard coral cover (proportion)"
    )
}

compute_cots_uncertainty_band <- function(tow_level, rhis_scope, log1p = FALSE) {
  tf <- if (log1p) base::log1p else identity

  bind_rows(
    tow_level |> mutate(v = tf(cots_observed), band = reef_band(reef_label)) |> group_by(band) |>
      summarise(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), n = sum(!is.na(v)), .groups = "drop") |>
      mutate(method = "manta"),
    tow_level |> mutate(v = tf(cots_reefscan_surface), band = reef_band(reef_label)) |> group_by(band) |>
      summarise(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), n = sum(!is.na(v)), .groups = "drop") |>
      mutate(method = "rs_surface"),
    tow_level |> mutate(v = tf(cots_reefscan_deep), band = reef_band(reef_label)) |> group_by(band) |>
      summarise(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), n = sum(!is.na(v)), .groups = "drop") |>
      mutate(method = "rs_deep"),
    rhis_scope |> filter(has_predator_survey) |> mutate(v = tf(cots_count_survey), band = reef_band(reef_label)) |> group_by(band) |>
      summarise(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE), n = n(), .groups = "drop") |>
      mutate(method = "rhis")
  ) |>
    filter(n > 0) |>
    mutate(
      method = factor(method, levels = uncertainty_method_levels, labels = uncertainty_method_labels),
      band_label = paste0(band, "-XXX"),
      quantity = "COTS count (per tow or survey)"
    )
}
