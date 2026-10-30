# edna_calibration_data_prep.R
# Data prep for the eDNA calibration models (reef level only; eDNA samples are not tied to
# tows, so all samples are pooled to the reef). Sourced WITHOUT local = TRUE: it builds the
# tables edna_reef, reef_manta and reef_surf that the eDNA model files read.
# Inputs (from the driver): edna_raw, xwalk (eDNA reef -> reef_label), reef_level.
#
# Response definitions (pooled over all samples in a reef):
# - Binary: n_success / n_trials = LOD-positive samples / samples with a non-NA LOD call.
# - Concentration (D-02 aggregation, at reef grain): PCR replicate concentrations are summed
#   within a sample (conc_row), then over the reef's samples (conc_reef_sum; n_assays_reef =
#   valid replicates), and rounded to an integer (integer_conc_reef). The count response suits
#   a NegBin with offset log(n_assays_reef). edna_conc_rate_reef = integer_conc_reef /
#   n_assays_reef is kept as a rate for use as a predictor.

edna_reef <- edna_raw |>
  mutate(
    n_valid_conc = (!is.na(pcr1_conc)) + (!is.na(pcr2_conc)),
    conc_row     = coalesce(pcr1_conc, 0) + coalesce(pcr2_conc, 0)
  ) |>
  left_join(xwalk, by = "edna_reef") |>
  filter(!is.na(reef_label)) |>
  group_by(reef_label) |>
  summarise(
    n_samples     = n(),
    n_success     = sum(lod_sample_positive, na.rm = TRUE),
    n_trials      = sum(!is.na(lod_sample_positive)),
    n_assays_reef = sum(n_valid_conc),
    conc_reef_sum = sum(conc_row, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(n_trials > 0, n_assays_reef > 0) |>
  mutate(
    edna_prop_reef      = n_success / n_trials,
    integer_conc_reef   = round(conc_reef_sum, digits = 0),
    edna_conc_rate_reef = integer_conc_reef / n_assays_reef
  )

# ---- Manta, reef grain ----
# manta_density_raw = (COTS observed + tows with scar) / n_tows; manta_density rescales it to
# COTS per ha. Used for eDNA-as-response models; manta-as-response models use the combined
# presence indicator instead.
manta_cots_per_ha_factor <- 68

reef_manta <- edna_reef |>
  inner_join(
    reef_level |> select(reef_label, cots_observed, feeding_scars, scar_count, n_tows, scar_rate, cots_rate,
                         manta_combined, manta_combined_rate),
    by = "reef_label"
  ) |>
  filter(n_tows > 0) |>
  mutate(
    manta_count       = cots_observed + scar_count,
    manta_density_raw = manta_count / n_tows,
    manta_density      = manta_density_raw * manta_cots_per_ha_factor
  )

# ---- RS Surface, reef grain ----
# surf_density = COTS per tow (cots_reefscan_surface / n_tows), as in the cross-method models.
# n_pts_surf_reef is only an imagery-coverage filter.
reef_surf <- edna_reef |>
  inner_join(
    reef_level |> select(reef_label, cots_reefscan_surface, n_pts_surf_reef, n_tows),
    by = "reef_label"
  ) |>
  # Also require a non-NA count: a reef can have RS Surface points but no matched manta tow,
  # leaving cots_reefscan_surface NA.
  filter(n_pts_surf_reef > 0, n_tows > 0, !is.na(cots_reefscan_surface)) |>
  mutate(surf_density = cots_reefscan_surface / n_tows)

cat("eDNA (reef level): edna_reef n =", nrow(edna_reef), "reefs\n")
cat("eDNA: reef_manta n =", nrow(reef_manta), "reefs\n")
cat("eDNA: reef_surf n =", nrow(reef_surf), "reefs\n")
