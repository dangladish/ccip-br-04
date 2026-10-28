# cots_crossmethod_data_prep.R
# Data prep for the cross-method COTS count models (manta, RS Surface, RS Deep) at tow, cull-site
# and reef level. Sourced WITHOUT local = TRUE: it builds the tables the model files read.
# Inputs (from the driver): tow_level, cull_site_level, reef_level.
# Outputs: ms_* (manta-surface), md_* (manta-deep), sd_* (surface-deep) tables at tow/cull/reef
#   level, plus ms_tow_bin / md_tow_bin (binary manta presence, tow level).
# Models use the other method's density as predictor, log(density + 0.01), and the response's
# own effort as the offset.
#
# Effort conventions:
# - RS Surface/Deep effort is in tows: 1 per row at tow level, and manta's n_tows at cull-site
#   and reef level (an approximation of their own effort), so densities are COTS per tow.
#   n_pts_* columns are only used to require that imagery exists.
# - Manta effort at tow level has two roles. effort_manta (= tow_distance_m) is the detection
#   offset when manta is the response. effort_manta_tow (= 1) is used when manta is the
#   predictor, so its density is per tow rather than per metre towed.
# - Tow-level manta-as-response models follow the D-02 binary template (presence of scar or
#   COTS, offset = log tow distance); cull-site/reef manta-as-response models have no D-02
#   template.

# ---- Tow level: count tables ----

ms_tow <- tow_level |> filter(!is.na(cots_observed), !is.na(cots_reefscan_surface), !is.na(tow_distance_m), !is.na(n_pts_surf)) |>
  transmute(tow_id, reef_label, y_manta = cots_observed, effort_manta = tow_distance_m, effort_manta_tow = 1, y_surf = cots_reefscan_surface, effort_surf = 1)
md_tow <- tow_level |> filter(!is.na(cots_observed), !is.na(cots_reefscan_deep), !is.na(tow_distance_m), !is.na(n_pts_deep)) |>
  transmute(tow_id, reef_label, y_manta = cots_observed, effort_manta = tow_distance_m, effort_manta_tow = 1, y_deep = cots_reefscan_deep, effort_deep = 1)
sd_tow <- tow_level |> filter(!is.na(cots_reefscan_surface), !is.na(cots_reefscan_deep), !is.na(n_pts_surf), !is.na(n_pts_deep)) |>
  transmute(tow_id, reef_label, y_surf = cots_reefscan_surface, effort_surf = 1, y_deep = cots_reefscan_deep, effort_deep = 1)

# ---- Tow level: binary tables (manta as response) ----
# has_manta = scar present OR >= 1 COTS observed. has_scar is kept because the manta-as-predictor
# models sum scar presence with the COTS count.
ms_tow_bin <- tow_level |>
  filter(!is.na(scar_present), !is.na(cots_reefscan_surface), !is.na(tow_distance_m),
         tow_distance_m > 0, !is.na(n_pts_surf)) |>
  transmute(tow_id, reef_label, has_manta = as.integer(scar_present == 1 | cots_observed > 0),
            has_scar = scar_present, tow_dist = tow_distance_m,
            surf_density = cots_reefscan_surface)
md_tow_bin <- tow_level |>
  filter(!is.na(scar_present), !is.na(cots_reefscan_deep), !is.na(tow_distance_m),
         tow_distance_m > 0, !is.na(n_pts_deep)) |>
  transmute(tow_id, reef_label, has_manta = as.integer(scar_present == 1 | cots_observed > 0),
            has_scar = scar_present, tow_dist = tow_distance_m,
            deep_density = cots_reefscan_deep)

# ---- Cull-site and reef tables ----
# Manta effort is n_tows here (a tow count); tow level uses tow_distance_m.
# feeding_scars = number of tows with scars present (manta records only a/p/c presence, not a
# tally); needed only where manta is involved.
# manta_combined / manta_combined_rate = tows with scar OR COTS observed (response). The manta
# predictor is (y_manta + feeding_scars) / effort_manta.
ms_cull <- cull_site_level |> filter(!is.na(cots_observed), !is.na(cots_reefscan_surface), !is.na(n_tows), !is.na(n_pts_surf)) |>
  transmute(cull_site, reef_label, y_manta = cots_observed, effort_manta = n_tows, y_surf = cots_reefscan_surface, effort_surf = n_tows, feeding_scars, manta_combined, manta_combined_rate)
md_cull <- cull_site_level |> filter(!is.na(cots_observed), !is.na(cots_reefscan_deep), !is.na(n_tows), !is.na(n_pts_deep)) |>
  transmute(cull_site, reef_label, y_manta = cots_observed, effort_manta = n_tows, y_deep = cots_reefscan_deep, effort_deep = n_tows, feeding_scars, manta_combined, manta_combined_rate)
sd_cull <- cull_site_level |> filter(!is.na(cots_reefscan_surface), !is.na(cots_reefscan_deep), !is.na(n_pts_surf), !is.na(n_pts_deep), !is.na(n_tows), n_tows > 0) |>
  transmute(cull_site, reef_label, y_surf = cots_reefscan_surface, effort_surf = n_tows, y_deep = cots_reefscan_deep, effort_deep = n_tows)

ms_reef <- reef_level |> filter(!is.na(cots_observed), !is.na(cots_reefscan_surface), !is.na(n_tows), !is.na(n_pts_surf_reef)) |>
  transmute(reef_label, y_manta = cots_observed, effort_manta = n_tows, y_surf = cots_reefscan_surface, effort_surf = n_tows, feeding_scars, manta_combined, manta_combined_rate)
md_reef <- reef_level |> filter(!is.na(cots_observed), !is.na(cots_reefscan_deep), !is.na(n_tows), !is.na(n_pts_deep_reef)) |>
  transmute(reef_label, y_manta = cots_observed, effort_manta = n_tows, y_deep = cots_reefscan_deep, effort_deep = n_tows, feeding_scars, manta_combined, manta_combined_rate)
sd_reef <- reef_level |> filter(!is.na(cots_reefscan_surface), !is.na(cots_reefscan_deep), !is.na(n_pts_surf_reef), !is.na(n_pts_deep_reef), !is.na(n_tows), n_tows > 0) |>
  transmute(reef_label, y_surf = cots_reefscan_surface, effort_surf = n_tows, y_deep = cots_reefscan_deep, effort_deep = n_tows)

cat("tow  ms n=", nrow(ms_tow), " md n=", nrow(md_tow), " sd n=", nrow(sd_tow), "\n")
cat("tow (binary, manta-as-response) ms n=", nrow(ms_tow_bin), " md n=", nrow(md_tow_bin), "\n")
cat("cull ms n=", nrow(ms_cull), " md n=", nrow(md_cull), " sd n=", nrow(sd_cull), "\n")
cat("reef ms n=", nrow(ms_reef), " md n=", nrow(md_reef), " sd n=", nrow(sd_reef), "\n")
