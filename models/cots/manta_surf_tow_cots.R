# manta_surf_tow_cots.R
# Manta <-> RS Surface COTS calibration at tow level, reef random effect (D-02 approach).
# Posteriors saved (via the fit_* helpers in R/cots_shared.R):
#   manta_from_surf_tow_cots: has_manta (scar present OR COTS > 0) ~ RS Surface density;
#     Binomial(n = 1), cloglog, offset log(tow distance); D-02 "Manta from X" template
#   surf_from_manta_tow_cots: RS Surface COTS count ~ manta (COTS + scars counted) per tow;
#     NegBin, offset log(effort_surf)
#   manta_from_surf_combined_tow_cots: manta count (COTS + scars counted) ~ RS Surface
#     density; NegBin, offset log(tow_distance_m)
# The scar code (a/p/c) carries no count, so the binary model is the only
# presence model at this level.
# Requires from the driver: ms_tow, ms_tow_bin (cots_crossmethod_data_prep.R).

post_manta_from_surf_tow_cots <- fit_binary_crossmethod_pair(
  ms_tow_bin$has_manta, ms_tow_bin$tow_dist, ms_tow_bin$surf_density, 0.01,
  ms_tow_bin$reef_label, "manta_from_surf_tow_cots")

# scar_count is on ms_tow_bin (different NA filter); join by tow_id, then sum COTS count and scars counted.
ms_tow_scar <- ms_tow |> left_join(ms_tow_bin |> select(tow_id, scar_count), by = "tow_id") |> filter(!is.na(scar_count)) |>
  mutate(manta_combined_tow = y_manta + scar_count)
post_surf_from_manta_tow_cots <- fit_crossmethod_pair(ms_tow_scar$y_surf, ms_tow_scar$effort_surf, ms_tow_scar$manta_combined_tow, ms_tow_scar$effort_manta_tow, ms_tow_scar$reef_label, 0.01, "surf_from_manta_tow_cots", "tow")
post_manta_from_surf_combined_tow_cots <- fit_crossmethod_pair(ms_tow_scar$manta_combined_tow, ms_tow_scar$effort_manta, ms_tow_scar$y_surf, ms_tow_scar$effort_surf, ms_tow_scar$reef_label, 0.01, "manta_from_surf_combined_tow_cots", "tow")

if (doPlot) {
  x_seq <- seq(0, max(ms_tow_bin$surf_density), length.out = 100)
  plot_calib_base(predict_binomial_curve(post_manta_from_surf_tow_cots, x_seq,
                                          offset = mean(log(ms_tow_bin$tow_dist)), eps = 0.01, link = "cloglog"),
                   ms_tow_bin$surf_density, ms_tow_bin$has_manta,
                   "RS Surface density", "P(manta tow has COTS or scar)", "Manta from RS Surface, tow level")
  plot_crossmethod_pair(post_surf_from_manta_tow_cots, ms_tow_scar$y_surf, ms_tow_scar$effort_surf, ms_tow_scar$manta_combined_tow, ms_tow_scar$effort_manta_tow, "Manta density (combined COTS+scar count / tow)", "RS Surface COTS count (per tow)", "RS Surface from Manta combined COTS+scar, tow level")
  plot_crossmethod_pair(post_manta_from_surf_combined_tow_cots, ms_tow_scar$manta_combined_tow, ms_tow_scar$effort_manta, ms_tow_scar$y_surf, ms_tow_scar$effort_surf, "RS Surface density", "Manta: total COTS + scars observed (per tow)", "Manta combined from RS Surface, tow level")
}
cat("\nmanta_surf_tow_cots.R complete -- 3 posteriors saved, tow level.\n")
