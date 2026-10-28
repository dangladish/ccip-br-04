# manta_deep_tow_cots.R
# Manta <-> RS Deep COTS calibration at tow level, reef random effect.
# Same structure as manta_surf_tow_cots.R with RS Deep in place of RS Surface.
# Posteriors saved (via the fit_* helpers in R/cots_shared.R):
#   manta_from_deep_tow_cots: has_manta (scar present OR COTS > 0) ~ RS Deep density;
#     Binomial(n = 1), cloglog, offset log(tow distance)
#   manta_from_deep_combined_tow_cots: manta count (COTS + scar presence) ~ RS Deep
#     density; NegBin, offset log(tow_distance_m)
#   deep_from_manta_tow_cots: RS Deep COTS count ~ manta (COTS + scar presence) per tow;
#     NegBin, offset log(effort_deep)
# Requires from the driver: md_tow, md_tow_bin (cots_crossmethod_data_prep.R).

post_manta_from_deep_tow_cots <- fit_binary_crossmethod_pair(
  md_tow_bin$has_manta, md_tow_bin$tow_dist, md_tow_bin$deep_density, 0.01,
  md_tow_bin$reef_label, "manta_from_deep_tow_cots")

md_tow_scar <- md_tow |> left_join(md_tow_bin |> select(tow_id, has_scar), by = "tow_id") |> filter(!is.na(has_scar)) |>
  mutate(manta_combined_tow = y_manta + has_scar)
post_manta_from_deep_combined_tow_cots <- fit_crossmethod_pair(md_tow_scar$manta_combined_tow, md_tow_scar$effort_manta, md_tow_scar$y_deep, md_tow_scar$effort_deep, md_tow_scar$reef_label, 0.01, "manta_from_deep_combined_tow_cots", "tow")
post_deep_from_manta_tow_cots <- fit_crossmethod_pair(md_tow_scar$y_deep, md_tow_scar$effort_deep, md_tow_scar$manta_combined_tow, md_tow_scar$effort_manta_tow, md_tow_scar$reef_label, 0.01, "deep_from_manta_tow_cots", "tow")

if (doPlot) {
  x_seq <- seq(0, max(md_tow_bin$deep_density, na.rm = TRUE), length.out = 100)
  plot_calib_base(predict_binomial_curve(post_manta_from_deep_tow_cots, x_seq,
                                          offset = mean(log(md_tow_bin$tow_dist)), eps = 0.01, link = "cloglog"),
                   md_tow_bin$deep_density, md_tow_bin$has_manta,
                   "RS Deep density", "P(manta tow has COTS or scar)", "Manta from RS Deep, tow level (exploratory)")
  plot_crossmethod_pair(post_manta_from_deep_combined_tow_cots, md_tow_scar$manta_combined_tow, md_tow_scar$effort_manta, md_tow_scar$y_deep, md_tow_scar$effort_deep, "RS Deep density", "Manta: total COTS + scars observed (per tow)", "Manta combined from RS Deep, tow level (exploratory)")
  plot_crossmethod_pair(post_deep_from_manta_tow_cots, md_tow_scar$y_deep, md_tow_scar$effort_deep, md_tow_scar$manta_combined_tow, md_tow_scar$effort_manta_tow, "Manta density (combined COTS+scar count / tow)", "RS Deep COTS count (per tow)", "RS Deep from Manta combined COTS+scar, tow level (exploratory)")
}
cat("\nmanta_deep_tow_cots.R complete -- 3 posteriors saved, tow level (exploratory).\n")
