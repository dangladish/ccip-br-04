# manta_surf_cull_cots.R
# Manta <-> RS Surface COTS calibration at cull-site level (follows the D-02 approach).
# Posteriors saved (via the fit_* helpers in R/cots_shared.R):
#   manta_from_surf_combined_cull_cots: manta_combined (tows with scar OR >= 1 COTS)
#     ~ RS Surface density; NegBin, offset log(n_tows), reef random effect
#   manta_from_surf_binary_cull_cots: manta_combined successes of n_tows
#     ~ RS Surface density; Binomial, logit, reef random effect
#   surf_from_manta_cull_cots: RS Surface COTS count ~ manta density
#     ((COTS + scar count) / n_tows, summed); NegBin, offset log(effort_surf), reef random effect
# Requires from the driver: ms_cull (cots_crossmethod_data_prep.R).

post_manta_from_surf_combined_cull_cots <- fit_crossmethod_pair(ms_cull$manta_combined, ms_cull$effort_manta, ms_cull$y_surf, ms_cull$effort_surf, ms_cull$reef_label, 0.01, "manta_from_surf_combined_cull_cots", "cull")

post_manta_from_surf_binary_cull_cots <- fit_binary_crossmethod_cull_reef_pair(ms_cull$manta_combined, ms_cull$effort_manta, ms_cull$y_surf / ms_cull$effort_surf, 0.01, ms_cull$reef_label, "manta_from_surf_binary_cull_cots", "cull")

post_surf_from_manta_cull_cots <- fit_crossmethod_pair(ms_cull$y_surf, ms_cull$effort_surf, ms_cull$y_manta + ms_cull$scar_count, ms_cull$effort_manta, ms_cull$reef_label, 0.01, "surf_from_manta_cull_cots", "cull")

if (doPlot) {
  plot_crossmethod_pair_density(post_manta_from_surf_combined_cull_cots, ms_cull$manta_combined, ms_cull$effort_manta, ms_cull$y_surf, ms_cull$effort_surf, "RS Surface density (COTS/tow)", "Manta combined COTS+scar density (per tow)", "Manta combined COTS+scar density from RS Surface, cull-site level")

  x_seq <- seq(0, max(ms_cull$y_surf / ms_cull$effort_surf), length.out = 100)
  size_rep <- round(median(ms_cull$effort_manta))
  plot_calib_base(predict_binomial_curve(post_manta_from_surf_binary_cull_cots, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   ms_cull$y_surf / ms_cull$effort_surf, ms_cull$manta_combined / ms_cull$effort_manta,
                   "RS Surface density", "Manta combined COTS+scar presence rate", "Manta combined COTS+scar presence rate from RS Surface, cull-site level")

  plot_crossmethod_pair_density(post_surf_from_manta_cull_cots, ms_cull$y_surf, ms_cull$effort_surf, ms_cull$y_manta + ms_cull$scar_count, ms_cull$effort_manta, "Manta density (combined COTS+scar count / tow)", "RS Surface COTS density (per tow)", "RS Surface density from Manta combined COTS+scar, cull-site level")
}
cat("\nmanta_surf_cull_cots.R complete -- 3 posteriors saved, cull-site level.\n")
