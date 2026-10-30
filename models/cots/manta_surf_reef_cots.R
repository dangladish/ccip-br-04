# manta_surf_reef_cots.R
# Manta <-> RS Surface COTS calibration at reef level (follows the D-02 approach).
# Posteriors saved (via the fit_* helpers in R/cots_shared.R):
#   manta_from_surf_combined_reef_cots: manta_combined ~ RS Surface density;
#     NegBin, offset log(n_tows), no reef random effect
#   manta_from_surf_binary_reef_cots: manta_combined successes of n_tows
#     ~ RS Surface density; Binomial, logit, observation-level reef random effect
#   surf_from_manta_reef_cots: RS Surface COTS count ~ manta density
#     ((COTS + scar count) / n_tows, summed); NegBin, offset log(effort_surf), no reef random effect
# No random effect in the NegBin models: one row per reef, so it cannot be
# separated from the NegBin dispersion.
# Requires from the driver: ms_reef (cots_crossmethod_data_prep.R).

post_manta_from_surf_combined_reef_cots <- fit_crossmethod_pair(ms_reef$manta_combined, ms_reef$effort_manta, ms_reef$y_surf, ms_reef$effort_surf, NULL, 0.01, "manta_from_surf_combined_reef_cots", "reef")

post_manta_from_surf_binary_reef_cots <- fit_binary_crossmethod_cull_reef_pair(ms_reef$manta_combined, ms_reef$effort_manta, ms_reef$y_surf / ms_reef$effort_surf, 0.01, ms_reef$reef_label, "manta_from_surf_binary_reef_cots", "reef")

post_surf_from_manta_reef_cots <- fit_crossmethod_pair(ms_reef$y_surf, ms_reef$effort_surf, ms_reef$y_manta + ms_reef$scar_count, ms_reef$effort_manta, NULL, 0.01, "surf_from_manta_reef_cots", "reef")

if (doPlot) {
  plot_crossmethod_pair_density(post_manta_from_surf_combined_reef_cots, ms_reef$manta_combined, ms_reef$effort_manta, ms_reef$y_surf, ms_reef$effort_surf, "RS Surface density (COTS/tow)", "Manta combined COTS+scar density (per tow)", "Manta combined COTS+scar density from RS Surface, reef level")

  x_seq <- seq(0, max(ms_reef$y_surf / ms_reef$effort_surf), length.out = 100)
  size_rep <- round(median(ms_reef$effort_manta))
  plot_calib_base(predict_binomial_curve(post_manta_from_surf_binary_reef_cots, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   ms_reef$y_surf / ms_reef$effort_surf, ms_reef$manta_combined / ms_reef$effort_manta,
                   "RS Surface density", "Manta combined COTS+scar presence rate", "Manta combined COTS+scar presence rate from RS Surface, reef level")

  plot_crossmethod_pair_density(post_surf_from_manta_reef_cots, ms_reef$y_surf, ms_reef$effort_surf, ms_reef$y_manta + ms_reef$scar_count, ms_reef$effort_manta, "Manta density (combined COTS+scar count / tow)", "RS Surface COTS density (per tow)", "RS Surface density from Manta combined COTS+scar, reef level")
}
cat("\nmanta_surf_reef_cots.R complete -- 3 posteriors saved, reef level.\n")
