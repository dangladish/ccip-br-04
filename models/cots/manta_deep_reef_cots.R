# manta_deep_reef_cots.R
# Manta <-> RS Deep COTS calibration at reef level (follows the D-02 approach).
# Posteriors saved (via the fit_* helpers in R/cots_shared.R):
#   manta_from_deep_combined_reef_cots: manta_combined ~ RS Deep density;
#     NegBin, offset log(n_tows), no reef random effect
#   manta_from_deep_binary_reef_cots: manta_combined successes of n_tows
#     ~ RS Deep density; Binomial, logit, observation-level reef random effect
#   deep_from_manta_reef_cots: RS Deep COTS count ~ manta density
#     ((COTS + scar count) / n_tows, summed); NegBin, offset log(effort_deep), no reef random effect
# No random effect in the NegBin models: one row per reef, so it cannot be
# separated from the NegBin dispersion.
# Requires from the driver: md_reef (cots_crossmethod_data_prep.R).

post_manta_from_deep_combined_reef_cots <- fit_crossmethod_pair(md_reef$manta_combined, md_reef$effort_manta, md_reef$y_deep, md_reef$effort_deep, NULL, 0.01, "manta_from_deep_combined_reef_cots", "reef")

post_manta_from_deep_binary_reef_cots <- fit_binary_crossmethod_cull_reef_pair(md_reef$manta_combined, md_reef$effort_manta, md_reef$y_deep / md_reef$effort_deep, 0.01, md_reef$reef_label, "manta_from_deep_binary_reef_cots", "reef")

post_deep_from_manta_reef_cots <- fit_crossmethod_pair(md_reef$y_deep, md_reef$effort_deep, md_reef$y_manta + md_reef$feeding_scars, md_reef$effort_manta, NULL, 0.01, "deep_from_manta_reef_cots", "reef")

if (doPlot) {
  plot_crossmethod_pair_density(post_manta_from_deep_combined_reef_cots, md_reef$manta_combined, md_reef$effort_manta, md_reef$y_deep, md_reef$effort_deep, "RS Deep density (COTS/tow)", "Manta combined COTS+scar density (per tow)", "Manta combined COTS+scar density from RS Deep, reef level")

  x_seq <- seq(0, max(md_reef$y_deep / md_reef$effort_deep), length.out = 100)
  size_rep <- round(median(md_reef$effort_manta))
  plot_calib_base(predict_binomial_curve(post_manta_from_deep_binary_reef_cots, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   md_reef$y_deep / md_reef$effort_deep, md_reef$manta_combined / md_reef$effort_manta,
                   "RS Deep density", "Manta combined COTS+scar presence rate", "Manta combined COTS+scar presence rate from RS Deep, reef level")

  plot_crossmethod_pair_density(post_deep_from_manta_reef_cots, md_reef$y_deep, md_reef$effort_deep, md_reef$y_manta + md_reef$feeding_scars, md_reef$effort_manta, "Manta density (combined COTS+scar count / tow)", "RS Deep COTS density (per tow)", "RS Deep density from Manta combined COTS+scar, reef level")
}
cat("\nmanta_deep_reef_cots.R complete -- 3 posteriors saved, reef level.\n")
