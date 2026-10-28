# manta_deep_cull_cots.R
# Manta <-> RS Deep COTS calibration at cull-site level (follows the D-02 approach).
# Posteriors saved (via the fit_* helpers in R/cots_shared.R):
#   manta_from_deep_combined_cull_cots: manta_combined (tows with scar OR >= 1 COTS)
#     ~ RS Deep density; NegBin, offset log(n_tows), reef random effect
#   manta_from_deep_binary_cull_cots: manta_combined successes of n_tows
#     ~ RS Deep density; Binomial, logit, reef random effect
#   deep_from_manta_cull_cots: RS Deep COTS count ~ manta density
#     ((COTS + scar count) / n_tows, summed); NegBin, offset log(effort_deep), reef random effect
# RS Deep is only modelled at cull-site/reef level here because of poor tow-level alignment.
# Requires from the driver: md_cull (cots_crossmethod_data_prep.R).

post_manta_from_deep_combined_cull_cots <- fit_crossmethod_pair(md_cull$manta_combined, md_cull$effort_manta, md_cull$y_deep, md_cull$effort_deep, md_cull$reef_label, 0.01, "manta_from_deep_combined_cull_cots", "cull")

post_manta_from_deep_binary_cull_cots <- fit_binary_crossmethod_cull_reef_pair(md_cull$manta_combined, md_cull$effort_manta, md_cull$y_deep / md_cull$effort_deep, 0.01, md_cull$reef_label, "manta_from_deep_binary_cull_cots", "cull")

post_deep_from_manta_cull_cots <- fit_crossmethod_pair(md_cull$y_deep, md_cull$effort_deep, md_cull$y_manta + md_cull$feeding_scars, md_cull$effort_manta, md_cull$reef_label, 0.01, "deep_from_manta_cull_cots", "cull")

if (doPlot) {
  plot_crossmethod_pair_density(post_manta_from_deep_combined_cull_cots, md_cull$manta_combined, md_cull$effort_manta, md_cull$y_deep, md_cull$effort_deep, "RS Deep density (COTS/tow)", "Manta combined COTS+scar density (per tow)", "Manta combined COTS+scar density from RS Deep, cull-site level")

  x_seq <- seq(0, max(md_cull$y_deep / md_cull$effort_deep), length.out = 100)
  size_rep <- round(median(md_cull$effort_manta))
  plot_calib_base(predict_binomial_curve(post_manta_from_deep_binary_cull_cots, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   md_cull$y_deep / md_cull$effort_deep, md_cull$manta_combined / md_cull$effort_manta,
                   "RS Deep density", "Manta combined COTS+scar presence rate", "Manta combined COTS+scar presence rate from RS Deep, cull-site level")

  plot_crossmethod_pair_density(post_deep_from_manta_cull_cots, md_cull$y_deep, md_cull$effort_deep, md_cull$y_manta + md_cull$feeding_scars, md_cull$effort_manta, "Manta density (combined COTS+scar count / tow)", "RS Deep COTS density (per tow)", "RS Deep density from Manta combined COTS+scar, cull-site level")
}
cat("\nmanta_deep_cull_cots.R complete -- 3 posteriors saved, cull-site level.\n")
