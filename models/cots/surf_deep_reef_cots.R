# surf_deep_reef_cots.R
# RS Surface <-> RS Deep COTS count calibration at reef level.
# Posteriors saved: surf_from_deep_reef_cots, deep_from_surf_reef_cots
# (NegBin, offset log(effort of the response method), no reef random effect:
# one row per reef).
# Requires from the driver: sd_reef (cots_crossmethod_data_prep.R).

post_surf_from_deep_reef_cots <- fit_crossmethod_pair(sd_reef$y_surf, sd_reef$effort_surf, sd_reef$y_deep, sd_reef$effort_deep, NULL, 0.01, "surf_from_deep_reef_cots", "reef")
post_deep_from_surf_reef_cots <- fit_crossmethod_pair(sd_reef$y_deep, sd_reef$effort_deep, sd_reef$y_surf, sd_reef$effort_surf, NULL, 0.01, "deep_from_surf_reef_cots", "reef")

if (doPlot) {
  plot_crossmethod_pair_density(post_surf_from_deep_reef_cots, sd_reef$y_surf, sd_reef$effort_surf, sd_reef$y_deep, sd_reef$effort_deep, "RS Deep density (COTS/tow)", "RS Surface COTS density (per tow)", "RS Surface density from RS Deep, reef level")
  plot_crossmethod_pair_density(post_deep_from_surf_reef_cots, sd_reef$y_deep, sd_reef$effort_deep, sd_reef$y_surf, sd_reef$effort_surf, "RS Surface density (COTS/tow)", "RS Deep COTS density (per tow)", "RS Deep density from RS Surface, reef level")
}
cat("\nCross-method driver complete.\n")
