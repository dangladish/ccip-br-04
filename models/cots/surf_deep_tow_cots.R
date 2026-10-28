# surf_deep_tow_cots.R
# RS Surface <-> RS Deep COTS count calibration at tow level.
# Posteriors saved: surf_from_deep_tow_cots, deep_from_surf_tow_cots
# (NegBin, offset log(effort of the response method), reef random effect).
# Requires from the driver: sd_tow (cots_crossmethod_data_prep.R).

post_surf_from_deep_tow_cots <- fit_crossmethod_pair(sd_tow$y_surf, sd_tow$effort_surf, sd_tow$y_deep, sd_tow$effort_deep, sd_tow$reef_label, 0.01, "surf_from_deep_tow_cots", "tow")
post_deep_from_surf_tow_cots <- fit_crossmethod_pair(sd_tow$y_deep, sd_tow$effort_deep, sd_tow$y_surf, sd_tow$effort_surf, sd_tow$reef_label, 0.01, "deep_from_surf_tow_cots", "tow")

if (doPlot) {
  plot_crossmethod_pair(post_surf_from_deep_tow_cots, sd_tow$y_surf, sd_tow$effort_surf, sd_tow$y_deep, sd_tow$effort_deep, "RS Deep density", "RS Surface COTS count (per tow)", "RS Surface from RS Deep, tow level (exploratory)")
  plot_crossmethod_pair(post_deep_from_surf_tow_cots, sd_tow$y_deep, sd_tow$effort_deep, sd_tow$y_surf, sd_tow$effort_surf, "RS Surface density", "RS Deep COTS count (per tow)", "RS Deep from RS Surface, tow level (exploratory)")
}
cat("\nsurf_deep_tow_cots.R complete -- 2 posteriors saved, tow level (exploratory).\n")
