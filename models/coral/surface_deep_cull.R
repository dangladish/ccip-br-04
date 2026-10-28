# surface_deep_cull.R
# Coral cover calibration: RS Surface <-> RS Deep, cull-site level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression, logit link, fit in JAGS via fit_jags_beta_pair().
# No reef random effect: too few reefs to estimate one.
# Inputs: cull_site_level (pct_HC_surf, pct_HC_deep).
# Saves: output/posteriors/{surf_from_deep_cull,deep_from_surf_cull}.RDS

calData <- cull_site_level |>
  select(cull_site, reef_label, y_surf = pct_HC_surf, y_deep = pct_HC_deep) |>
  filter(!is.na(y_surf), !is.na(y_deep))

post_surf_from_deep_cull <- fit_jags_beta_pair(calData$y_surf, calData$y_deep, "surf_from_deep_cull")
post_deep_from_surf_cull <- fit_jags_beta_pair(calData$y_deep, calData$y_surf, "deep_from_surf_cull")

if (doPlot) {
  x_seq <- seq(min(calData$y_deep), max(calData$y_deep), length.out = 100)
  plot_calib_base(predict_beta_curve(post_surf_from_deep_cull, x_seq), calData$y_deep, calData$y_surf,
                   "RS Deep HC %", "RS Surface HC %", "RS Surface from RS Deep, cull-site level")
  x_seq <- seq(min(calData$y_surf), max(calData$y_surf), length.out = 100)
  plot_calib_base(predict_beta_curve(post_deep_from_surf_cull, x_seq), calData$y_surf, calData$y_deep,
                   "RS Surface HC %", "RS Deep HC %", "RS Deep from RS Surface, cull-site level")
}
