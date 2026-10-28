# surface_deep_reef.R
# Coral cover calibration: RS Surface <-> RS Deep, reef level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression, logit link, fit in JAGS via fit_jags_beta_pair(); one row per reef, no random effect.
# Inputs: reef_level (mean_HC_surf, mean_HC_deep).
# Saves: output/posteriors/{surf_from_deep_reef,deep_from_surf_reef}.RDS

calData <- reef_level |> filter(has_surf, has_deep, !is.na(mean_HC_surf), !is.na(mean_HC_deep))

post_surf_from_deep_reef <- fit_jags_beta_pair(calData$mean_HC_surf, calData$mean_HC_deep, "surf_from_deep_reef")
post_deep_from_surf_reef <- fit_jags_beta_pair(calData$mean_HC_deep, calData$mean_HC_surf, "deep_from_surf_reef")

if (doPlot) {
  x_seq <- seq(min(calData$mean_HC_deep), max(calData$mean_HC_deep), length.out = 100)
  plot_calib_base(predict_beta_curve(post_surf_from_deep_reef, x_seq), calData$mean_HC_deep, calData$mean_HC_surf,
                   "RS Deep HC %", "RS Surface HC %", "RS Surface from RS Deep, reef level")
  x_seq <- seq(min(calData$mean_HC_surf), max(calData$mean_HC_surf), length.out = 100)
  plot_calib_base(predict_beta_curve(post_deep_from_surf_reef, x_seq), calData$mean_HC_surf, calData$mean_HC_deep,
                   "RS Surface HC %", "RS Deep HC %", "RS Deep from RS Surface, reef level")
}
