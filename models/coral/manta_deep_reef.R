# manta_deep_reef.R
# Coral cover calibration: manta tow <-> RS Deep, reef level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression, logit link, fit in JAGS via fit_jags_beta_pair(); one row per reef, no random effect.
# Inputs: reef_level (mean_HC_manta, mean_HC_deep).
# Saves: output/posteriors/{manta_from_deep_reef,deep_from_manta_reef}.RDS

calData <- reef_level |> filter(has_deep, !is.na(mean_HC_manta), !is.na(mean_HC_deep))

post_manta_from_deep_reef <- fit_jags_beta_pair(calData$mean_HC_manta, calData$mean_HC_deep, "manta_from_deep_reef")
post_deep_from_manta_reef <- fit_jags_beta_pair(calData$mean_HC_deep, calData$mean_HC_manta, "deep_from_manta_reef")

if (doPlot) {
  x_seq <- seq(min(calData$mean_HC_deep), max(calData$mean_HC_deep), length.out = 100)
  plot_calib_base(predict_beta_curve(post_manta_from_deep_reef, x_seq), calData$mean_HC_deep, calData$mean_HC_manta,
                   "RS Deep HC %", "Manta HC %", "Manta from RS Deep, reef level")
  x_seq <- seq(min(calData$mean_HC_manta), max(calData$mean_HC_manta), length.out = 100)
  plot_calib_base(predict_beta_curve(post_deep_from_manta_reef, x_seq), calData$mean_HC_manta, calData$mean_HC_deep,
                   "Manta HC %", "RS Deep HC %", "RS Deep from Manta, reef level")
}
