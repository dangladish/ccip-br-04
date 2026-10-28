# manta_deep_cull.R
# Coral cover calibration: manta tow <-> RS Deep, cull-site level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression, logit link, fit in JAGS via fit_jags_beta_pair().
# No reef random effect: too few reefs to estimate one.
# Inputs: cull_site_level (perc_hard, pct_HC_deep).
# Saves: output/posteriors/{manta_from_deep_cull,deep_from_manta_cull}.RDS

calData <- cull_site_level |>
  select(cull_site, reef_label, y_manta = perc_hard, y_deep = pct_HC_deep) |>
  filter(!is.na(y_manta), !is.na(y_deep))

post_manta_from_deep_cull <- fit_jags_beta_pair(calData$y_manta, calData$y_deep, "manta_from_deep_cull")
post_deep_from_manta_cull <- fit_jags_beta_pair(calData$y_deep, calData$y_manta, "deep_from_manta_cull")

if (doPlot) {
  x_seq <- seq(min(calData$y_deep), max(calData$y_deep), length.out = 100)
  plot_calib_base(predict_beta_curve(post_manta_from_deep_cull, x_seq), calData$y_deep, calData$y_manta,
                   "RS Deep HC %", "Manta HC %", "Manta from RS Deep, cull-site level")
  x_seq <- seq(min(calData$y_manta), max(calData$y_manta), length.out = 100)
  plot_calib_base(predict_beta_curve(post_deep_from_manta_cull, x_seq), calData$y_manta, calData$y_deep,
                   "Manta HC %", "RS Deep HC %", "RS Deep from Manta, cull-site level")
}
