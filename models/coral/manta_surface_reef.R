# manta_surface_reef.R
# Coral cover calibration: manta tow <-> RS Surface, reef level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression, fit in JAGS via fit_jags_beta_pair(); one row per reef, no random effect.
# manta_from_surf uses a loglog link (written directly into the JAGS model, so alpha/beta are
# already on the loglog(y) scale); surf_from_manta uses logit.
# Inputs: reef_level (mean_HC_manta, mean_HC_surf).
# Saves: output/posteriors/{manta_from_surf_reef,surf_from_manta_reef}.RDS

calData <- reef_level |> filter(has_surf, !is.na(mean_HC_manta), !is.na(mean_HC_surf))

post_manta_from_surf_reef <- fit_jags_beta_pair(calData$mean_HC_manta, calData$mean_HC_surf, "manta_from_surf_reef", link = "loglog")
post_surf_from_manta_reef <- fit_jags_beta_pair(calData$mean_HC_surf, calData$mean_HC_manta, "surf_from_manta_reef")

if (doPlot) {
  x_seq <- seq(min(calData$mean_HC_surf), max(calData$mean_HC_surf), length.out = 100)
  plot_calib_base(predict_beta_curve(post_manta_from_surf_reef, x_seq), calData$mean_HC_surf, calData$mean_HC_manta,
                   "RS Surface HC %", "Manta HC %", "Manta from RS Surface, reef level")
  x_seq <- seq(min(calData$mean_HC_manta), max(calData$mean_HC_manta), length.out = 100)
  plot_calib_base(predict_beta_curve(post_surf_from_manta_reef, x_seq), calData$mean_HC_manta, calData$mean_HC_surf,
                   "Manta HC %", "RS Surface HC %", "RS Surface from Manta, reef level")
}
