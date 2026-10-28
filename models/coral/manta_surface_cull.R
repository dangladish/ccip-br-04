# manta_surface_cull.R
# Coral cover calibration: manta tow <-> RS Surface, cull-site level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression (brms/Stan), reef random effect (1 | reef_label), D-02 prior on phi.
# manta_from_surf uses a loglog link (see manta_surface_tow.R); surf_from_manta uses logit.
# Inputs: cull_site_level (perc_hard, pct_HC_surf).
# Saves: output/posteriors/{manta_from_surf_cull,surf_from_manta_cull}.RDS

calData <- cull_site_level |>
  select(cull_site, reef_label, y_manta = perc_hard, y_surf = pct_HC_surf) |>
  filter(!is.na(y_manta), !is.na(y_surf))

fit_manta_from_surf_cull <- brm(
  bf(clip01(1 - y_manta) ~ y_surf + (1 | reef_label)), family = Beta(link = "cloglog"),
  data = calData, prior = c(prior(normal(0, 1), class = "sd"), d02_prior_coral_phi),
  chains = 4, iter = 4000, warmup = 1000, seed = coral_seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_manta_from_surf_cull <- as.matrix(as_draws_matrix(fit_manta_from_surf_cull))
colnames(post_manta_from_surf_cull)[colnames(post_manta_from_surf_cull) == "b_Intercept"] <- "alpha"
colnames(post_manta_from_surf_cull)[colnames(post_manta_from_surf_cull) == "b_y_surf"]    <- "beta"
# Fitted as cloglog on 1 - y; flip signs to get the loglog(y) scale.
post_manta_from_surf_cull[, "alpha"] <- -post_manta_from_surf_cull[, "alpha"]
post_manta_from_surf_cull[, "beta"]  <- -post_manta_from_surf_cull[, "beta"]
attr(post_manta_from_surf_cull, "link") <- "loglog"
saveRDS(post_manta_from_surf_cull, "output/posteriors/manta_from_surf_cull.RDS")
rm(fit_manta_from_surf_cull)

fit_surf_from_manta_cull <- brm(
  bf(clip01(y_surf) ~ y_manta + (1 | reef_label)), family = Beta(link = "logit"),
  data = calData, prior = c(prior(normal(0, 1), class = "sd"), d02_prior_coral_phi),
  chains = 4, iter = 4000, warmup = 1000, seed = coral_seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_surf_from_manta_cull <- as.matrix(as_draws_matrix(fit_surf_from_manta_cull))
colnames(post_surf_from_manta_cull)[colnames(post_surf_from_manta_cull) == "b_Intercept"] <- "alpha"
colnames(post_surf_from_manta_cull)[colnames(post_surf_from_manta_cull) == "b_y_manta"]   <- "beta"
saveRDS(post_surf_from_manta_cull, "output/posteriors/surf_from_manta_cull.RDS")
rm(fit_surf_from_manta_cull)

if (doPlot) {
  x_seq <- seq(min(calData$y_surf), max(calData$y_surf), length.out = 100)
  plot_calib_base(predict_beta_curve(post_manta_from_surf_cull, x_seq), calData$y_surf, calData$y_manta,
                   "RS Surface HC %", "Manta HC %", "Manta from RS Surface, cull-site level")
  x_seq <- seq(min(calData$y_manta), max(calData$y_manta), length.out = 100)
  plot_calib_base(predict_beta_curve(post_surf_from_manta_cull, x_seq), calData$y_manta, calData$y_surf,
                   "Manta HC %", "RS Surface HC %", "RS Surface from Manta, cull-site level")
}
