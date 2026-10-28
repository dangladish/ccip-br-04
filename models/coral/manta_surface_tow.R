# manta_surface_tow.R
# Coral cover calibration: manta tow <-> RS Surface, tow level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression (brms/Stan), reef random effect (1 | reef_label), D-02 prior on phi.
# manta_from_surf uses a loglog link; surf_from_manta uses logit.
# brms has no loglog link for Beta, so manta_from_surf is fitted with cloglog on 1 - y and
# alpha/beta are sign-flipped to the loglog(y) scale (phi is unchanged). The saved posterior
# carries attr(, "link") = "loglog", which predict_beta_curve() reads.
# Inputs: tow_level (perc_hard, pct_HC_surf).
# Saves: output/posteriors/{manta_from_surf_tow,surf_from_manta_tow}.RDS

calData <- tow_level |> filter(!is.na(pct_HC_surf), !is.na(perc_hard)) |>
  select(tow_id, reef_label, y_manta = perc_hard, y_surf = pct_HC_surf)

fit_manta_from_surf_tow <- brm(
  bf(clip01(1 - y_manta) ~ y_surf + (1 | reef_label)), family = Beta(link = "cloglog"),
  data = calData, prior = c(prior(normal(0, 1), class = "sd"), d02_prior_coral_phi),
  chains = 4, iter = 4000, warmup = 1000, seed = coral_seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_manta_from_surf_tow <- as.matrix(as_draws_matrix(fit_manta_from_surf_tow))
colnames(post_manta_from_surf_tow)[colnames(post_manta_from_surf_tow) == "b_Intercept"] <- "alpha"
colnames(post_manta_from_surf_tow)[colnames(post_manta_from_surf_tow) == "b_y_surf"]    <- "beta"
# Fitted as cloglog on 1 - y; flip signs to get the loglog(y) scale.
post_manta_from_surf_tow[, "alpha"] <- -post_manta_from_surf_tow[, "alpha"]
post_manta_from_surf_tow[, "beta"]  <- -post_manta_from_surf_tow[, "beta"]
attr(post_manta_from_surf_tow, "link") <- "loglog"
saveRDS(post_manta_from_surf_tow, "output/posteriors/manta_from_surf_tow.RDS")
rm(fit_manta_from_surf_tow)

fit_surf_from_manta_tow <- brm(
  bf(clip01(y_surf) ~ y_manta + (1 | reef_label)), family = Beta(link = "logit"),
  data = calData, prior = c(prior(normal(0, 1), class = "sd"), d02_prior_coral_phi),
  chains = 4, iter = 4000, warmup = 1000, seed = coral_seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_surf_from_manta_tow <- as.matrix(as_draws_matrix(fit_surf_from_manta_tow))
colnames(post_surf_from_manta_tow)[colnames(post_surf_from_manta_tow) == "b_Intercept"] <- "alpha"
colnames(post_surf_from_manta_tow)[colnames(post_surf_from_manta_tow) == "b_y_manta"]   <- "beta"
saveRDS(post_surf_from_manta_tow, "output/posteriors/surf_from_manta_tow.RDS")
rm(fit_surf_from_manta_tow)

if (doPlot) {
  x_seq <- seq(min(calData$y_surf), max(calData$y_surf), length.out = 100)
  plot_calib_base(predict_beta_curve(post_manta_from_surf_tow, x_seq), calData$y_surf, calData$y_manta,
                   "RS Surface HC %", "Manta HC %", "Manta from RS Surface, tow level")
  x_seq <- seq(min(calData$y_manta), max(calData$y_manta), length.out = 100)
  plot_calib_base(predict_beta_curve(post_surf_from_manta_tow, x_seq), calData$y_manta, calData$y_surf,
                   "Manta HC %", "RS Surface HC %", "RS Surface from Manta, tow level")
}
