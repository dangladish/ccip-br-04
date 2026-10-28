# manta_deep_tow.R
# Coral cover calibration: manta tow <-> RS Deep, tow level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression (brms/Stan), logit link, reef random effect (1 | reef_label), D-02 prior on phi.
# Inputs: tow_level (perc_hard, pct_HC_deep).
# Saves: output/posteriors/{manta_from_deep_tow,deep_from_manta_tow}.RDS

calData <- tow_level |> filter(!is.na(pct_HC_deep), !is.na(perc_hard)) |>
  select(tow_id, reef_label, y_manta = perc_hard, y_deep = pct_HC_deep)

fit_manta_from_deep_tow <- brm(
  bf(clip01(y_manta) ~ y_deep + (1 | reef_label)), family = Beta(link = "logit"),
  data = calData, prior = c(prior(normal(0, 1), class = "sd"), d02_prior_coral_phi),
  chains = 4, iter = 2000, warmup = 1000, seed = coral_seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_manta_from_deep_tow <- as.matrix(as_draws_matrix(fit_manta_from_deep_tow))
colnames(post_manta_from_deep_tow)[colnames(post_manta_from_deep_tow) == "b_Intercept"] <- "alpha"
colnames(post_manta_from_deep_tow)[colnames(post_manta_from_deep_tow) == "b_y_deep"]    <- "beta"
saveRDS(post_manta_from_deep_tow, "output/posteriors/manta_from_deep_tow.RDS")
rm(fit_manta_from_deep_tow)

fit_deep_from_manta_tow <- brm(
  bf(clip01(y_deep) ~ y_manta + (1 | reef_label)), family = Beta(link = "logit"),
  data = calData, prior = c(prior(normal(0, 1), class = "sd"), d02_prior_coral_phi),
  chains = 4, iter = 2000, warmup = 1000, seed = coral_seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_deep_from_manta_tow <- as.matrix(as_draws_matrix(fit_deep_from_manta_tow))
colnames(post_deep_from_manta_tow)[colnames(post_deep_from_manta_tow) == "b_Intercept"] <- "alpha"
colnames(post_deep_from_manta_tow)[colnames(post_deep_from_manta_tow) == "b_y_manta"]   <- "beta"
saveRDS(post_deep_from_manta_tow, "output/posteriors/deep_from_manta_tow.RDS")
rm(fit_deep_from_manta_tow)

if (doPlot) {
  x_seq <- seq(min(calData$y_deep), max(calData$y_deep), length.out = 100)
  plot_calib_base(predict_beta_curve(post_manta_from_deep_tow, x_seq), calData$y_deep, calData$y_manta,
                   "RS Deep HC %", "Manta HC %", "Manta from RS Deep, tow level")
  x_seq <- seq(min(calData$y_manta), max(calData$y_manta), length.out = 100)
  plot_calib_base(predict_beta_curve(post_deep_from_manta_tow, x_seq), calData$y_manta, calData$y_deep,
                   "Manta HC %", "RS Deep HC %", "RS Deep from Manta, tow level")
}
