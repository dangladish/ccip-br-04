# surface_deep_tow.R
# Coral cover calibration: RS Surface <-> RS Deep, tow level (both directions).
# Response ~ predictor: hard coral cover (proportion) of one method ~ the other's.
# Beta regression (brms/Stan), logit link, reef random effect (1 | reef_label), D-02 prior on phi.
# Inputs: tow_level (pct_HC_surf, pct_HC_deep).
# Saves: output/posteriors/{surf_from_deep_tow,deep_from_surf_tow}.RDS

calData <- tow_level |> filter(!is.na(pct_HC_surf), !is.na(pct_HC_deep)) |>
  select(tow_id, reef_label, y_surf = pct_HC_surf, y_deep = pct_HC_deep)

fit_surf_from_deep_tow <- brm(
  bf(clip01(y_surf) ~ y_deep + (1 | reef_label)), family = Beta(link = "logit"),
  data = calData, prior = c(prior(normal(0, 1), class = "sd"), d02_prior_coral_phi),
  chains = 4, iter = 2000, warmup = 1000, seed = coral_seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_surf_from_deep_tow <- as.matrix(as_draws_matrix(fit_surf_from_deep_tow))
colnames(post_surf_from_deep_tow)[colnames(post_surf_from_deep_tow) == "b_Intercept"] <- "alpha"
colnames(post_surf_from_deep_tow)[colnames(post_surf_from_deep_tow) == "b_y_deep"]    <- "beta"
saveRDS(post_surf_from_deep_tow, "output/posteriors/surf_from_deep_tow.RDS")
rm(fit_surf_from_deep_tow)

fit_deep_from_surf_tow <- brm(
  bf(clip01(y_deep) ~ y_surf + (1 | reef_label)), family = Beta(link = "logit"),
  data = calData, prior = c(prior(normal(0, 1), class = "sd"), d02_prior_coral_phi),
  chains = 4, iter = 2000, warmup = 1000, seed = coral_seed,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_deep_from_surf_tow <- as.matrix(as_draws_matrix(fit_deep_from_surf_tow))
colnames(post_deep_from_surf_tow)[colnames(post_deep_from_surf_tow) == "b_Intercept"] <- "alpha"
colnames(post_deep_from_surf_tow)[colnames(post_deep_from_surf_tow) == "b_y_surf"]    <- "beta"
saveRDS(post_deep_from_surf_tow, "output/posteriors/deep_from_surf_tow.RDS")
rm(fit_deep_from_surf_tow)

if (doPlot) {
  x_seq <- seq(min(calData$y_deep), max(calData$y_deep), length.out = 100)
  plot_calib_base(predict_beta_curve(post_surf_from_deep_tow, x_seq), calData$y_deep, calData$y_surf,
                   "RS Deep HC %", "RS Surface HC %", "RS Surface from RS Deep, tow level")
  x_seq <- seq(min(calData$y_surf), max(calData$y_surf), length.out = 100)
  plot_calib_base(predict_beta_curve(post_deep_from_surf_tow, x_seq), calData$y_surf, calData$y_deep,
                   "RS Surface HC %", "RS Deep HC %", "RS Deep from RS Surface, tow level")
}
