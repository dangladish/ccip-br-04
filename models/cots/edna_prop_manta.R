# edna_prop_manta.R
# eDNA proportion positive <-> manta tow, reef level (D-02 approach, at reef grain).
#   eDNAprop_from_manta_reef: n_success | trials(n_trials) ~ log(manta_density_raw + 0.01);
#     Binomial (logit), observation-level reef random effect.
#     Saves output/posteriors/eDNAprop_from_manta_reef.RDS
#   manta_from_eDNAprop_binary_reef: manta_combined | trials(n_tows) ~
#     log(edna_prop_reef + 0.01); Binomial (logit), observation-level reef random effect.
#     manta_combined = tows with scar present OR >= 1 COTS observed.
#     Saves output/posteriors/manta_from_eDNAprop_binary_reef.RDS
# Observation-level reef random effect: identifiable for a Binomial with many trials per row,
# and absorbs the extra-binomial variation. Prior: d02_priors_manta_binary.
# Sourced with local = TRUE. Requires (from the driver): reef_manta, d02_priors_manta_binary,
# cots_seed.

# ---- eDNAprop_from_manta_reef ----
fit_eDNAprop_from_manta_reef <- brm(
  bf(n_success | trials(n_trials) ~ log(manta_density_raw + 0.01) + (1 | reef_label)),
  family = binomial(link = "logit"), data = reef_manta, prior = d02_priors_manta_binary,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_eDNAprop_from_manta_reef <- as.matrix(as_draws_matrix(fit_eDNAprop_from_manta_reef))
colnames(post_eDNAprop_from_manta_reef)[colnames(post_eDNAprop_from_manta_reef) == "b_Intercept"] <- "alpha"
colnames(post_eDNAprop_from_manta_reef)[startsWith(colnames(post_eDNAprop_from_manta_reef), "b_") & colnames(post_eDNAprop_from_manta_reef) != "b_Intercept"] <- "beta"
attr(post_eDNAprop_from_manta_reef, "link") <- "logit"
attr(post_eDNAprop_from_manta_reef, "obs_re_sd") <- "sd_reef_label__Intercept"
saveRDS(post_eDNAprop_from_manta_reef, "output/posteriors/eDNAprop_from_manta_reef.RDS")
rm(fit_eDNAprop_from_manta_reef)

# ---- manta_from_eDNAprop_binary_reef ----
fit_manta_from_eDNAprop_binary_reef <- brm(
  bf(manta_combined | trials(n_tows) ~ log(edna_prop_reef + 0.01) + (1 | reef_label)),
  family = binomial(link = "logit"), data = reef_manta, prior = d02_priors_manta_binary,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_manta_from_eDNAprop_binary_reef <- as.matrix(as_draws_matrix(fit_manta_from_eDNAprop_binary_reef))
colnames(post_manta_from_eDNAprop_binary_reef)[colnames(post_manta_from_eDNAprop_binary_reef) == "b_Intercept"] <- "alpha"
colnames(post_manta_from_eDNAprop_binary_reef)[startsWith(colnames(post_manta_from_eDNAprop_binary_reef), "b_") & colnames(post_manta_from_eDNAprop_binary_reef) != "b_Intercept"] <- "beta"
attr(post_manta_from_eDNAprop_binary_reef, "link") <- "logit"
attr(post_manta_from_eDNAprop_binary_reef, "obs_re_sd") <- "sd_reef_label__Intercept"
saveRDS(post_manta_from_eDNAprop_binary_reef, "output/posteriors/manta_from_eDNAprop_binary_reef.RDS")
rm(fit_manta_from_eDNAprop_binary_reef)

if (doPlot) {
  x_seq <- seq(0, max(reef_manta$manta_density_raw), length.out = 100)
  size_rep <- round(median(reef_manta$n_trials))
  plot_calib_base(predict_binomial_curve(post_eDNAprop_from_manta_reef, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   reef_manta$manta_density_raw, reef_manta$n_success / reef_manta$n_trials, "Manta scars+observed per tow", "eDNA proportion positive (reef)", "eDNA proportion from Manta, reef level")

  x_seq <- seq(0, max(reef_manta$edna_prop_reef), length.out = 100)
  size_rep <- round(median(reef_manta$n_tows))
  plot_calib_base(predict_binomial_curve(post_manta_from_eDNAprop_binary_reef, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   reef_manta$edna_prop_reef, reef_manta$manta_combined_rate, "eDNA proportion positive", "Manta combined COTS+scar presence rate", "Manta combined COTS+scar presence rate from eDNA proportion, reef level")
}
cat("\nedna_prop_manta.R complete -- 2 posteriors saved (1 eDNA-as-response Binomial + 1 manta-as-response Binomial, reef level).\n")
