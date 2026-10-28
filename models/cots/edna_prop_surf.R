# edna_prop_surf.R
# eDNA proportion positive <-> RS Surface, reef level.
#   eDNAprop_from_surf_reef: n_success | trials(n_trials) ~ log(surf_density + 0.01); Binomial
#     (logit), observation-level reef random effect, d02_priors_manta_binary.
#     surf_density = RS Surface COTS per tow.
#     Saves output/posteriors/eDNAprop_from_surf_reef.RDS
#   surf_from_eDNAprop_reef: cots_reefscan_surface ~ log(edna_prop_reef + 0.01); NegBin (log
#     link), offset log(effort_surf) with effort_surf = n_tows, no random effect (one row per
#     reef), d02_priors_eDNAconc_negbin_reef (fairly informative shape prior, since few reefs
#     cannot identify the dispersion).
#     Saves output/posteriors/surf_from_eDNAprop_reef.RDS
# Sourced with local = TRUE. Requires (from the driver): reef_surf, d02_priors_manta_binary,
# d02_priors_eDNAconc_negbin_reef, cots_seed.

# ---- eDNAprop_from_surf_reef ----
# Observation-level reef random effect: identifiable because each row has many trials.
fit_eDNAprop_from_surf_reef <- brm(
  bf(n_success | trials(n_trials) ~ log(surf_density + 0.01) + (1 | reef_label)),
  family = binomial(link = "logit"), data = reef_surf, prior = d02_priors_manta_binary,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_eDNAprop_from_surf_reef <- as.matrix(as_draws_matrix(fit_eDNAprop_from_surf_reef))
colnames(post_eDNAprop_from_surf_reef)[colnames(post_eDNAprop_from_surf_reef) == "b_Intercept"] <- "alpha"
colnames(post_eDNAprop_from_surf_reef)[startsWith(colnames(post_eDNAprop_from_surf_reef), "b_") & colnames(post_eDNAprop_from_surf_reef) != "b_Intercept"] <- "beta"
attr(post_eDNAprop_from_surf_reef, "link") <- "logit"
attr(post_eDNAprop_from_surf_reef, "obs_re_sd") <- "sd_reef_label__Intercept"
saveRDS(post_eDNAprop_from_surf_reef, "output/posteriors/eDNAprop_from_surf_reef.RDS")
rm(fit_eDNAprop_from_surf_reef)

# ---- surf_from_eDNAprop_reef ----
sepr <- reef_surf |> mutate(effort_surf = n_tows)

fit_surf_from_eDNAprop_reef <- brm(
  bf(cots_reefscan_surface ~ log(edna_prop_reef + 0.01) + offset(log(effort_surf))),
  family = negbinomial(link = "log"), data = sepr, prior = d02_priors_eDNAconc_negbin_reef,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_surf_from_eDNAprop_reef <- as.matrix(as_draws_matrix(fit_surf_from_eDNAprop_reef))
colnames(post_surf_from_eDNAprop_reef)[colnames(post_surf_from_eDNAprop_reef) == "b_Intercept"] <- "alpha"
colnames(post_surf_from_eDNAprop_reef)[startsWith(colnames(post_surf_from_eDNAprop_reef), "b_") & colnames(post_surf_from_eDNAprop_reef) != "b_Intercept"] <- "beta"
colnames(post_surf_from_eDNAprop_reef)[colnames(post_surf_from_eDNAprop_reef) == "shape"] <- "r"
saveRDS(post_surf_from_eDNAprop_reef, "output/posteriors/surf_from_eDNAprop_reef.RDS")
rm(fit_surf_from_eDNAprop_reef)

if (doPlot) {
  x_seq <- seq(0, max(reef_surf$surf_density), length.out = 100)
  size_rep <- round(median(reef_surf$n_trials))
  plot_calib_base(predict_binomial_curve(post_eDNAprop_from_surf_reef, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   reef_surf$surf_density, reef_surf$n_success / reef_surf$n_trials, "RS Surface density (COTS/tow)", "eDNA proportion positive (reef)", "eDNA proportion from RS Surface, reef level")

  x_seq <- seq(0, max(sepr$edna_prop_reef), length.out = 100)
  effort_rep <- mean(sepr$effort_surf)
  plot_calib_base(predict_negbin_curve_density(post_surf_from_eDNAprop_reef, x_seq, effort_rep, eps = 0.01),
                   sepr$edna_prop_reef, sepr$cots_reefscan_surface / sepr$effort_surf, "eDNA proportion positive", "RS Surface COTS density (per tow)", "RS Surface from eDNA proportion, reef level")
}
cat("\nedna_prop_surf.R complete -- 2 posteriors saved, reef level.\n")
