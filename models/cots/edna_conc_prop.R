# edna_conc_prop.R
# eDNA concentration <-> eDNA proportion positive, reef level (both from the same eDNA samples;
# follows the D-02 self-comparison of its two eDNA response scales).
#   eDNAprop_from_eDNAconc_reef: n_success | trials(n_trials) ~ log(edna_conc_rate_reef + 0.01);
#     Binomial (logit), observation-level reef random effect.
#     Saves output/posteriors/eDNAprop_from_eDNAconc_reef.RDS
#   eDNAconc_from_eDNAprop_reef: integer_conc_reef ~ log(edna_prop_reef + 0.01); NegBin (log
#     link), offset log(n_assays_reef), no random effect (one row per reef).
#     Saves output/posteriors/eDNAconc_from_eDNAprop_reef.RDS
# Sourced with local = TRUE. Requires (from the driver): edna_reef, d02_priors_manta_binary,
# d02_priors_eDNAconc_negbin_reef, cots_seed.

# ---- eDNAprop_from_eDNAconc_reef ----
# Observation-level reef random effect: identifiable because each row has many trials.
fit_eDNAprop_from_eDNAconc_reef <- brm(
  bf(n_success | trials(n_trials) ~ log(edna_conc_rate_reef + 0.01) + (1 | reef_label)),
  family = binomial(link = "logit"), data = edna_reef, prior = d02_priors_manta_binary,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_eDNAprop_from_eDNAconc_reef <- as.matrix(as_draws_matrix(fit_eDNAprop_from_eDNAconc_reef))
colnames(post_eDNAprop_from_eDNAconc_reef)[colnames(post_eDNAprop_from_eDNAconc_reef) == "b_Intercept"] <- "alpha"
colnames(post_eDNAprop_from_eDNAconc_reef)[startsWith(colnames(post_eDNAprop_from_eDNAconc_reef), "b_") & colnames(post_eDNAprop_from_eDNAconc_reef) != "b_Intercept"] <- "beta"
attr(post_eDNAprop_from_eDNAconc_reef, "link") <- "logit"
attr(post_eDNAprop_from_eDNAconc_reef, "obs_re_sd") <- "sd_reef_label__Intercept"
saveRDS(post_eDNAprop_from_eDNAconc_reef, "output/posteriors/eDNAprop_from_eDNAconc_reef.RDS")
rm(fit_eDNAprop_from_eDNAconc_reef)

# ---- eDNAconc_from_eDNAprop_reef ----
fit_eDNAconc_from_eDNAprop_reef <- brm(
  bf(integer_conc_reef ~ log(edna_prop_reef + 0.01) + offset(log(n_assays_reef))),
  family = negbinomial(link = "log"), data = edna_reef, prior = d02_priors_eDNAconc_negbin_reef,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_eDNAconc_from_eDNAprop_reef <- as.matrix(as_draws_matrix(fit_eDNAconc_from_eDNAprop_reef))
colnames(post_eDNAconc_from_eDNAprop_reef)[colnames(post_eDNAconc_from_eDNAprop_reef) == "b_Intercept"] <- "alpha"
colnames(post_eDNAconc_from_eDNAprop_reef)[startsWith(colnames(post_eDNAconc_from_eDNAprop_reef), "b_") & colnames(post_eDNAconc_from_eDNAprop_reef) != "b_Intercept"] <- "beta"
colnames(post_eDNAconc_from_eDNAprop_reef)[colnames(post_eDNAconc_from_eDNAprop_reef) == "shape"] <- "r"
saveRDS(post_eDNAconc_from_eDNAprop_reef, "output/posteriors/eDNAconc_from_eDNAprop_reef.RDS")
rm(fit_eDNAconc_from_eDNAprop_reef)

if (doPlot) {
  x_seq <- seq(0, max(edna_reef$edna_conc_rate_reef), length.out = 100)
  size_rep <- round(median(edna_reef$n_trials))
  plot_calib_base(predict_binomial_curve(post_eDNAprop_from_eDNAconc_reef, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   edna_reef$edna_conc_rate_reef, edna_reef$edna_prop_reef, "eDNA conc. rate (conc/assay)", "eDNA proportion positive", "eDNA proportion from eDNA concentration, reef level")

  x_seq <- seq(0, max(edna_reef$edna_prop_reef), length.out = 100)
  plot_calib_base(predict_negbin_curve(post_eDNAconc_from_eDNAprop_reef, x_seq, log_offset = 0, eps = 0.01),
                   edna_reef$edna_prop_reef, edna_reef$edna_conc_rate_reef, "eDNA proportion positive", "eDNA concentration (mean per assay)", "eDNA concentration from eDNA proportion, reef level")
}
cat("\nedna_conc_prop.R complete -- 2 posteriors saved, reef level.\n")
