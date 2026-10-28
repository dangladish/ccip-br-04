# edna_conc_manta.R
# eDNA concentration <-> manta tow, reef level (D-02 approach for the eDNA-as-response model).
#   eDNAconc_from_manta_reef: integer_conc_reef ~ log(manta_density_raw + 0.01); NegBin (log
#     link), offset log(n_assays_reef), no random effect (one row per reef, so a reef effect
#     cannot be separated from the NegBin dispersion).
#     Saves output/posteriors/eDNAconc_from_manta_reef.RDS
#   manta_from_eDNAconc_binary_reef: manta_combined | trials(n_tows) ~
#     log(edna_conc_rate_reef + 0.01); Binomial (logit), observation-level reef random effect.
#     manta_combined = tows with scar present OR >= 1 COTS observed.
#     Saves output/posteriors/manta_from_eDNAconc_binary_reef.RDS
# Priors: d02_priors_eDNAconc_negbin_reef, d02_priors_manta_binary.
# The NegBin curve is plotted on the per-assay scale (log_offset = 0, edna_conc_rate_reef).
# Sourced with local = TRUE. Requires (from the driver): reef_manta, the priors, cots_seed.

# ---- eDNAconc_from_manta_reef ----
fit_eDNAconc_from_manta_reef <- brm(
  bf(integer_conc_reef ~ log(manta_density_raw + 0.01) + offset(log(n_assays_reef))),
  family = negbinomial(link = "log"), data = reef_manta, prior = d02_priors_eDNAconc_negbin_reef,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_eDNAconc_from_manta_reef <- as.matrix(as_draws_matrix(fit_eDNAconc_from_manta_reef))
colnames(post_eDNAconc_from_manta_reef)[colnames(post_eDNAconc_from_manta_reef) == "b_Intercept"] <- "alpha"
colnames(post_eDNAconc_from_manta_reef)[startsWith(colnames(post_eDNAconc_from_manta_reef), "b_") & colnames(post_eDNAconc_from_manta_reef) != "b_Intercept"] <- "beta"
colnames(post_eDNAconc_from_manta_reef)[colnames(post_eDNAconc_from_manta_reef) == "shape"] <- "r"
saveRDS(post_eDNAconc_from_manta_reef, "output/posteriors/eDNAconc_from_manta_reef.RDS")
rm(fit_eDNAconc_from_manta_reef)

# ---- manta_from_eDNAconc_binary_reef ----
# Observation-level reef random effect: identifiable because each row has many trials.
fit_manta_from_eDNAconc_binary_reef <- brm(
  bf(manta_combined | trials(n_tows) ~ log(edna_conc_rate_reef + 0.01) + (1 | reef_label)),
  family = binomial(link = "logit"), data = reef_manta, prior = d02_priors_manta_binary,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_manta_from_eDNAconc_binary_reef <- as.matrix(as_draws_matrix(fit_manta_from_eDNAconc_binary_reef))
colnames(post_manta_from_eDNAconc_binary_reef)[colnames(post_manta_from_eDNAconc_binary_reef) == "b_Intercept"] <- "alpha"
colnames(post_manta_from_eDNAconc_binary_reef)[startsWith(colnames(post_manta_from_eDNAconc_binary_reef), "b_") & colnames(post_manta_from_eDNAconc_binary_reef) != "b_Intercept"] <- "beta"
attr(post_manta_from_eDNAconc_binary_reef, "link") <- "logit"
attr(post_manta_from_eDNAconc_binary_reef, "obs_re_sd") <- "sd_reef_label__Intercept"
saveRDS(post_manta_from_eDNAconc_binary_reef, "output/posteriors/manta_from_eDNAconc_binary_reef.RDS")
rm(fit_manta_from_eDNAconc_binary_reef)

if (doPlot) {
  x_seq <- seq(0, max(reef_manta$manta_density_raw), length.out = 100)
  plot_calib_base(predict_negbin_curve(post_eDNAconc_from_manta_reef, x_seq, log_offset = 0, eps = 0.01),
                   reef_manta$manta_density_raw, reef_manta$edna_conc_rate_reef, "Manta scars+observed per tow", "eDNA concentration (mean per assay)", "eDNA concentration from Manta, reef level")

  x_seq <- seq(0, max(reef_manta$edna_conc_rate_reef), length.out = 100)
  size_rep <- round(median(reef_manta$n_tows))
  plot_calib_base(predict_binomial_curve(post_manta_from_eDNAconc_binary_reef, x_seq, eps = 0.01, link = "logit", size = size_rep),
                   reef_manta$edna_conc_rate_reef, reef_manta$manta_combined_rate, "eDNA conc. rate (conc/assay)", "Manta combined COTS+scar presence rate", "Manta combined COTS+scar presence rate from eDNA concentration, reef level")
}
cat("\nedna_conc_manta.R complete -- 2 posteriors saved (1 eDNA-as-response NegBin + 1 manta-as-response Binomial, reef level).\n")
