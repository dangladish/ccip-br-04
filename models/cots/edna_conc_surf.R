# edna_conc_surf.R
# eDNA concentration <-> RS Surface, reef level.
#   eDNAconc_from_surf_reef: integer_conc_reef ~ log(surf_density + 0.01); NegBin (log link),
#     offset log(n_assays_reef), no random effect. surf_density = RS Surface COTS per tow.
#     Saves output/posteriors/eDNAconc_from_surf_reef.RDS
#   surf_from_eDNAconc_reef: cots_reefscan_surface ~ log(edna_conc_rate_reef + 0.01); NegBin
#     (log link), offset log(effort_surf) with effort_surf = n_tows, no random effect.
#     Saves output/posteriors/surf_from_eDNAconc_reef.RDS
# No reef random effect: one row per reef, so it cannot be separated from the NegBin dispersion.
# Prior d02_priors_eDNAconc_negbin_reef has a fairly informative shape prior because the small
# number of reefs cannot identify the dispersion.
# Sourced with local = TRUE. Requires (from the driver): reef_surf,
# d02_priors_eDNAconc_negbin_reef.

# ---- eDNAconc_from_surf_reef ----
fit_eDNAconc_from_surf_reef <- brm(
  bf(integer_conc_reef ~ log(surf_density + 0.01) + offset(log(n_assays_reef))),
  family = negbinomial(link = "log"), data = reef_surf, prior = d02_priors_eDNAconc_negbin_reef,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_eDNAconc_from_surf_reef <- as.matrix(as_draws_matrix(fit_eDNAconc_from_surf_reef))
colnames(post_eDNAconc_from_surf_reef)[colnames(post_eDNAconc_from_surf_reef) == "b_Intercept"] <- "alpha"
colnames(post_eDNAconc_from_surf_reef)[startsWith(colnames(post_eDNAconc_from_surf_reef), "b_") & colnames(post_eDNAconc_from_surf_reef) != "b_Intercept"] <- "beta"
colnames(post_eDNAconc_from_surf_reef)[colnames(post_eDNAconc_from_surf_reef) == "shape"] <- "r"
saveRDS(post_eDNAconc_from_surf_reef, "output/posteriors/eDNAconc_from_surf_reef.RDS")
rm(fit_eDNAconc_from_surf_reef)

# ---- surf_from_eDNAconc_reef ----
ser <- reef_surf |> mutate(effort_surf = n_tows)

fit_surf_from_eDNAconc_reef <- brm(
  bf(cots_reefscan_surface ~ log(edna_conc_rate_reef + 0.01) + offset(log(effort_surf))),
  family = negbinomial(link = "log"), data = ser, prior = d02_priors_eDNAconc_negbin_reef,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
post_surf_from_eDNAconc_reef <- as.matrix(as_draws_matrix(fit_surf_from_eDNAconc_reef))
colnames(post_surf_from_eDNAconc_reef)[colnames(post_surf_from_eDNAconc_reef) == "b_Intercept"] <- "alpha"
colnames(post_surf_from_eDNAconc_reef)[startsWith(colnames(post_surf_from_eDNAconc_reef), "b_") & colnames(post_surf_from_eDNAconc_reef) != "b_Intercept"] <- "beta"
colnames(post_surf_from_eDNAconc_reef)[colnames(post_surf_from_eDNAconc_reef) == "shape"] <- "r"
saveRDS(post_surf_from_eDNAconc_reef, "output/posteriors/surf_from_eDNAconc_reef.RDS")
rm(fit_surf_from_eDNAconc_reef)

if (doPlot) {
  x_seq <- seq(0, max(reef_surf$surf_density), length.out = 100)
  plot_calib_base(predict_negbin_curve(post_eDNAconc_from_surf_reef, x_seq, log_offset = 0, eps = 0.01),
                   reef_surf$surf_density, reef_surf$edna_conc_rate_reef, "RS Surface density (COTS/tow)", "eDNA concentration (mean per assay)", "eDNA concentration from RS Surface, reef level")

  x_seq <- seq(0, max(ser$edna_conc_rate_reef), length.out = 100)
  effort_rep <- mean(ser$effort_surf)
  plot_calib_base(predict_negbin_curve_density(post_surf_from_eDNAconc_reef, x_seq, effort_rep, eps = 0.01),
                   ser$edna_conc_rate_reef, ser$cots_reefscan_surface / ser$effort_surf, "eDNA conc. rate (conc/assay)", "RS Surface COTS density (per tow)", "RS Surface from eDNA concentration, reef level")
}
cat("\nedna_conc_surf.R complete -- 2 posteriors saved, reef level.\n")
