# surf_deep_cull_cots.R
# RS Surface <-> RS Deep COTS count calibration at cull-site level.
# Posteriors saved: surf_from_deep_cull_cots, deep_from_surf_cull_cots
# (NegBin, offset log(effort of the response method), reef random effect).
# deep_from_surf_cull_cots uses d02_priors_eDNAconc_negbin (R/cots_shared.R) for the
# NegBin shape: sparse, mostly-zero RS Deep counts do not constrain the default prior.
# Requires from the driver: sd_cull (cots_crossmethod_data_prep.R).
## d02_priors_eDNAconc_negbin (R/cots_shared.R).

post_surf_from_deep_cull_cots <- fit_crossmethod_pair(sd_cull$y_surf, sd_cull$effort_surf, sd_cull$y_deep, sd_cull$effort_deep, sd_cull$reef_label, 0.01, "surf_from_deep_cull_cots", "cull")
post_deep_from_surf_cull_cots <- fit_crossmethod_pair(sd_cull$y_deep, sd_cull$effort_deep, sd_cull$y_surf, sd_cull$effort_surf, sd_cull$reef_label, 0.01, "deep_from_surf_cull_cots", "cull", prior_override = d02_priors_eDNAconc_negbin)

if (doPlot) {
  plot_crossmethod_pair_density(post_surf_from_deep_cull_cots, sd_cull$y_surf, sd_cull$effort_surf, sd_cull$y_deep, sd_cull$effort_deep, "RS Deep density (COTS/tow)", "RS Surface COTS density (per tow)", "RS Surface density from RS Deep, cull-site level")
  plot_crossmethod_pair_density(post_deep_from_surf_cull_cots, sd_cull$y_deep, sd_cull$effort_deep, sd_cull$y_surf, sd_cull$effort_surf, "RS Surface density (COTS/tow)", "RS Deep COTS density (per tow)", "RS Deep density from RS Surface, cull-site level")
}
