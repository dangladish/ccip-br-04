## plot_report_figures.R
##
## Combined calibration figures (one multi-panel PNG each, no refitting):
##   1. Coral cover, manta <-> RS Surface, 3 levels x 2 directions
##   2. COTS cross-method, cull-site level, 3 pairs x 2 directions
##      (manta-as-response panels use the NegBin "combined" models)
##   3. eDNA <-> manta / RS Surface, reef level, 3 pairs x 2 directions
## Inputs: output/posteriors/*.RDS, data/derived/*.RDS. Outputs: PNGs in output/figures/.
## Panel data are built as in each model file's own plotting block.
## Both shared files define plot_calib_base() with different scaling, so neither
## is called here; the local plot_row() draws all panels (coral predictors
## and curves come from predict_*_curve()).

library(tidyverse)
library(brms)

fig_dir  <- "output/figures"
post_dir <- "output/posteriors"

load_post <- function(label) readRDS(file.path(post_dir, paste0(label, ".RDS")))

## Draw one calibration panel; scale100 converts proportions to % (coral only).
plot_row <- function(pred, obs_x, obs_y, xlab, ylab, main, scale100 = FALSE, ylim = NULL) {
  mult <- if (scale100) 100 else 1
  cols <- cbind(pred$mean, pred$lo95, pred$hi95, pred$pi_lo95, pred$pi_hi95) * mult
  x <- pred$x * mult; ox <- obs_x * mult; oy <- obs_y * mult
  if (is.null(ylim)) ylim <- range(c(cols, oy), na.rm = TRUE); xlim <- range(c(x, ox), na.rm = TRUE)
  # Coral cover panels start the x-axis at 0
  if (scale100) xlim[1] <- min(0, xlim[1])
  matplot(x, cols, type = "l", lty = c(1, 2, 2, 3, 3), lwd = c(2, 1, 1, 1, 1),
          col = c("black", "blue", "blue", "grey45", "grey45"),
          xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim)
  points(ox, oy, pch = 1)
  legend("topleft", bty = "n", cex = 0.7,
         legend = c("Mean", "95% CrI", "95% PI", "Observations"),
         lty = c(1, 2, 3, NA), pch = c(NA, NA, NA, 1),
         col = c("black", "blue", "grey45", "black"), lwd = 2)
}

# ---- coral cover, manta <-> RS Surface ----

tow_level       <- readRDS("data/derived/tow_level.RDS")
cull_site_level <- readRDS("data/derived/cull_site_level.RDS")
reef_level      <- readRDS("data/derived/reef_level.RDS")

source("R/coral_shared.R")  # predict_beta_curve()

## Calibration data, as in the coral model files
calData_tow <- tow_level |> filter(!is.na(pct_HC_surf), !is.na(perc_hard)) |>
  select(tow_id, reef_label, y_manta = perc_hard, y_surf = pct_HC_surf)

calData_cull <- cull_site_level |>
  select(cull_site, reef_label, y_manta = perc_hard, y_surf = pct_HC_surf) |>
  filter(!is.na(y_manta), !is.na(y_surf))

calData_reef <- reef_level |> filter(has_surf, !is.na(mean_HC_manta), !is.na(mean_HC_surf)) |>
  transmute(y_manta = mean_HC_manta, y_surf = mean_HC_surf)

levels_list <- list(
  tow  = list(d = calData_tow,  word = "tow",        lab_m = "manta_from_surf_tow",  lab_s = "surf_from_manta_tow"),
  cull = list(d = calData_cull, word = "cull-site",   lab_m = "manta_from_surf_cull", lab_s = "surf_from_manta_cull"),
  reef = list(d = calData_reef, word = "reef",        lab_m = "manta_from_surf_reef", lab_s = "surf_from_manta_reef")
)

panel_letters <- LETTERS[1:6]
panel_i <- 0

png(file.path(fig_dir, "report_coral_manta_surface_all_levels.png"), width = 1800, height = 1950, res = 120)
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2.5, 1), cex = 1.3, cex.main = 1.15)
for (lvl in levels_list) {
  d <- lvl$d
  n <- nrow(d)

  panel_i <- panel_i + 1
  post_m <- load_post(lvl$lab_m)
  x_seq  <- seq(min(d$y_surf), max(d$y_surf), length.out = 100)
  plot_row(predict_beta_curve(post_m, x_seq), d$y_surf, d$y_manta,
           "RS Surface HC %", "Manta HC %",
           sprintf("%s) Manta from RS Surface, %s level (n=%d)", panel_letters[panel_i], lvl$word, n),
           scale100 = TRUE)

  panel_i <- panel_i + 1
  post_s <- load_post(lvl$lab_s)
  x_seq  <- seq(min(d$y_manta), max(d$y_manta), length.out = 100)
  plot_row(predict_beta_curve(post_s, x_seq), d$y_manta, d$y_surf,
           "Manta HC %", "RS Surface HC %",
           sprintf("%s) RS Surface from Manta, %s level (n=%d)", panel_letters[panel_i], lvl$word, n),
           scale100 = TRUE)
}
dev.off()
cat("Saved", file.path(fig_dir, "report_coral_manta_surface_all_levels.png"), "\n")

# ---- COTS cross-method, cull-site level ----

source("R/cots_shared.R")                     # predict_negbin_curve_density()
source("models/cots/cots_crossmethod_data_prep.R")  # builds ms_cull, md_cull, sd_cull

cots_panel_i <- 0

png(file.path(fig_dir, "report_cots_crossmethod_cull_all_pairs.png"), width = 1800, height = 1950, res = 120)
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2.5, 1), cex = 1.3, cex.main = 1.15)

## Row 1: Manta <-> RS Surface
cots_panel_i <- cots_panel_i + 1
post <- load_post("manta_from_surf_combined_cull_cots")
x_seq <- seq(0, max(ms_cull$y_surf / ms_cull$effort_surf), length.out = 100)
plot_row(predict_negbin_curve_density(post, x_seq, mean(ms_cull$effort_manta), eps = 0.01),
         ms_cull$y_surf / ms_cull$effort_surf, ms_cull$manta_combined / ms_cull$effort_manta,
         "RS Surface density (COTS/tow)", "Manta combined COTS+scar density (per tow)",
         sprintf("%s) Manta combined density from RS Surface (n=%d)", LETTERS[cots_panel_i], nrow(ms_cull)))

cots_panel_i <- cots_panel_i + 1
post <- load_post("surf_from_manta_cull_cots")
dens_manta_comb_s <- (ms_cull$y_manta + ms_cull$scar_count) / ms_cull$effort_manta
x_seq <- seq(0, max(dens_manta_comb_s), length.out = 100)
plot_row(predict_negbin_curve_density(post, x_seq, mean(ms_cull$effort_surf), eps = 0.01),
         dens_manta_comb_s, ms_cull$y_surf / ms_cull$effort_surf,
         "Manta density (combined COTS+scar count / tow)", "RS Surface COTS density (per tow)",
         sprintf("%s) RS Surface density from Manta combined (n=%d)", LETTERS[cots_panel_i], nrow(ms_cull)))

## Row 2: Manta <-> RS Deep
cots_panel_i <- cots_panel_i + 1
post <- load_post("manta_from_deep_combined_cull_cots")
x_seq <- seq(0, max(md_cull$y_deep / md_cull$effort_deep), length.out = 100)
plot_row(predict_negbin_curve_density(post, x_seq, mean(md_cull$effort_manta), eps = 0.01),
         md_cull$y_deep / md_cull$effort_deep, md_cull$manta_combined / md_cull$effort_manta,
         "RS Deep density (COTS/tow)", "Manta combined COTS+scar density (per tow)",
         sprintf("%s) Manta combined density from RS Deep (n=%d)", LETTERS[cots_panel_i], nrow(md_cull)))

cots_panel_i <- cots_panel_i + 1
post <- load_post("deep_from_manta_cull_cots")
dens_manta_comb_d <- (md_cull$y_manta + md_cull$scar_count) / md_cull$effort_manta
x_seq <- seq(0, max(dens_manta_comb_d), length.out = 100)
plot_row(predict_negbin_curve_density(post, x_seq, mean(md_cull$effort_deep), eps = 0.01),
         dens_manta_comb_d, md_cull$y_deep / md_cull$effort_deep,
         "Manta density (combined COTS+scar count / tow)", "RS Deep COTS density (per tow)",
         sprintf("%s) RS Deep density from Manta combined (n=%d)", LETTERS[cots_panel_i], nrow(md_cull)))

## Row 3: RS Surface <-> RS Deep
cots_panel_i <- cots_panel_i + 1
post <- load_post("surf_from_deep_cull_cots")
x_seq <- seq(0, max(sd_cull$y_deep / sd_cull$effort_deep), length.out = 100)
plot_row(predict_negbin_curve_density(post, x_seq, mean(sd_cull$effort_surf), eps = 0.01),
         sd_cull$y_deep / sd_cull$effort_deep, sd_cull$y_surf / sd_cull$effort_surf,
         "RS Deep density (COTS/tow)", "RS Surface COTS density (per tow)",
         sprintf("%s) RS Surface density from RS Deep (n=%d)", LETTERS[cots_panel_i], nrow(sd_cull)))

cots_panel_i <- cots_panel_i + 1
post <- load_post("deep_from_surf_cull_cots")
x_seq <- seq(0, max(sd_cull$y_surf / sd_cull$effort_surf), length.out = 100)
plot_row(predict_negbin_curve_density(post, x_seq, mean(sd_cull$effort_deep), eps = 0.01),
         sd_cull$y_surf / sd_cull$effort_surf, sd_cull$y_deep / sd_cull$effort_deep,
         "RS Surface density (COTS/tow)", "RS Deep COTS density (per tow)",
         sprintf("%s) RS Deep density from RS Surface (n=%d)", LETTERS[cots_panel_i], nrow(sd_cull)))

dev.off()
cat("Saved", file.path(fig_dir, "report_cots_crossmethod_cull_all_pairs.png"), "\n")

# ---- eDNA, reef level ----
## eDNAconc_from_manta_reef is predicted with log_offset = 0 (no assay offset in that model),
## unlike the density-scale curves above.
edna_raw <- readRDS("data/derived/edna_raw.RDS")
reef_xwalk <- readRDS("data/derived/reef_xwalk.RDS")
xwalk    <- reef_xwalk |> filter(has_edna) |> select(edna_reef, reef_label)

source("models/cots/edna_calibration_data_prep.R")  # builds edna_reef, reef_manta, reef_surf

edna_panel_i <- 0

png(file.path(fig_dir, "report_edna_reef_pairs.png"), width = 1800, height = 1950, res = 120)
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2.5, 1), cex = 1.3, cex.main = 1.15)

## Row 1: eDNA proportion <-> manta
edna_panel_i <- edna_panel_i + 1
post <- load_post("eDNAprop_from_manta_reef")
x_seq <- seq(0, max(reef_manta$manta_density_raw), length.out = 100)
size_rep <- round(median(reef_manta$n_trials))
plot_row(predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep),
         reef_manta$manta_density_raw, reef_manta$n_success / reef_manta$n_trials,
         "Manta scars+observed per tow", "eDNA proportion positive",
         sprintf("%s) eDNA proportion from Manta (n=%d)", LETTERS[edna_panel_i], nrow(reef_manta)))

edna_panel_i <- edna_panel_i + 1
post <- load_post("manta_from_eDNAprop_binary_reef")
x_seq <- seq(0, max(reef_manta$edna_prop_reef), length.out = 100)
size_rep <- round(median(reef_manta$n_tows))
plot_row(predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep),
         reef_manta$edna_prop_reef, reef_manta$manta_combined_rate,
         "eDNA proportion positive", "Manta combined presence rate",
         sprintf("%s) Manta combined from eDNA proportion (n=%d)", LETTERS[edna_panel_i], nrow(reef_manta)))

## Row 2: eDNA proportion <-> RS Surface
edna_panel_i <- edna_panel_i + 1
post <- load_post("eDNAprop_from_surf_reef")
x_seq <- seq(0, max(reef_surf$surf_density), length.out = 100)
size_rep <- round(median(reef_surf$n_trials))
plot_row(predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep),
         reef_surf$surf_density, reef_surf$n_success / reef_surf$n_trials,
         "RS Surface density (COTS/tow)", "eDNA proportion positive",
         sprintf("%s) eDNA proportion from RS Surface (n=%d)", LETTERS[edna_panel_i], nrow(reef_surf)))

edna_panel_i <- edna_panel_i + 1
sepr <- reef_surf |> mutate(effort_surf = n_tows)
post <- load_post("surf_from_eDNAprop_reef")
x_seq <- seq(0, max(sepr$edna_prop_reef), length.out = 100)
effort_rep <- mean(sepr$effort_surf)
plot_row(predict_negbin_curve_density(post, x_seq, effort_rep, eps = 0.01),
         sepr$edna_prop_reef, sepr$cots_reefscan_surface / sepr$effort_surf,
         "eDNA proportion positive", "RS Surface density (COTS/tow)",
         sprintf("%s) RS Surface from eDNA proportion (n=%d)", LETTERS[edna_panel_i], nrow(sepr)),
         ylim = c(0, 1))

## Row 3: eDNA concentration <-> manta
edna_panel_i <- edna_panel_i + 1
post <- load_post("eDNAconc_from_manta_reef")
x_seq <- seq(0, max(reef_manta$manta_density_raw), length.out = 100)
plot_row(predict_negbin_curve(post, x_seq, log_offset = 0, eps = 0.01),
         reef_manta$manta_density_raw, reef_manta$edna_conc_rate_reef,
         "Manta scars+observed per tow", "eDNA concentration (per assay)",
         sprintf("%s) eDNA concentration from Manta (n=%d)", LETTERS[edna_panel_i], nrow(reef_manta)))

edna_panel_i <- edna_panel_i + 1
## Manta combined count ~ eDNA concentration (NegBin, no group term: one row per reef).
## This model is fitted here, for this figure only.
fit_manta_from_eDNAconc_combined_reef <- brm(
  bf(manta_combined ~ log(edna_conc_rate_reef + 0.01) + offset(log(n_tows))),
  family = negbinomial(link = "log"), data = reef_manta, prior = d02_priors_eDNAconc_negbin_reef,
  chains = 4, iter = 4000, warmup = 1000, seed = cots_seed, control = list(adapt_delta = 0.95, max_treedepth = 12)
)
cat("manta_from_eDNAconc_combined_reef: max Rhat =", max(brms::rhat(fit_manta_from_eDNAconc_combined_reef)), "\n")
post_manta_from_eDNAconc_combined_reef <- as.matrix(as_draws_matrix(fit_manta_from_eDNAconc_combined_reef))
colnames(post_manta_from_eDNAconc_combined_reef)[colnames(post_manta_from_eDNAconc_combined_reef) == "b_Intercept"] <- "alpha"
colnames(post_manta_from_eDNAconc_combined_reef)[startsWith(colnames(post_manta_from_eDNAconc_combined_reef), "b_") & colnames(post_manta_from_eDNAconc_combined_reef) != "b_Intercept"] <- "beta"
colnames(post_manta_from_eDNAconc_combined_reef)[colnames(post_manta_from_eDNAconc_combined_reef) == "shape"] <- "r"
saveRDS(post_manta_from_eDNAconc_combined_reef, "output/posteriors/manta_from_eDNAconc_combined_reef.RDS")
rm(fit_manta_from_eDNAconc_combined_reef)

x_seq <- seq(0, max(reef_manta$edna_conc_rate_reef), length.out = 100)
effort_rep <- mean(reef_manta$n_tows)
plot_row(predict_negbin_curve_density(post_manta_from_eDNAconc_combined_reef, x_seq, effort_rep, eps = 0.01),
         reef_manta$edna_conc_rate_reef, reef_manta$manta_combined / reef_manta$n_tows,
         "eDNA concentration (per assay)", "Manta combined COTS+scar density (per tow)",
sprintf("%s) Manta density from eDNA concentration (n=%d)", LETTERS[edna_panel_i], nrow(reef_manta)))

dev.off()
cat("Saved", file.path(fig_dir, "report_edna_reef_pairs.png"), "\n")
