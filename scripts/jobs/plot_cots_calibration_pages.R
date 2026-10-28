## plot_cots_calibration_pages.R
##
## COTS calibration figure pages (one page per response tool, one panel per
## predictor tool, stacked by spatial level). Pages: manta rate, manta presence,
## RS Surface, RS Deep, eDNA concentration, eDNA proportion (the last two are
## reef level only). Each predictor tool keeps a fixed column on every page;
## rows are tow / cull-site / reef level. Mean = black, lower CI/PI = green,
## upper CI/PI = blue (CI thin, PI thick).
## Inputs: output/posteriors/*.RDS and data/derived/*.RDS. Outputs: one PNG per
## page in output/figures/. No refitting.

library(tidyverse)

fig_dir  <- "output/figures"
post_dir <- "output/posteriors"
load_post <- function(label) readRDS(file.path(post_dir, paste0(label, ".RDS")))

tow_level       <- readRDS("data/derived/tow_level.RDS")
cull_site_level <- readRDS("data/derived/cull_site_level.RDS")
reef_level      <- readRDS("data/derived/reef_level.RDS")
reef_xwalk      <- readRDS("data/derived/reef_xwalk.RDS")
edna_raw        <- readRDS("data/derived/edna_raw.RDS")
xwalk           <- reef_xwalk |> filter(has_edna) |> select(edna_reef, reef_label)

source("R/cots_shared.R")                     # prediction-curve functions
source("models/cots/cots_crossmethod_data_prep.R")  # ms_*/md_*/sd_* tables
source("models/cots/edna_calibration_data_prep.R")  # reef_manta, reef_surf

# ---- Drawing primitives ----

## One predictor panel: mean, CI (thin) and PI (thick); `label` is printed inside
## the box. ylab is blank because each page has one shared y label.
d02_panel <- function(pred, obs_x, obs_y, xlab, label, ylab = "") {
  ylim <- range(c(pred$mean, pred$lo95, pred$hi95, pred$pi_lo95, pred$pi_hi95, obs_y), na.rm = TRUE)
  xlim <- range(c(pred$x, obs_x), na.rm = TRUE)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  lines(pred$x, pred$pi_hi95, col = "blue", lty = 4, lwd = 2.5)
  lines(pred$x, pred$hi95,    col = "blue", lty = 4, lwd = 1)
  lines(pred$x, pred$mean,    col = "black", lty = 1, lwd = 2.5)
  lines(pred$x, pred$lo95,    col = "forestgreen", lty = 3, lwd = 1)
  lines(pred$x, pred$pi_lo95, col = "forestgreen", lty = 3, lwd = 2.5)
  points(obs_x, obs_y, pch = 1)
  mtext(label, side = 3, line = -1.8, adj = 0.04, cex = 1.15)
}

## Empty cell (pair not fitted, or level not available)
d02_blank <- function() plot.new()

d02_legend_cell <- function() {
  plot.new()
  legend("center", bty = "n", cex = 1.15,
         legend = c("Mean", "Lower PI", "Upper PI", "Lower CI", "Upper CI", "Observations"),
         col    = c("black", "forestgreen", "blue", "forestgreen", "blue", "black"),
         lty    = c(1, 3, 4, 3, 4, NA), lwd = c(2.5, 2.5, 2.5, 1, 1, NA),
         pch    = c(NA, NA, NA, NA, NA, 1))
}

## Row-label cell (stacked layout): level name, rotated
d02_row_label <- function(label) {
  op <- par(mar = c(0, 0, 0, 0), xpd = NA)
  on.exit(par(op))
  plot.new()
  text(0.5, 0.5, label, srt = 90, cex = 1.15, font = 2)
}

d02_page_close <- function(title, ylab) {
  mtext(ylab, side = 2, outer = TRUE, line = 1.8, cex = 1.3)
  mtext(title, side = 3, outer = TRUE, line = 1, cex = 1.8, font = 2)
}

## Simple 2x3 page, used for the reef-only eDNA pages
d02_page_open <- function() {
  par(mfrow = c(2, 3), mar = c(4.5, 4.8, 1.8, 1), oma = c(0, 4.5, 3.5, 0), cex = 1.3)
}
build_page <- function(filename, title, ylab, panels) {
  png(file.path(fig_dir, filename), width = 2100, height = 1450, res = 130)
  d02_page_open()
  for (p in panels) if (is.null(p)) d02_blank() else p()
  d02_page_close(title, ylab)
  dev.off()
  cat("Saved", file.path(fig_dir, filename), "\n")
}

## Stacked page: one row per level, ncol_data predictor columns plus a leading
## row-label column. `cells` is a row-major list: NULL (blank), "LEGEND", or a
## function that draws one panel.
build_stacked_page <- function(filename, title, ylab, row_labels, ncol_data, cells) {
  nr <- length(row_labels)
  stopifnot(length(cells) == nr * ncol_data)
  png(file.path(fig_dir, filename), width = 650 + 500 * ncol_data, height = 480 * nr + 250, res = 130)
  mat <- matrix(seq_len(nr * (ncol_data + 1)), nrow = nr, byrow = TRUE)
  layout(mat, widths = c(0.42, rep(1, ncol_data)))
  par(mar = c(4.4, 4, 1, 0.3), oma = c(0, 4.3, 3.2, 0), cex = 1.25)
  idx <- 1
  for (r in seq_len(nr)) {
    d02_row_label(row_labels[r])
    for (cc in seq_len(ncol_data)) {
      cell <- cells[[idx]]; idx <- idx + 1
      if (is.null(cell)) d02_blank()
      else if (identical(cell, "LEGEND")) d02_legend_cell()
      else cell()
    }
  }
  d02_page_close(title, ylab)
  dev.off()
  cat("Saved", file.path(fig_dir, filename), "\n")
}

# ---- Shared data: predictor densities per level ----

ms_tow_scar <- ms_tow |> left_join(ms_tow_bin |> select(tow_id, has_scar), by = "tow_id") |>
  filter(!is.na(has_scar)) |> mutate(manta_combined_tow = y_manta + has_scar)
dens_surf_tow  <- ms_tow_scar$y_surf / ms_tow_scar$effort_surf
dens_manta_tow <- ms_tow_scar$manta_combined_tow / ms_tow_scar$effort_manta_tow

dens_surf_cull <- ms_cull$y_surf / ms_cull$effort_surf
dens_deep_cull <- md_cull$y_deep / md_cull$effort_deep
dens_manta_ms_cull <- (ms_cull$y_manta + ms_cull$feeding_scars) / ms_cull$effort_manta
dens_manta_md_cull <- (md_cull$y_manta + md_cull$feeding_scars) / md_cull$effort_manta

dens_surf_reef <- ms_reef$y_surf / ms_reef$effort_surf
dens_deep_reef <- md_reef$y_deep / md_reef$effort_deep
dens_manta_ms_reef <- (ms_reef$y_manta + ms_reef$feeding_scars) / ms_reef$effort_manta
dens_manta_md_reef <- (md_reef$y_manta + md_reef$feeding_scars) / md_reef$effort_manta

lvl_labels <- c("Tow level", "Cull-site level", "Reef level")

# ---- Page: predicting manta rate (columns: RS Surface, RS Deep, eDNA conc, eDNA prop) ----
build_stacked_page("cots_calibration_predicting_manta_rate.png", "Predicting Manta Rate",
  "Manta combined COTS+scar density (per tow)", lvl_labels, 4,
  list(
    ## -- Tow level --
    function() {
      post <- load_post("manta_from_surf_combined_tow_cots")
      pred <- predict_negbin_curve(post, seq(0, max(dens_surf_tow, na.rm = TRUE), length.out = 100),
                                    log_offset = mean(log(ms_tow_scar$effort_manta)), eps = 0.01)
      d02_panel(pred, dens_surf_tow, ms_tow_scar$manta_combined_tow, "RS Surface density (COTS/tow)", "RS Surface")
    },
    NULL, NULL, NULL,
    ## -- Cull-site level --
    function() {
      post <- load_post("manta_from_surf_combined_cull_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(dens_surf_cull, na.rm = TRUE), length.out = 100), mean(ms_cull$effort_manta), eps = 0.01)
      d02_panel(pred, dens_surf_cull, ms_cull$manta_combined / ms_cull$effort_manta, "RS Surface density (COTS/tow)", "RS Surface")
    },
    function() {
      post <- load_post("manta_from_deep_combined_cull_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(dens_deep_cull, na.rm = TRUE), length.out = 100), mean(md_cull$effort_manta), eps = 0.01)
      d02_panel(pred, dens_deep_cull, md_cull$manta_combined / md_cull$effort_manta, "RS Deep density (COTS/tow)", "RS Deep")
    },
    NULL, NULL,
    ## -- Reef level --
    function() {
      post <- load_post("manta_from_surf_combined_reef_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(dens_surf_reef, na.rm = TRUE), length.out = 100), mean(ms_reef$effort_manta), eps = 0.01)
      d02_panel(pred, dens_surf_reef, ms_reef$manta_combined / ms_reef$effort_manta, "RS Surface density (COTS/tow)", "RS Surface")
    },
    function() {
      post <- load_post("manta_from_deep_combined_reef_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(dens_deep_reef, na.rm = TRUE), length.out = 100), mean(md_reef$effort_manta), eps = 0.01)
      d02_panel(pred, dens_deep_reef, md_reef$manta_combined / md_reef$effort_manta, "RS Deep density (COTS/tow)", "RS Deep")
    },
    function() {
      post <- load_post("manta_from_eDNAconc_combined_reef")
      pred <- predict_negbin_curve_density(post, seq(0, max(reef_manta$edna_conc_rate_reef), length.out = 100), mean(reef_manta$n_tows), eps = 0.01)
      d02_panel(pred, reef_manta$edna_conc_rate_reef, reef_manta$manta_combined / reef_manta$n_tows, "eDNA concentration (per assay)", "eDNA conc.")
    },
    "LEGEND"
  ))

# ---- Page: predicting manta presence ----
build_stacked_page("cots_calibration_predicting_manta_presence.png", "Predicting Manta Presence",
  "Manta combined presence rate", lvl_labels, 4,
  list(
    ## -- Tow level -- (Binomial, n = 1, cloglog link)
    function() {
      post <- load_post("manta_from_surf_tow_cots")
      pred <- predict_binomial_curve(post, seq(0, max(ms_tow_bin$surf_density, na.rm = TRUE), length.out = 100),
                                      offset = mean(log(ms_tow_bin$tow_dist)), eps = 0.01, link = "cloglog")
      d02_panel(pred, ms_tow_bin$surf_density, ms_tow_bin$has_manta, "RS Surface density (COTS/tow)", "RS Surface")
    },
    NULL, NULL, "LEGEND",
    ## -- Cull-site level --
    function() {
      post <- load_post("manta_from_surf_binary_cull_cots")
      size_rep <- round(median(ms_cull$effort_manta))
      pred <- predict_binomial_curve(post, seq(0, max(dens_surf_cull, na.rm = TRUE), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, dens_surf_cull, ms_cull$manta_combined / ms_cull$effort_manta, "RS Surface density (COTS/tow)", "RS Surface")
    },
    function() {
      post <- load_post("manta_from_deep_binary_cull_cots")
      size_rep <- round(median(md_cull$effort_manta))
      pred <- predict_binomial_curve(post, seq(0, max(dens_deep_cull, na.rm = TRUE), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, dens_deep_cull, md_cull$manta_combined / md_cull$effort_manta, "RS Deep density (COTS/tow)", "RS Deep")
    },
    NULL, NULL,
    ## -- Reef level --
    function() {
      post <- load_post("manta_from_surf_binary_reef_cots")
      size_rep <- round(median(ms_reef$effort_manta))
      pred <- predict_binomial_curve(post, seq(0, max(dens_surf_reef, na.rm = TRUE), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, dens_surf_reef, ms_reef$manta_combined / ms_reef$effort_manta, "RS Surface density (COTS/tow)", "RS Surface")
    },
    function() {
      post <- load_post("manta_from_deep_binary_reef_cots")
      size_rep <- round(median(md_reef$effort_manta))
      pred <- predict_binomial_curve(post, seq(0, max(dens_deep_reef, na.rm = TRUE), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, dens_deep_reef, md_reef$manta_combined / md_reef$effort_manta, "RS Deep density (COTS/tow)", "RS Deep")
    },
    function() {
      post <- load_post("manta_from_eDNAconc_binary_reef")
      size_rep <- round(median(reef_manta$n_tows))
      pred <- predict_binomial_curve(post, seq(0, max(reef_manta$edna_conc_rate_reef), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, reef_manta$edna_conc_rate_reef, reef_manta$manta_combined_rate, "eDNA concentration (per assay)", "eDNA conc.")
    },
    function() {
      post <- load_post("manta_from_eDNAprop_binary_reef")
      size_rep <- round(median(reef_manta$n_tows))
      pred <- predict_binomial_curve(post, seq(0, max(reef_manta$edna_prop_reef), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, reef_manta$edna_prop_reef, reef_manta$manta_combined_rate, "eDNA proportion positive", "eDNA prop.")
    }
  ))

# ---- Page: predicting RS Surface (columns: Manta, RS Deep, eDNA conc, eDNA prop) ----
build_stacked_page("cots_calibration_predicting_surf.png", "Predicting RS Surface Rate",
  "RS Surface COTS density (per tow)", lvl_labels, 4,
  list(
    ## -- Tow level --
    function() {
      post <- load_post("surf_from_manta_tow_cots")
      pred <- predict_negbin_curve(post, seq(0, max(dens_manta_tow, na.rm = TRUE), length.out = 100),
                                    log_offset = mean(log(ms_tow_scar$effort_surf)), eps = 0.01)
      d02_panel(pred, dens_manta_tow, ms_tow_scar$y_surf, "Manta combined density (per tow)", "Manta")
    },
    NULL, NULL, "LEGEND",
    ## -- Cull-site level --
    function() {
      post <- load_post("surf_from_manta_cull_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(dens_manta_ms_cull, na.rm = TRUE), length.out = 100), mean(ms_cull$effort_surf), eps = 0.01)
      d02_panel(pred, dens_manta_ms_cull, ms_cull$y_surf / ms_cull$effort_surf, "Manta combined density (per tow)", "Manta")
    },
    function() {
      post <- load_post("surf_from_deep_cull_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(sd_cull$y_deep / sd_cull$effort_deep, na.rm = TRUE), length.out = 100), mean(sd_cull$effort_surf), eps = 0.01)
      d02_panel(pred, sd_cull$y_deep / sd_cull$effort_deep, sd_cull$y_surf / sd_cull$effort_surf, "RS Deep density (COTS/tow)", "RS Deep")
    },
    NULL, NULL,
    ## -- Reef level --
    function() {
      post <- load_post("surf_from_manta_reef_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(dens_manta_ms_reef, na.rm = TRUE), length.out = 100), mean(ms_reef$effort_surf), eps = 0.01)
      d02_panel(pred, dens_manta_ms_reef, ms_reef$y_surf / ms_reef$effort_surf, "Manta combined density (per tow)", "Manta")
    },
    function() {
      post <- load_post("surf_from_deep_reef_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(sd_reef$y_deep / sd_reef$effort_deep, na.rm = TRUE), length.out = 100), mean(sd_reef$effort_surf), eps = 0.01)
      d02_panel(pred, sd_reef$y_deep / sd_reef$effort_deep, sd_reef$y_surf / sd_reef$effort_surf, "RS Deep density (COTS/tow)", "RS Deep")
    },
    function() {
      ser <- reef_surf |> mutate(effort_surf = n_tows)
      post <- load_post("surf_from_eDNAconc_reef")
      pred <- predict_negbin_curve_density(post, seq(0, max(ser$edna_conc_rate_reef), length.out = 100), mean(ser$effort_surf), eps = 0.01)
      d02_panel(pred, ser$edna_conc_rate_reef, ser$cots_reefscan_surface / ser$effort_surf, "eDNA concentration (per assay)", "eDNA conc.")
    },
    function() {
      sepr <- reef_surf |> mutate(effort_surf = n_tows)
      post <- load_post("surf_from_eDNAprop_reef")
      pred <- predict_negbin_curve_density(post, seq(0, max(sepr$edna_prop_reef), length.out = 100), mean(sepr$effort_surf), eps = 0.01)
      d02_panel(pred, sepr$edna_prop_reef, sepr$cots_reefscan_surface / sepr$effort_surf, "eDNA proportion positive", "eDNA prop.")
    }
  ))

# ---- Page: predicting RS Deep (cull-site and reef level only; columns: Manta, RS Surface) ----
build_stacked_page("cots_calibration_predicting_deep.png", "Predicting RS Deep Rate",
  "RS Deep COTS density (per tow)", c("Cull-site level", "Reef level"), 3,
  list(
    ## -- Cull-site level --
    function() {
      post <- load_post("deep_from_manta_cull_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(dens_manta_md_cull, na.rm = TRUE), length.out = 100), mean(md_cull$effort_deep), eps = 0.01)
      d02_panel(pred, dens_manta_md_cull, md_cull$y_deep / md_cull$effort_deep, "Manta combined density (per tow)", "Manta")
    },
    function() {
      post <- load_post("deep_from_surf_cull_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(sd_cull$y_surf / sd_cull$effort_surf, na.rm = TRUE), length.out = 100), mean(sd_cull$effort_deep), eps = 0.01)
      d02_panel(pred, sd_cull$y_surf / sd_cull$effort_surf, sd_cull$y_deep / sd_cull$effort_deep, "RS Surface density (COTS/tow)", "RS Surface")
    },
    NULL,
    ## -- Reef level --
    function() {
      post <- load_post("deep_from_manta_reef_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(dens_manta_md_reef, na.rm = TRUE), length.out = 100), mean(md_reef$effort_deep), eps = 0.01)
      d02_panel(pred, dens_manta_md_reef, md_reef$y_deep / md_reef$effort_deep, "Manta combined density (per tow)", "Manta")
    },
    function() {
      post <- load_post("deep_from_surf_reef_cots")
      pred <- predict_negbin_curve_density(post, seq(0, max(sd_reef$y_surf / sd_reef$effort_surf, na.rm = TRUE), length.out = 100), mean(sd_reef$effort_deep), eps = 0.01)
      d02_panel(pred, sd_reef$y_surf / sd_reef$effort_surf, sd_reef$y_deep / sd_reef$effort_deep, "RS Surface density (COTS/tow)", "RS Surface")
    },
    "LEGEND"
  ))

# ---- Page: predicting eDNA concentration (reef level) ----
build_page("cots_calibration_predicting_ednaconc_reef.png", "Predicting eDNA Concentration (reef level)",
  "eDNA concentration (per assay)",
  list(
    function() {
      post <- load_post("eDNAconc_from_manta_reef")
      pred <- predict_negbin_curve(post, seq(0, max(reef_manta$manta_density_raw), length.out = 100), log_offset = 0, eps = 0.01)
      d02_panel(pred, reef_manta$manta_density_raw, reef_manta$edna_conc_rate_reef, "Manta scars+observed per tow", "Manta")
    },
    function() {
      post <- load_post("eDNAconc_from_surf_reef")
      pred <- predict_negbin_curve(post, seq(0, max(reef_surf$surf_density), length.out = 100), log_offset = 0, eps = 0.01)
      d02_panel(pred, reef_surf$surf_density, reef_surf$edna_conc_rate_reef, "RS Surface density (COTS/tow)", "RS Surface")
    },
    NULL, NULL,
    function() {
      post <- load_post("eDNAconc_from_eDNAprop_reef")
      pred <- predict_negbin_curve(post, seq(0, max(edna_reef$edna_prop_reef), length.out = 100), log_offset = 0, eps = 0.01)
      d02_panel(pred, edna_reef$edna_prop_reef, edna_reef$edna_conc_rate_reef, "eDNA proportion positive", "eDNA prop.")
    },
    d02_legend_cell
  ))

# ---- Page: predicting eDNA proportion (reef level) ----
build_page("cots_calibration_predicting_ednaprop_reef.png", "Predicting eDNA Proportion (reef level)",
  "eDNA proportion positive",
  list(
    function() {
      post <- load_post("eDNAprop_from_manta_reef")
      size_rep <- round(median(reef_manta$n_trials))
      pred <- predict_binomial_curve(post, seq(0, max(reef_manta$manta_density_raw), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, reef_manta$manta_density_raw, reef_manta$n_success / reef_manta$n_trials, "Manta scars+observed per tow", "Manta")
    },
    function() {
      post <- load_post("eDNAprop_from_surf_reef")
      size_rep <- round(median(reef_surf$n_trials))
      pred <- predict_binomial_curve(post, seq(0, max(reef_surf$surf_density), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, reef_surf$surf_density, reef_surf$n_success / reef_surf$n_trials, "RS Surface density (COTS/tow)", "RS Surface")
    },
    NULL,
    function() {
      post <- load_post("eDNAprop_from_eDNAconc_reef")
      size_rep <- round(median(edna_reef$n_trials))
      pred <- predict_binomial_curve(post, seq(0, max(edna_reef$edna_conc_rate_reef), length.out = 100), eps = 0.01, link = "logit", size = size_rep)
      d02_panel(pred, edna_reef$edna_conc_rate_reef, edna_reef$edna_prop_reef, "eDNA concentration (per assay)", "eDNA conc.")
    },
    NULL, d02_legend_cell
  ))
