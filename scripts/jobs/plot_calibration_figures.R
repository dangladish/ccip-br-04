## plot_calibration_figures.R
##
## Draws one calibration figure per model direction (coral cover, COTS
## cross-method, eDNA) plus two combined PDFs. No models are fitted.
## Inputs:  output/posteriors/*.RDS, data/derived/*.RDS.
## Outputs: output/figures/{coral,cots}_calibration.pdf and individual PNGs
##          in output/figures/{coral,cots}_calibration_individual/.
##
## Data-prep tables are rebuilt by sourcing the shared and data-prep files,
## then each direction's saved posterior is passed to the matching
## predict_*_curve() and drawn with the local plot_row() (the two shared
## files define different plot_calib_base() versions, so neither is used).
## Titles follow "X from Y, Z level" (X = response, Y = predictor).
## Posterior labels and predictor/response construction mirror
## build_shiny_app_data.R.

library(tidyverse)

fig_dir       <- "output/figures"
post_dir      <- "output/posteriors"
coral_ind_dir <- file.path(fig_dir, "coral_calibration_individual")
cots_ind_dir  <- file.path(fig_dir, "cots_calibration_individual")
dir.create(coral_ind_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cots_ind_dir,  showWarnings = FALSE, recursive = TRUE)

# ---- Load canonical data and rebuild model prep tables ----

tow_level       <- readRDS("data/derived/tow_level.RDS")
cull_site_level <- readRDS("data/derived/cull_site_level.RDS")
reef_level      <- readRDS("data/derived/reef_level.RDS")
# needed by edna_calibration_data_prep.R
mt              <- readRDS("data/derived/mt.RDS")
edna_raw        <- readRDS("data/derived/edna_raw.RDS")
reef_xwalk      <- readRDS("data/derived/reef_xwalk.RDS")
xwalk           <- reef_xwalk |> filter(has_edna) |> select(edna_reef, reef_label)

source("R/coral_shared.R")
source("R/cots_shared.R")
source("models/cots/edna_calibration_data_prep.R")
source("models/cots/cots_crossmethod_data_prep.R")

load_post <- function(label) readRDS(file.path(post_dir, paste0(label, ".RDS")))

# ---- Plotting helper ----
# scale100 = TRUE for coral cover (proportions shown as %). Curves carry a 95% CrI
# and a 95% PI (predict_binomial_curve() also returns PI bounds).

plot_row <- function(pred, obs_x, obs_y, xlab, ylab, main, scale100 = FALSE, has_pi = TRUE) {
  mult <- if (scale100) 100 else 1
  cols <- if (has_pi) cbind(pred$mean, pred$lo95, pred$hi95, pred$pi_lo95, pred$pi_hi95) * mult
          else        cbind(pred$mean, pred$lo95, pred$hi95) * mult
  lty  <- if (has_pi) c(1, 2, 2, 3, 3) else c(1, 2, 2)
  col  <- if (has_pi) c("black", "blue", "blue", "grey45", "grey45") else c("black", "blue", "blue")
  x <- pred$x * mult; ox <- obs_x * mult; oy <- obs_y * mult
  ylim <- range(c(cols, oy), na.rm = TRUE)
  xlim <- range(c(x, ox), na.rm = TRUE)
  matplot(x, cols, type = "l", lty = lty, lwd = c(2, rep(1, length(lty) - 1)), col = col,
          xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim)
  points(ox, oy, pch = 1)
  legend("topleft", bty = "n", cex = 0.7,
         legend = c("Mean", "95% CrI", if (has_pi) "95% PI", "Observations"),
         lty = c(1, 2, if (has_pi) 3, NA), pch = c(NA, NA, if (has_pi) NA, 1),
         col = c("black", "blue", if (has_pi) "grey45", "black"), lwd = 2)
}

# Write an individual PNG and redraw the plot into the open combined PDF;
# written_pngs records each PNG written.
written_pngs <- character()
emit <- function(label, ind_dir, plot_call) {
  # Restore the PDF device after dev.off() so the redraw goes to the PDF
  pdf_dev <- dev.cur()
  png(file.path(ind_dir, paste0(label, ".png")), width = 900, height = 650, res = 120)
  plot_call()
  dev.off()
  if (pdf_dev > 1) dev.set(pdf_dev)
  plot_call()
  written_pngs[[length(written_pngs) + 1]] <<- file.path(ind_dir, paste0(label, ".png"))
}

lvl_word <- c(tow = "tow", cull = "cull-site", reef = "reef")

# ---- Coral cover: 3 pairs x 3 levels x 2 directions ----

coral_pairs <- list(
  list(A = "Manta", B = "RS Surface",
       tow  = tow_level |> filter(!is.na(pct_HC_surf), !is.na(perc_hard)) |> select(y_manta = perc_hard, y_surf = pct_HC_surf),
       cull = cull_site_level |> select(y_manta = perc_hard, y_surf = pct_HC_surf) |> filter(!is.na(y_manta), !is.na(y_surf)),
       reef = reef_level |> filter(has_surf, !is.na(mean_HC_manta), !is.na(mean_HC_surf)) |> transmute(y_manta = mean_HC_manta, y_surf = mean_HC_surf),
       colA = "y_manta", colB = "y_surf",
       lab_from_B = c(tow = "manta_from_surf_tow", cull = "manta_from_surf_cull", reef = "manta_from_surf_reef"),
       lab_from_A = c(tow = "surf_from_manta_tow", cull = "surf_from_manta_cull", reef = "surf_from_manta_reef")),
  list(A = "Manta", B = "RS Deep",
       tow  = tow_level |> filter(!is.na(pct_HC_deep), !is.na(perc_hard)) |> select(y_manta = perc_hard, y_deep = pct_HC_deep),
       cull = cull_site_level |> select(y_manta = perc_hard, y_deep = pct_HC_deep) |> filter(!is.na(y_manta), !is.na(y_deep)),
       reef = reef_level |> filter(has_deep, !is.na(mean_HC_manta), !is.na(mean_HC_deep)) |> transmute(y_manta = mean_HC_manta, y_deep = mean_HC_deep),
       colA = "y_manta", colB = "y_deep",
       lab_from_B = c(tow = "manta_from_deep_tow", cull = "manta_from_deep_cull", reef = "manta_from_deep_reef"),
       lab_from_A = c(tow = "deep_from_manta_tow", cull = "deep_from_manta_cull", reef = "deep_from_manta_reef")),
  list(A = "RS Surface", B = "RS Deep",
       tow  = tow_level |> filter(!is.na(pct_HC_surf), !is.na(pct_HC_deep)) |> select(y_surf = pct_HC_surf, y_deep = pct_HC_deep),
       cull = cull_site_level |> select(y_surf = pct_HC_surf, y_deep = pct_HC_deep) |> filter(!is.na(y_surf), !is.na(y_deep)),
       reef = reef_level |> filter(has_surf, has_deep, !is.na(mean_HC_surf), !is.na(mean_HC_deep)) |> transmute(y_surf = mean_HC_surf, y_deep = mean_HC_deep),
       colA = "y_surf", colB = "y_deep",
       lab_from_B = c(tow = "surf_from_deep_tow", cull = "surf_from_deep_cull", reef = "surf_from_deep_reef"),
       lab_from_A = c(tow = "deep_from_surf_tow", cull = "deep_from_surf_cull", reef = "deep_from_surf_reef"))
)

pdf(file.path(fig_dir, "coral_calibration.pdf"), onefile = TRUE)
for (p in coral_pairs) {
  for (lvl in c("tow", "cull", "reef")) {
    d <- p[[lvl]]; xA <- d[[p$colA]]; xB <- d[[p$colB]]
    lw <- lvl_word[[lvl]]

    lab <- p$lab_from_B[[lvl]]
    post <- load_post(lab)
    x_seq <- seq(min(xB), max(xB), length.out = 100)
    emit(lab, coral_ind_dir, function() plot_row(
      predict_beta_curve(post, x_seq), xB, xA,
      paste(p$B, "HC %"), paste(p$A, "HC %"),
      sprintf("%s from %s, %s level", p$A, p$B, lw), scale100 = TRUE))

    lab <- p$lab_from_A[[lvl]]
    post <- load_post(lab)
    x_seq <- seq(min(xA), max(xA), length.out = 100)
    emit(lab, coral_ind_dir, function() plot_row(
      predict_beta_curve(post, x_seq), xA, xB,
      paste(p$A, "HC %"), paste(p$B, "HC %"),
      sprintf("%s from %s, %s level", p$B, p$A, lw), scale100 = TRUE))
  }
}
dev.off()
cat("Coral cover: 18 figures written to", coral_ind_dir, "+ coral_calibration.pdf\n")

# ---- COTS cross-method ----

msmd_pairs <- list(
  list(prefix = "surf", label = "RS Surface", cull = ms_cull, reef = ms_reef,
       y_col = "y_surf", eff_col = "effort_surf"),
  list(prefix = "deep", label = "RS Deep", cull = md_cull, reef = md_reef,
       y_col = "y_deep", eff_col = "effort_deep")
)

pdf(file.path(fig_dir, "cots_calibration.pdf"), onefile = TRUE)

for (p in msmd_pairs) {
  for (lvl in c("cull", "reef")) {
    d <- p[[lvl]]; lw <- lvl_word[[lvl]]
    y_manta <- d$y_manta; eff_manta <- d$effort_manta
    yB <- d[[p$y_col]]; effB <- d[[p$eff_col]]
    combined <- d$manta_combined
    # Manta predictor: (COTS count + scar count) / effort, summed (not the OR indicator)
    dens_B <- yB / effB; dens_manta_comb <- (y_manta + d$feeding_scars) / eff_manta

    lab <- sprintf("manta_from_%s_combined_%s_cots", p$prefix, lvl)
    post <- load_post(lab)
    x_seq <- seq(0, max(dens_B, na.rm = TRUE), length.out = 100)
    eff_manta_rep <- mean(eff_manta)
    emit(lab, cots_ind_dir, function() plot_row(
      predict_negbin_curve_density(post, x_seq, eff_manta_rep, eps = 0.01), dens_B, combined / eff_manta,
      paste(p$label, "density (COTS/tow)"), "Manta combined COTS+scar density (per tow)",
      sprintf("Manta combined density from %s, %s level", p$label, lw)))

    lab <- sprintf("manta_from_%s_binary_%s_cots", p$prefix, lvl)
    post <- load_post(lab)
    size_rep <- round(median(eff_manta))
    emit(lab, cots_ind_dir, function() plot_row(
      predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep), dens_B, combined / eff_manta,
      paste(p$label, "density (COTS/tow)"), "Manta combined presence rate",
      sprintf("Manta combined presence rate from %s, %s level", p$label, lw)))

    lab <- sprintf("%s_from_manta_%s_cots", p$prefix, lvl)
    post <- load_post(lab)
    x_seq2 <- seq(0, max(dens_manta_comb, na.rm = TRUE), length.out = 100)
    effB_rep <- mean(effB)
    emit(lab, cots_ind_dir, function() plot_row(
      predict_negbin_curve_density(post, x_seq2, effB_rep, eps = 0.01), dens_manta_comb, yB / effB,
      "Manta density (combined COTS+scar count / tow)", paste(p$label, "COTS density (per tow)"),
      sprintf("%s density from Manta combined COTS+scar, %s level", p$label, lw)))
  }
}

# RS Surface <-> RS Deep, cull and reef level (NegBin both ways)
for (lvl in c("cull", "reef")) {
  d <- if (lvl == "cull") sd_cull else sd_reef; lw <- lvl_word[[lvl]]
  yA <- d$y_surf; effA <- d$effort_surf; yB <- d$y_deep; effB <- d$effort_deep
  dens_A <- yA / effA; dens_B <- yB / effB

  lab <- sprintf("surf_from_deep_%s_cots", lvl); effA_rep <- mean(effA)
  post <- load_post(lab)
  x_seq <- seq(0, max(dens_B, na.rm = TRUE), length.out = 100)
  emit(lab, cots_ind_dir, function() plot_row(
    predict_negbin_curve_density(post, x_seq, effA_rep, eps = 0.01), dens_B, yA / effA,
    "RS Deep density (COTS/tow)", "RS Surface COTS density (per tow)",
    sprintf("RS Surface density from RS Deep, %s level", lw)))

  lab <- sprintf("deep_from_surf_%s_cots", lvl); effB_rep2 <- mean(effB)
  post <- load_post(lab)
  x_seq <- seq(0, max(dens_A, na.rm = TRUE), length.out = 100)
  emit(lab, cots_ind_dir, function() plot_row(
    predict_negbin_curve_density(post, x_seq, effB_rep2, eps = 0.01), dens_A, yB / effB,
    "RS Surface density (COTS/tow)", "RS Deep COTS density (per tow)",
    sprintf("RS Deep density from RS Surface, %s level", lw)))
}

# Manta <-> RS Surface, tow level. Manta-as-response is binary (cloglog);
# RS Surface-as-response uses the combined manta predictor (y_manta + has_scar).
ms_tow_scar <- ms_tow |> left_join(ms_tow_bin |> select(tow_id, has_scar), by = "tow_id") |> filter(!is.na(has_scar)) |>
  mutate(manta_combined_tow = y_manta + has_scar)

lab <- "manta_from_surf_tow_cots"
post <- load_post(lab)
x_seq <- seq(0, max(ms_tow_bin$surf_density, na.rm = TRUE), length.out = 100)
emit(lab, cots_ind_dir, function() plot_row(
  predict_binomial_curve(post, x_seq, offset = mean(log(ms_tow_bin$tow_dist)), eps = 0.01, link = "cloglog"),
  ms_tow_bin$surf_density, ms_tow_bin$has_manta, "RS Surface density (COTS/tow)", "P(manta tow has COTS or scar)",
  "Manta from RS Surface, tow level"))

dens_manta_tow <- ms_tow_scar$manta_combined_tow / ms_tow_scar$effort_manta_tow
lab <- "surf_from_manta_tow_cots"
post <- load_post(lab)
x_seq <- seq(0, max(dens_manta_tow, na.rm = TRUE), length.out = 100)
emit(lab, cots_ind_dir, function() plot_row(
  predict_negbin_curve(post, x_seq, log_offset = mean(log(ms_tow_scar$effort_surf)), eps = 0.01), dens_manta_tow, ms_tow_scar$y_surf,
  "Manta density (combined COTS+scar count / tow)", "RS Surface COTS count (per tow)",
  "RS Surface from Manta combined COTS+scar, tow level"))

# Remaining tow-level directions. The curve is computed before the PNG device
# opens, inside tryCatch, so a failing direction is skipped with a message.
try_emit <- function(lab, pred_fun, obs_x, obs_y, xlab, ylab, main) {
  pred <- tryCatch(pred_fun(load_post(lab)), error = function(e) { message("SKIPPED ", lab, ": ", conditionMessage(e)); NULL })
  if (is.null(pred)) return(invisible(FALSE))
  emit(lab, cots_ind_dir, function() plot_row(pred, obs_x, obs_y, xlab, ylab, main))
  invisible(TRUE)
}
md_tow_scar <- md_tow |> left_join(md_tow_bin |> select(tow_id, has_scar), by = "tow_id") |> filter(!is.na(has_scar)) |>
  mutate(manta_combined_tow = y_manta + has_scar)
n_explor <- 0L

dens_surf_tow <- ms_tow_scar$y_surf / ms_tow_scar$effort_surf
n_explor <- n_explor + try_emit("manta_from_surf_combined_tow_cots",
  function(post) predict_negbin_curve(post, seq(0, max(dens_surf_tow, na.rm = TRUE), length.out = 100), log_offset = mean(log(ms_tow_scar$effort_manta)), eps = 0.01),
  dens_surf_tow, ms_tow_scar$manta_combined_tow, "RS Surface density (COTS/tow)", "Manta: total COTS + scars observed (per tow)",
  "Manta combined count from RS Surface, tow level")

n_explor <- n_explor + try_emit("manta_from_deep_tow_cots",
  function(post) predict_binomial_curve(post, seq(0, max(md_tow_bin$deep_density, na.rm = TRUE), length.out = 100), offset = mean(log(md_tow_bin$tow_dist)), eps = 0.01, link = "cloglog"),
  md_tow_bin$deep_density, md_tow_bin$has_manta, "RS Deep density (COTS/tow)", "P(manta tow has COTS or scar)",
  "Manta from RS Deep, tow level (exploratory)")

dens_deep_tow <- md_tow_scar$y_deep / md_tow_scar$effort_deep
n_explor <- n_explor + try_emit("manta_from_deep_combined_tow_cots",
  function(post) predict_negbin_curve(post, seq(0, max(dens_deep_tow, na.rm = TRUE), length.out = 100), log_offset = mean(log(md_tow_scar$effort_manta)), eps = 0.01),
  dens_deep_tow, md_tow_scar$manta_combined_tow, "RS Deep density (COTS/tow)", "Manta: total COTS + scars observed (per tow)",
  "Manta combined count from RS Deep, tow level (exploratory)")

dens_manta_tow_d <- md_tow_scar$manta_combined_tow / md_tow_scar$effort_manta_tow
n_explor <- n_explor + try_emit("deep_from_manta_tow_cots",
  function(post) predict_negbin_curve(post, seq(0, max(dens_manta_tow_d, na.rm = TRUE), length.out = 100), log_offset = mean(log(md_tow_scar$effort_deep)), eps = 0.01),
  dens_manta_tow_d, md_tow_scar$y_deep, "Manta density (combined COTS+scar count / tow)", "RS Deep COTS count (per tow)",
  "RS Deep from Manta combined COTS+scar, tow level (exploratory)")

dens_sd_deep <- sd_tow$y_deep / sd_tow$effort_deep; dens_sd_surf <- sd_tow$y_surf / sd_tow$effort_surf
n_explor <- n_explor + try_emit("surf_from_deep_tow_cots",
  function(post) predict_negbin_curve(post, seq(0, max(dens_sd_deep, na.rm = TRUE), length.out = 100), log_offset = mean(log(sd_tow$effort_surf)), eps = 0.01),
  dens_sd_deep, sd_tow$y_surf, "RS Deep density (COTS/tow)", "RS Surface COTS count (per tow)",
  "RS Surface from RS Deep, tow level (exploratory)")

n_explor <- n_explor + try_emit("deep_from_surf_tow_cots",
  function(post) predict_negbin_curve(post, seq(0, max(dens_sd_surf, na.rm = TRUE), length.out = 100), log_offset = mean(log(sd_tow$effort_deep)), eps = 0.01),
  dens_sd_surf, sd_tow$y_deep, "RS Surface density (COTS/tow)", "RS Deep COTS count (per tow)",
  "RS Deep from RS Surface, tow level (exploratory)")

cat("COTS cross-method:", 18 + n_explor, "figures written to", cots_ind_dir, "(expected 24 = 18 + 6 exploratory tow-level; pdf still open)\n")

# ---- eDNA (reef level; written to the same cots_calibration.pdf) ----

reef_surf_eff <- reef_surf |> mutate(effort_surf = n_tows)

lab <- "eDNAconc_from_manta_reef"
post <- load_post(lab)
x_seq <- seq(0, max(reef_manta$manta_density_raw, na.rm = TRUE), length.out = 100)
emit(lab, cots_ind_dir, function() plot_row(
  predict_negbin_curve(post, x_seq, log_offset = 0, eps = 0.01),
  reef_manta$manta_density_raw, reef_manta$edna_conc_rate_reef,
  "Manta scars+observed per tow", "eDNA concentration (mean per assay)", "eDNA concentration from Manta, reef level"))

lab <- "manta_from_eDNAconc_binary_reef"
post <- load_post(lab)
x_seq <- seq(0, max(reef_manta$edna_conc_rate_reef, na.rm = TRUE), length.out = 100)
size_rep <- round(median(reef_manta$n_tows))
emit(lab, cots_ind_dir, function() plot_row(
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep),
  reef_manta$edna_conc_rate_reef, reef_manta$manta_combined_rate,
  "eDNA conc. rate (conc/assay)", "Manta combined presence rate", "Manta combined presence rate from eDNA concentration, reef level"))

# eps must match the model fit
lab <- "eDNAconc_from_surf_reef"
post <- load_post(lab)
x_seq <- seq(0, max(reef_surf$surf_density, na.rm = TRUE), length.out = 100)
emit(lab, cots_ind_dir, function() plot_row(
  predict_negbin_curve(post, x_seq, log_offset = 0, eps = 0.01),
  reef_surf$surf_density, reef_surf$edna_conc_rate_reef,
  "RS Surface density (COTS/tow)", "eDNA concentration (mean per assay)", "eDNA concentration from RS Surface, reef level"))

lab <- "surf_from_eDNAconc_reef"; surf_eff_rep1 <- mean(reef_surf_eff$effort_surf)
post <- load_post(lab)
x_seq <- seq(0, max(reef_surf_eff$edna_conc_rate_reef, na.rm = TRUE), length.out = 100)
emit(lab, cots_ind_dir, function() plot_row(
  predict_negbin_curve_density(post, x_seq, surf_eff_rep1, eps = 0.01),
  reef_surf_eff$edna_conc_rate_reef, reef_surf_eff$cots_reefscan_surface / reef_surf_eff$effort_surf,
  "eDNA conc. rate (conc/assay)", "RS Surface COTS density (per tow)", "RS Surface from eDNA concentration, reef level"))

lab <- "eDNAprop_from_eDNAconc_reef"
post <- load_post(lab)
x_seq <- seq(0, max(edna_reef$edna_conc_rate_reef, na.rm = TRUE), length.out = 100)
size_rep <- round(median(edna_reef$n_trials))
emit(lab, cots_ind_dir, function() plot_row(
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep),
  edna_reef$edna_conc_rate_reef, edna_reef$edna_prop_reef,
  "eDNA conc. rate (conc/assay)", "eDNA proportion positive", "eDNA proportion from eDNA concentration, reef level"))

lab <- "eDNAconc_from_eDNAprop_reef"
post <- load_post(lab)
x_seq <- seq(0, max(edna_reef$edna_prop_reef, na.rm = TRUE), length.out = 100)
emit(lab, cots_ind_dir, function() plot_row(
  predict_negbin_curve(post, x_seq, log_offset = 0, eps = 0.01),
  edna_reef$edna_prop_reef, edna_reef$edna_conc_rate_reef,
  "eDNA proportion positive", "eDNA concentration (mean per assay)", "eDNA concentration from eDNA proportion, reef level"))

lab <- "eDNAprop_from_manta_reef"
post <- load_post(lab)
x_seq <- seq(0, max(reef_manta$manta_density_raw, na.rm = TRUE), length.out = 100)
size_rep <- round(median(reef_manta$n_trials))
emit(lab, cots_ind_dir, function() plot_row(
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep),
  reef_manta$manta_density_raw, reef_manta$n_success / reef_manta$n_trials,
  "Manta scars+observed per tow", "eDNA proportion positive (reef)", "eDNA proportion from Manta, reef level"))

lab <- "manta_from_eDNAprop_binary_reef"
post <- load_post(lab)
x_seq <- seq(0, max(reef_manta$edna_prop_reef, na.rm = TRUE), length.out = 100)
size_rep <- round(median(reef_manta$n_tows))
emit(lab, cots_ind_dir, function() plot_row(
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep),
  reef_manta$edna_prop_reef, reef_manta$manta_combined_rate,
  "eDNA proportion positive", "Manta combined presence rate", "Manta combined presence rate from eDNA proportion, reef level"))

lab <- "eDNAprop_from_surf_reef"
post <- load_post(lab)
x_seq <- seq(0, max(reef_surf$surf_density, na.rm = TRUE), length.out = 100)
size_rep <- round(median(reef_surf$n_trials))
emit(lab, cots_ind_dir, function() plot_row(
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep),
  reef_surf$surf_density, reef_surf$n_success / reef_surf$n_trials,
  "RS Surface density (COTS/tow)", "eDNA proportion positive (reef)", "eDNA proportion from RS Surface, reef level"))

lab <- "surf_from_eDNAprop_reef"; surf_eff_rep2 <- mean(reef_surf_eff$effort_surf)
post <- load_post(lab)
x_seq <- seq(0, max(reef_surf_eff$edna_prop_reef, na.rm = TRUE), length.out = 100)
emit(lab, cots_ind_dir, function() plot_row(
  predict_negbin_curve_density(post, x_seq, surf_eff_rep2, eps = 0.01),
  reef_surf_eff$edna_prop_reef, reef_surf_eff$cots_reefscan_surface / reef_surf_eff$effort_surf,
  "eDNA proportion positive", "RS Surface COTS density (per tow)", "RS Surface from eDNA proportion, reef level"))

dev.off()
cat("eDNA: 10 figures written to", cots_ind_dir, "\n")
cat("\nAll", 46 + n_explor, "figures regenerated (expected 52 = 18 coral + 24 cross-method + 10 eDNA).\n")
cat("Combined PDFs: ", file.path(fig_dir, "coral_calibration.pdf"), ", ",
    file.path(fig_dir, "cots_calibration.pdf"), "\n", sep = "")
