## build_shiny_app_data.R
##
## Precomputes the data bundle for the Shiny calibration explorer:
## prediction curves + observations for every model direction, plus the
## uncertainty comparison tables.
## Inputs:  output/posteriors/*.RDS, data/derived/*.RDS.
## Output:  shiny_app/shiny_data.RDS.
##
## The data-prep tables are rebuilt by sourcing the shared and data-prep
## files (R/, models/cots/data-prep), so they match what the models were
## fit on; no models are fitted here. Re-run whenever the posteriors or
## canonical tables change.
##
## `pred_curves` has one row per x-grid point per model direction:
## coral cover (Beta), COTS cross-method (NegBin/Binomial) and eDNA
## (family_app = "eDNA"; the app UI keys on this exact label).
## `observations` holds the matching scatter data.
## Every direction has a credible interval (lo95/hi95) and a prediction
## interval (pi_lo95/pi_hi95). Binomial PIs use a representative trial
## size (`size_rep`, the median n_tows or n_trials), since trials vary
## by row.
## `is_probability = TRUE` flags presence-rate directions (for y-axis labels).


library(tidyverse)
library(arrow)
library(sf)

app_dir  <- "shiny_app"
post_dir <- "output/posteriors"
dir.create(app_dir, showWarnings = FALSE, recursive = TRUE)

# Seed for the random draws used in prediction intervals (reproducible bundle)
set.seed(3947)

# ---- Load canonical data and rebuild model prep tables ----

tow_level       <- readRDS("data/derived/tow_level.RDS")
cull_site_level <- readRDS("data/derived/cull_site_level.RDS")
reef_level      <- readRDS("data/derived/reef_level.RDS")
mt              <- readRDS("data/derived/mt.RDS")
edna_raw        <- readRDS("data/derived/edna_raw.RDS")
reef_xwalk      <- readRDS("data/derived/reef_xwalk.RDS")
rhis            <- readRDS("data/derived/rhis.RDS")
xwalk           <- reef_xwalk |> filter(has_edna) |> select(edna_reef, reef_label)

source("R/coral_shared.R")
source("R/cots_shared.R")
source("models/cots/edna_calibration_data_prep.R")
source("models/cots/cots_crossmethod_data_prep.R")
source("R/uncertainty_shared.R")


# ---- Helpers: load posterior, build one curve + observation set ----

load_post <- function(label) readRDS(file.path(post_dir, paste0(label, ".RDS")))

curve_row <- function(id_bits, x, y, curve, n, beta, is_probability = FALSE, ...) {
  extra <- list(...)
  row <- curve |>
    mutate(
      id = paste(id_bits, collapse = "||"),
      family_app = id_bits[1], level = id_bits[2], from = id_bits[3], to = id_bits[4],
      n = n, beta_mean = mean(beta), beta_lo95 = quantile(beta, 0.025), beta_hi95 = quantile(beta, 0.975),
      is_probability = is_probability
    )
  for (nm in names(extra)) row[[nm]] <- extra[[nm]]
  list(
    curve = row,
    obs = tibble(id = paste(id_bits, collapse = "||"),
                 family_app = id_bits[1], level = id_bits[2], from = id_bits[3], to = id_bits[4],
                 x = x, y = y)
  )
}

pred_list <- list(); obs_list <- list()
add_result <- function(res) {
  pred_list[[length(pred_list) + 1]] <<- res$curve
  obs_list[[length(obs_list) + 1]]   <<- res$obs
}

# ---- Coral cover: 3 pairs x 3 levels x 2 directions ----
coral_pairs <- list(
  list(A = "Manta tow", B = "ReefScan Surface",
       tow  = tow_level |> filter(!is.na(pct_HC_surf), !is.na(perc_hard)) |> select(y_manta = perc_hard, y_surf = pct_HC_surf),
       cull = cull_site_level |> select(y_manta = perc_hard, y_surf = pct_HC_surf) |> filter(!is.na(y_manta), !is.na(y_surf)),
       reef = reef_level |> filter(has_surf, !is.na(mean_HC_manta), !is.na(mean_HC_surf)) |> transmute(y_manta = mean_HC_manta, y_surf = mean_HC_surf),
       colA = "y_manta", colB = "y_surf",
       lab_from_B = c(tow = "manta_from_surf_tow", cull = "manta_from_surf_cull", reef = "manta_from_surf_reef"),
       lab_from_A = c(tow = "surf_from_manta_tow", cull = "surf_from_manta_cull", reef = "surf_from_manta_reef")),
  list(A = "Manta tow", B = "ReefScan Deep",
       tow  = tow_level |> filter(!is.na(pct_HC_deep), !is.na(perc_hard)) |> select(y_manta = perc_hard, y_deep = pct_HC_deep),
       cull = cull_site_level |> select(y_manta = perc_hard, y_deep = pct_HC_deep) |> filter(!is.na(y_manta), !is.na(y_deep)),
       reef = reef_level |> filter(has_deep, !is.na(mean_HC_manta), !is.na(mean_HC_deep)) |> transmute(y_manta = mean_HC_manta, y_deep = mean_HC_deep),
       colA = "y_manta", colB = "y_deep",
       lab_from_B = c(tow = "manta_from_deep_tow", cull = "manta_from_deep_cull", reef = "manta_from_deep_reef"),
       lab_from_A = c(tow = "deep_from_manta_tow", cull = "deep_from_manta_cull", reef = "deep_from_manta_reef")),
  list(A = "ReefScan Surface", B = "ReefScan Deep",
       tow  = tow_level |> filter(!is.na(pct_HC_surf), !is.na(pct_HC_deep)) |> select(y_surf = pct_HC_surf, y_deep = pct_HC_deep),
       cull = cull_site_level |> select(y_surf = pct_HC_surf, y_deep = pct_HC_deep) |> filter(!is.na(y_surf), !is.na(y_deep)),
       reef = reef_level |> filter(has_surf, has_deep, !is.na(mean_HC_surf), !is.na(mean_HC_deep)) |> transmute(y_surf = mean_HC_surf, y_deep = mean_HC_deep),
       colA = "y_surf", colB = "y_deep",
       lab_from_B = c(tow = "surf_from_deep_tow", cull = "surf_from_deep_cull", reef = "surf_from_deep_reef"),
       lab_from_A = c(tow = "deep_from_surf_tow", cull = "deep_from_surf_cull", reef = "deep_from_surf_reef"))
)
lvl_ui <- c(tow = "Tow", cull = "Cull-site", reef = "Reef")

for (p in coral_pairs) {
  for (lvl in c("tow", "cull", "reef")) {
    d <- p[[lvl]]
    xA <- d[[p$colA]]; xB <- d[[p$colB]]

    post_AfromB <- load_post(p$lab_from_B[[lvl]])
    x_seq <- seq(min(xB), max(xB), length.out = 100)
    add_result(curve_row(
      c("Coral cover", lvl_ui[[lvl]], p$B, p$A), xB, xA,
      predict_beta_curve(post_AfromB, x_seq), length(xA), post_AfromB[, "beta"],
      xlab = paste(p$B, "hard coral cover (%)"), ylab = paste(p$A, "hard coral cover (%)"), pct_axes = TRUE
    ))

    post_BfromA <- load_post(p$lab_from_A[[lvl]])
    x_seq <- seq(min(xA), max(xA), length.out = 100)
    add_result(curve_row(
      c("Coral cover", lvl_ui[[lvl]], p$A, p$B), xA, xB,
      predict_beta_curve(post_BfromA, x_seq), length(xB), post_BfromA[, "beta"],
      xlab = paste(p$A, "hard coral cover (%)"), ylab = paste(p$B, "hard coral cover (%)"), pct_axes = TRUE
    ))
  }
}

# ---- COTS cross-method: manta / RS Surface / RS Deep ----
msmd_pairs <- list(
  list(prefix = "surf", label = "ReefScan Surface", cull = ms_cull, reef = ms_reef,
       y_col = "y_surf", eff_col = "effort_surf"),
  list(prefix = "deep", label = "ReefScan Deep", cull = md_cull, reef = md_reef,
       y_col = "y_deep", eff_col = "effort_deep")
)

for (p in msmd_pairs) {
  for (lvl in c("cull", "reef")) {
    d <- p[[lvl]]
    y_manta <- d$y_manta; eff_manta <- d$effort_manta
    yB <- d[[p$y_col]]; effB <- d[[p$eff_col]]
    combined <- d$manta_combined
    # Manta predictor density: COTS count + scar count, summed (not the OR indicator)
    dens_B <- yB / effB; dens_manta_comb <- (y_manta + d$feeding_scars) / eff_manta

    post <- load_post(sprintf("manta_from_%s_combined_%s_cots", p$prefix, lvl))
    x_seq <- seq(0, max(dens_B, na.rm = TRUE), length.out = 100)
    eff_manta_rep <- mean(eff_manta)
    add_result(curve_row(
      c("COTS cross-method", lvl_ui[[lvl]], p$label, "Manta tow -- combined density"), dens_B, combined / eff_manta,
      predict_negbin_curve_density(post, x_seq, eff_manta_rep, eps = 0.01), length(combined), post[, "beta"],
      xlab = paste(p$label, "density (COTS/tow)"), ylab = "Manta combined COTS+scar density (per tow)"
    ))

    post <- load_post(sprintf("manta_from_%s_binary_%s_cots", p$prefix, lvl))
    size_rep <- round(median(eff_manta))
    add_result(curve_row(
      c("COTS cross-method", lvl_ui[[lvl]], p$label, "Manta tow -- presence rate"), dens_B, combined / eff_manta,
      predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep), length(combined), post[, "beta"], is_probability = TRUE,
      xlab = paste(p$label, "density (COTS/tow)"), ylab = "Manta combined presence rate",
      pi_note = if (lvl == "reef") sprintf("95%% PI: a new reef with %d tows (median observed), incl. between-reef variation", size_rep)
                else sprintf("95%% PI assumes %d tows (median observed) at this cull site (average reef)", size_rep)
    ))

    post <- load_post(sprintf("%s_from_manta_%s_cots", p$prefix, lvl))
    x_seq2 <- seq(0, max(dens_manta_comb, na.rm = TRUE), length.out = 100); effB_rep <- mean(effB)
    add_result(curve_row(
      c("COTS cross-method", lvl_ui[[lvl]], "Manta tow", p$label), dens_manta_comb, yB / effB,
      predict_negbin_curve_density(post, x_seq2, effB_rep, eps = 0.01), length(yB), post[, "beta"],
      xlab = "Manta density (combined COTS+scar count / tow)", ylab = paste(p$label, "COTS density (per tow)")
    ))
  }
}

# RS Surface <-> RS Deep, cull and reef level (NegBin both ways)
for (lvl in c("cull", "reef")) {
  d <- if (lvl == "cull") sd_cull else sd_reef
  yA <- d$y_surf; effA <- d$effort_surf; yB <- d$y_deep; effB <- d$effort_deep
  dens_A <- yA / effA; dens_B <- yB / effB

  post_AfromB <- load_post(sprintf("surf_from_deep_%s_cots", lvl))
  x_seq <- seq(0, max(dens_B, na.rm = TRUE), length.out = 100); effA_rep <- mean(effA)
  add_result(curve_row(
    c("COTS cross-method", lvl_ui[[lvl]], "ReefScan Deep", "ReefScan Surface"), dens_B, yA / effA,
    predict_negbin_curve_density(post_AfromB, x_seq, effA_rep, eps = 0.01), length(yA), post_AfromB[, "beta"],
    xlab = "RS Deep density (COTS/tow)", ylab = "RS Surface COTS density (per tow)"
  ))

  post_BfromA <- load_post(sprintf("deep_from_surf_%s_cots", lvl))
  x_seq <- seq(0, max(dens_A, na.rm = TRUE), length.out = 100); effB_rep2 <- mean(effB)
  add_result(curve_row(
    c("COTS cross-method", lvl_ui[[lvl]], "ReefScan Surface", "ReefScan Deep"), dens_A, yB / effB,
    predict_negbin_curve_density(post_BfromA, x_seq, effB_rep2, eps = 0.01), length(yB), post_BfromA[, "beta"],
    xlab = "RS Surface density (COTS/tow)", ylab = "RS Deep COTS density (per tow)"
  ))
}

# Manta <-> RS Surface, tow level. Manta-as-response is binary (cloglog);
# RS Surface-as-response uses the combined manta predictor (y_manta + has_scar).
ms_tow_scar <- ms_tow |> left_join(ms_tow_bin |> select(tow_id, has_scar), by = "tow_id") |> filter(!is.na(has_scar)) |>
  mutate(manta_combined_tow = y_manta + has_scar)

post <- load_post("manta_from_surf_tow_cots")
x_seq <- seq(0, max(ms_tow_bin$surf_density, na.rm = TRUE), length.out = 100)
add_result(curve_row(
  c("COTS cross-method", "Tow", "ReefScan Surface", "Manta tow"), ms_tow_bin$surf_density, ms_tow_bin$has_manta,
  predict_binomial_curve(post, x_seq, offset = mean(log(ms_tow_bin$tow_dist)), eps = 0.01, link = "cloglog"),
  length(ms_tow_bin$has_manta), post[, "beta"], is_probability = TRUE,
  xlab = "RS Surface density (COTS/tow)", ylab = "P(Manta tow has COTS or scar)"
))

dens_manta_tow <- ms_tow_scar$manta_combined_tow / ms_tow_scar$effort_manta_tow
post <- load_post("surf_from_manta_tow_cots")
x_seq <- seq(0, max(dens_manta_tow, na.rm = TRUE), length.out = 100)
add_result(curve_row(
  c("COTS cross-method", "Tow", "Manta tow", "ReefScan Surface"), dens_manta_tow, ms_tow_scar$y_surf,
  predict_negbin_curve(post, x_seq, log_offset = mean(log(ms_tow_scar$effort_surf)), eps = 0.01), length(ms_tow_scar$y_surf), post[, "beta"],
  xlab = "Manta density (combined COTS+scar count / tow)", ylab = "RS Surface COTS count (per tow)"
))

# Remaining tow-level directions (manta combined count ~ RS Surface, and the
# RS Deep tow models). try_add() skips a direction with a message if its curve fails.
try_add <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) { message("SKIPPED ", label, ": ", conditionMessage(e)); NULL })
  if (!is.null(res)) add_result(res)
}
md_tow_scar <- md_tow |> left_join(md_tow_bin |> select(tow_id, has_scar), by = "tow_id") |> filter(!is.na(has_scar)) |>
  mutate(manta_combined_tow = y_manta + has_scar)

try_add("manta_from_surf_combined_tow_cots", {
  post <- load_post("manta_from_surf_combined_tow_cots")
  dens_surf_tow <- ms_tow_scar$y_surf / ms_tow_scar$effort_surf
  x_seq <- seq(0, max(dens_surf_tow, na.rm = TRUE), length.out = 100)
  curve_row(
    c("COTS cross-method", "Tow", "ReefScan Surface", "Manta tow -- combined count"), dens_surf_tow, ms_tow_scar$manta_combined_tow,
    predict_negbin_curve(post, x_seq, log_offset = mean(log(ms_tow_scar$effort_manta)), eps = 0.01), length(ms_tow_scar$manta_combined_tow), post[, "beta"],
    xlab = "RS Surface density (COTS/tow)", ylab = "Manta: total COTS + scars observed (per tow)"
  )
})

try_add("manta_from_deep_tow_cots", {
  post <- load_post("manta_from_deep_tow_cots")
  x_seq <- seq(0, max(md_tow_bin$deep_density, na.rm = TRUE), length.out = 100)
  curve_row(
    c("COTS cross-method", "Tow", "ReefScan Deep", "Manta tow"), md_tow_bin$deep_density, md_tow_bin$has_manta,
    predict_binomial_curve(post, x_seq, offset = mean(log(md_tow_bin$tow_dist)), eps = 0.01, link = "cloglog"),
    length(md_tow_bin$has_manta), post[, "beta"], is_probability = TRUE,
    xlab = "RS Deep density (COTS/tow)", ylab = "P(Manta tow has COTS or scar)"
  )
})

dens_deep_tow <- md_tow_scar$y_deep / md_tow_scar$effort_deep
try_add("manta_from_deep_combined_tow_cots", {
  post <- load_post("manta_from_deep_combined_tow_cots")
  x_seq <- seq(0, max(dens_deep_tow, na.rm = TRUE), length.out = 100)
  curve_row(
    c("COTS cross-method", "Tow", "ReefScan Deep", "Manta tow -- combined count"), dens_deep_tow, md_tow_scar$manta_combined_tow,
    predict_negbin_curve(post, x_seq, log_offset = mean(log(md_tow_scar$effort_manta)), eps = 0.01), length(md_tow_scar$manta_combined_tow), post[, "beta"],
    xlab = "RS Deep density (COTS/tow)", ylab = "Manta: total COTS + scars observed (per tow)"
  )
})

try_add("deep_from_manta_tow_cots", {
  post <- load_post("deep_from_manta_tow_cots")
  dens_manta_tow_d <- md_tow_scar$manta_combined_tow / md_tow_scar$effort_manta_tow
  x_seq <- seq(0, max(dens_manta_tow_d, na.rm = TRUE), length.out = 100)
  curve_row(
    c("COTS cross-method", "Tow", "Manta tow", "ReefScan Deep"), dens_manta_tow_d, md_tow_scar$y_deep,
    predict_negbin_curve(post, x_seq, log_offset = mean(log(md_tow_scar$effort_deep)), eps = 0.01), length(md_tow_scar$y_deep), post[, "beta"],
    xlab = "Manta density (combined COTS+scar count / tow)", ylab = "RS Deep COTS count (per tow)"
  )
})

try_add("surf_from_deep_tow_cots", {
  post <- load_post("surf_from_deep_tow_cots")
  dens_B <- sd_tow$y_deep / sd_tow$effort_deep
  x_seq <- seq(0, max(dens_B, na.rm = TRUE), length.out = 100)
  curve_row(
    c("COTS cross-method", "Tow", "ReefScan Deep", "ReefScan Surface"), dens_B, sd_tow$y_surf,
    predict_negbin_curve(post, x_seq, log_offset = mean(log(sd_tow$effort_surf)), eps = 0.01), length(sd_tow$y_surf), post[, "beta"],
    xlab = "RS Deep density (COTS/tow)", ylab = "RS Surface COTS count (per tow)"
  )
})

try_add("deep_from_surf_tow_cots", {
  post <- load_post("deep_from_surf_tow_cots")
  dens_A <- sd_tow$y_surf / sd_tow$effort_surf
  x_seq <- seq(0, max(dens_A, na.rm = TRUE), length.out = 100)
  curve_row(
    c("COTS cross-method", "Tow", "ReefScan Surface", "ReefScan Deep"), dens_A, sd_tow$y_deep,
    predict_negbin_curve(post, x_seq, log_offset = mean(log(sd_tow$effort_deep)), eps = 0.01), length(sd_tow$y_deep), post[, "beta"],
    xlab = "RS Surface density (COTS/tow)", ylab = "RS Deep COTS count (per tow)"
  )
})

# ---- eDNA (reef level only, so level = NA) ----
reef_surf_eff <- reef_surf |> mutate(effort_surf = n_tows)

post <- load_post("eDNAconc_from_manta_reef")
x_seq <- seq(0, max(reef_manta$manta_density_raw, na.rm = TRUE), length.out = 100)
add_result(curve_row(c("eDNA", NA_character_, "Manta tow", "eDNA concentration"), reef_manta$manta_density_raw, reef_manta$edna_conc_rate_reef,
  predict_negbin_curve(post, x_seq, log_offset = 0, eps = 0.01), nrow(reef_manta), post[, "beta"],
  xlab = "Manta scars+observed per tow", ylab = "eDNA concentration (mean per assay)"))

post <- load_post("manta_from_eDNAconc_binary_reef")
x_seq <- seq(0, max(reef_manta$edna_conc_rate_reef, na.rm = TRUE), length.out = 100)
size_rep <- round(median(reef_manta$n_tows))
add_result(curve_row(c("eDNA", NA_character_, "eDNA concentration", "Manta tow -- presence rate"), reef_manta$edna_conc_rate_reef, reef_manta$manta_combined_rate,
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep), nrow(reef_manta), post[, "beta"], is_probability = TRUE,
  xlab = "eDNA conc. rate (conc/assay)", ylab = "Manta combined presence rate",
  pi_note = sprintf("95%% PI: a new reef with %d tows (median observed), incl. between-reef variation", size_rep)))

# eps must match the model fit
post <- load_post("eDNAconc_from_surf_reef")
x_seq <- seq(0, max(reef_surf$surf_density, na.rm = TRUE), length.out = 100)
add_result(curve_row(c("eDNA", NA_character_, "ReefScan Surface", "eDNA concentration"), reef_surf$surf_density, reef_surf$edna_conc_rate_reef,
  predict_negbin_curve(post, x_seq, log_offset = 0, eps = 0.01), nrow(reef_surf), post[, "beta"],
  xlab = "RS Surface density (COTS/tow)", ylab = "eDNA concentration (mean per assay)"))

post <- load_post("surf_from_eDNAconc_reef")
x_seq <- seq(0, max(reef_surf_eff$edna_conc_rate_reef, na.rm = TRUE), length.out = 100); surf_eff_rep1 <- mean(reef_surf_eff$effort_surf)
add_result(curve_row(c("eDNA", NA_character_, "eDNA concentration", "ReefScan Surface"), reef_surf_eff$edna_conc_rate_reef, reef_surf_eff$cots_reefscan_surface / reef_surf_eff$effort_surf,
  predict_negbin_curve_density(post, x_seq, surf_eff_rep1, eps = 0.01), nrow(reef_surf_eff), post[, "beta"],
  xlab = "eDNA conc. rate (conc/assay)", ylab = "RS Surface COTS density (per tow)"))

post <- load_post("eDNAprop_from_eDNAconc_reef")
x_seq <- seq(0, max(edna_reef$edna_conc_rate_reef, na.rm = TRUE), length.out = 100)
size_rep <- round(median(edna_reef$n_trials))
add_result(curve_row(c("eDNA", NA_character_, "eDNA concentration", "eDNA proportion"), edna_reef$edna_conc_rate_reef, edna_reef$edna_prop_reef,
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep), nrow(edna_reef), post[, "beta"], is_probability = TRUE,
  xlab = "eDNA conc. rate (conc/assay)", ylab = "eDNA proportion positive",
  pi_note = sprintf("95%% PI: a new reef with %d samples (median observed), incl. between-reef variation", size_rep)))

post <- load_post("eDNAconc_from_eDNAprop_reef")
x_seq <- seq(0, max(edna_reef$edna_prop_reef, na.rm = TRUE), length.out = 100)
add_result(curve_row(c("eDNA", NA_character_, "eDNA proportion", "eDNA concentration"), edna_reef$edna_prop_reef, edna_reef$edna_conc_rate_reef,
  predict_negbin_curve(post, x_seq, log_offset = 0, eps = 0.01), nrow(edna_reef), post[, "beta"],
  xlab = "eDNA proportion positive", ylab = "eDNA concentration (mean per assay)"))

post <- load_post("eDNAprop_from_manta_reef")
x_seq <- seq(0, max(reef_manta$manta_density_raw, na.rm = TRUE), length.out = 100)
size_rep <- round(median(reef_manta$n_trials))
add_result(curve_row(c("eDNA", NA_character_, "Manta tow", "eDNA proportion"), reef_manta$manta_density_raw, reef_manta$n_success / reef_manta$n_trials,
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep), nrow(reef_manta), post[, "beta"], is_probability = TRUE,
  xlab = "Manta scars+observed per tow", ylab = "eDNA proportion positive",
  pi_note = sprintf("95%% PI: a new reef with %d samples (median observed), incl. between-reef variation", size_rep)))

post <- load_post("manta_from_eDNAprop_binary_reef")
x_seq <- seq(0, max(reef_manta$edna_prop_reef, na.rm = TRUE), length.out = 100)
size_rep <- round(median(reef_manta$n_tows))
add_result(curve_row(c("eDNA", NA_character_, "eDNA proportion", "Manta tow -- presence rate"), reef_manta$edna_prop_reef, reef_manta$manta_combined_rate,
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep), nrow(reef_manta), post[, "beta"], is_probability = TRUE,
  xlab = "eDNA proportion positive", ylab = "Manta combined presence rate",
  pi_note = sprintf("95%% PI: a new reef with %d tows (median observed), incl. between-reef variation", size_rep)))

post <- load_post("eDNAprop_from_surf_reef")
x_seq <- seq(0, max(reef_surf$surf_density, na.rm = TRUE), length.out = 100)
size_rep <- round(median(reef_surf$n_trials))
add_result(curve_row(c("eDNA", NA_character_, "ReefScan Surface", "eDNA proportion"), reef_surf$surf_density, reef_surf$n_success / reef_surf$n_trials,
  predict_binomial_curve(post, x_seq, eps = 0.01, link = "logit", size = size_rep), nrow(reef_surf), post[, "beta"], is_probability = TRUE,
  xlab = "RS Surface density (COTS/tow)", ylab = "eDNA proportion positive",
  pi_note = sprintf("95%% PI: a new reef with %d samples (median observed), incl. between-reef variation", size_rep)))

post <- load_post("surf_from_eDNAprop_reef")
x_seq <- seq(0, max(reef_surf_eff$edna_prop_reef, na.rm = TRUE), length.out = 100); surf_eff_rep2 <- mean(reef_surf_eff$effort_surf)
add_result(curve_row(c("eDNA", NA_character_, "eDNA proportion", "ReefScan Surface"), reef_surf_eff$edna_prop_reef, reef_surf_eff$cots_reefscan_surface / reef_surf_eff$effort_surf,
  predict_negbin_curve_density(post, x_seq, surf_eff_rep2, eps = 0.01), nrow(reef_surf_eff), post[, "beta"],
  xlab = "eDNA proportion positive", ylab = "RS Surface COTS density (per tow)"))
pred_curves  <- bind_rows(pred_list)
observations <- bind_rows(obs_list)

cat("\npred_curves:", nrow(pred_curves), "rows across", n_distinct(pred_curves$id), "model directions (expected 52 = 18 coral + 24 COTS cross-method + 10 eDNA)\n")
cat("observations:", nrow(observations), "rows\n")

# ---- Uncertainty layer 1: analytic SE columns on cull-site/reef observations ----
# Attaches the matching `_se` column for each axis where one exists; NA means
# no error bar (tow level, eDNA, density predictors).
se_lookup_cull <- cull_site_level |>
  transmute(perc_hard, perc_hard_se, pct_HC_surf, pct_HC_surf_se, pct_HC_deep, pct_HC_deep_se,
            cots_observed, cots_observed_se, cots_reefscan_surface, cots_reefscan_surface_se,
            cots_reefscan_deep, cots_reefscan_deep_se)
se_lookup_reef <- reef_level |>
  transmute(mean_HC_manta, mean_HC_manta_se, mean_HC_surf, mean_HC_surf_se, mean_HC_deep, mean_HC_deep_se,
            cots_observed, cots_observed_se, cots_reefscan_surface, cots_reefscan_surface_se,
            cots_reefscan_deep, cots_reefscan_deep_se)

# Look up a row's SE by matching its value against the level's aggregate table
attach_se <- function(obs) {
  obs$x_se <- NA_real_
  obs$y_se <- NA_real_
  for (lvl in c("Cull-site", "Reef")) {
    tbl <- if (lvl == "Cull-site") se_lookup_cull else se_lookup_reef
    rows <- which(obs$level == lvl)
    if (length(rows) == 0) next
    for (i in rows) {
      fam <- obs$family_app[i]
      is_coral <- fam == "Coral cover"
      match_se <- function(val, colnames_val, colnames_se) {
        for (k in seq_along(colnames_val)) {
          v <- tbl[[colnames_val[k]]]
          hit <- which(abs(v - val) < 1e-9)
          if (length(hit) > 0) return(tbl[[colnames_se[k]]][hit[1]])
        }
        NA_real_
      }
      if (is_coral) {
        pairs <- list(c("perc_hard", "perc_hard_se"), c("pct_HC_surf", "pct_HC_surf_se"),
                      c("pct_HC_deep", "pct_HC_deep_se"), c("mean_HC_manta", "mean_HC_manta_se"),
                      c("mean_HC_surf", "mean_HC_surf_se"), c("mean_HC_deep", "mean_HC_deep_se"))
        pairs <- Filter(function(pr) pr[1] %in% names(tbl) && pr[2] %in% names(tbl), pairs)
        val_cols <- vapply(pairs, `[`, character(1), 1); se_cols <- vapply(pairs, `[`, character(1), 2)
        obs$x_se[i] <- match_se(obs$x[i], val_cols, se_cols)
        obs$y_se[i] <- match_se(obs$y[i], val_cols, se_cols)
      } else if (fam == "COTS cross-method") {
        val_cols <- c("cots_observed", "cots_reefscan_surface", "cots_reefscan_deep")
        se_cols  <- c("cots_observed_se", "cots_reefscan_surface_se", "cots_reefscan_deep_se")
        obs$y_se[i] <- match_se(obs$y[i], val_cols, se_cols)  # x here is a density ratio, not a raw SE'd quantity
      }
    }
  }
  obs
}
observations <- attach_se(observations)
cat("observations with a usable y_se: ", sum(!is.na(observations$y_se)), " of ", nrow(observations), "\n", sep = "")

# ---- Uncertainty layer 2: per-reef mean +/- SD across sub-reef units ----
# Same design as 03_uncertainty_report.R; logic is in R/uncertainty_shared.R.

rs_deep_raw <- read_csv(
  "data/raw/reefscan_deep/reefcloud-point-classification-GRB-ReefScan-Deep-2026-08-13-filtered.csv",
  show_col_types = FALSE
)

# Reload RS Deep with the voyage exclusion applied (the read above is unfiltered)
rs_deep_raw <- load_rs_deep_raw_3voyage_excluded()

deep_cell_summ_all <- compute_deep_cell_summary(rs_deep_raw, reef_xwalk)

uncertainty_coral <- compute_coral_uncertainty_reef(tow_level, deep_cell_summ_all, rhis, reef_xwalk)

uncertainty_cots <- compute_cots_uncertainty_reef(tow_level, rhis, reef_xwalk)

uncertainty_comparison <- bind_rows(uncertainty_coral, uncertainty_cots)
cat("uncertainty_comparison:", nrow(uncertainty_comparison), "rows,",
    n_distinct(uncertainty_comparison$reef_label), "reefs,",
    n_distinct(uncertainty_comparison$method), "methods\n")

# COTS on the log1p scale: mean/SD computed on log1p(count) per unit, kept as its own table
uncertainty_cots_log1p <- compute_cots_uncertainty_reef(tow_level, rhis, reef_xwalk, log1p = TRUE)

cat("uncertainty_cots_log1p:", nrow(uncertainty_cots_log1p), "rows,",
    n_distinct(uncertainty_cots_log1p$reef_label), "reefs\n")

# ---- Uncertainty layer 2b: pooled across all reefs, one row per method ----
# Pools every sub-reef unit directly (not a mean of reef means).
rhis_scope <- rhis |> filter(reef_label %in% reef_xwalk$reef_label)

uncertainty_pooled <- bind_rows(
  compute_coral_uncertainty_pooled(tow_level, deep_cell_summ_all, rhis_scope),
  compute_cots_uncertainty_pooled(tow_level, rhis_scope)
)
cat("uncertainty_pooled:", nrow(uncertainty_pooled), "rows (one per method x quantity, pooled across all reefs)\n")

uncertainty_pooled_cots_log1p <- compute_cots_uncertainty_pooled(tow_level, rhis_scope, log1p = TRUE)

# ---- Uncertainty layer 2c: pooled by reef ID band (prefix of reef_label, e.g. "21-XXX") ----
uncertainty_band           <- compute_coral_uncertainty_band(tow_level, deep_cell_summ_all, rhis_scope)
uncertainty_cots_band      <- compute_cots_uncertainty_band(tow_level, rhis_scope)
uncertainty_band_comparison <- bind_rows(uncertainty_band, uncertainty_cots_band)
uncertainty_band_cots_log1p <- compute_cots_uncertainty_band(tow_level, rhis_scope, log1p = TRUE)
cat("uncertainty_band_comparison:", nrow(uncertainty_band_comparison), "rows,",
    n_distinct(uncertainty_band_comparison$band), "bands\n")

# ---- Save bundle ----

coral_midpoints <- c("0" = 0, "1-" = 0.025, "1+" = 0.075, "2-" = 0.15, "2+" = 0.25,
                      "3-" = 0.35, "3+" = 0.45, "4-" = 0.55, "4+" = 0.65,
                      "5-" = 0.75, "5+" = 0.875)
coral_cat_breaks <- c(-Inf, 0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, Inf)
stopifnot(length(coral_cat_breaks) == length(coral_midpoints) + 1)

app_data <- list(
  pred_curves  = pred_curves,
  observations = observations,
  uncertainty_comparison = uncertainty_comparison,
  uncertainty_cots_log1p = uncertainty_cots_log1p,
  uncertainty_pooled = uncertainty_pooled,
  uncertainty_pooled_cots_log1p = uncertainty_pooled_cots_log1p,
  uncertainty_band_comparison = uncertainty_band_comparison,
  uncertainty_band_cots_log1p = uncertainty_band_cots_log1p,
  cots_log1p_breaks = c(0, 1, 2, 5, 10, 20, 50),
  coral_cat_mid    = coral_midpoints,
  coral_cat_breaks = coral_cat_breaks,
  built_at = Sys.time(),
  source   = "ccip-br-04"
)

saveRDS(app_data, file.path(app_dir, "shiny_data.RDS"))
cat("\nSaved:", file.path(app_dir, "shiny_data.RDS"), "\n")
