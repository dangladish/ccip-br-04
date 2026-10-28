## 02_cots_calibration.R
##
## Driver for the COTS calibration models (CCIP D-02 approach, fitted with
## brms/Stan). One file per method pair + spatial level under models/cots/.
##  - eDNA models: reef level only (eDNA has no tow-level grain); eDNA vs manta
##    and eDNA vs RS Surface.
##  - Cross-method models: manta <-> RS Surface <-> RS Deep at tow, cull-site
##    and reef level.
## Inputs:  data/derived/*.RDS (from 00_read_and_build_data.R)
## Outputs: output/posteriors/*.RDS, output/figures/cots_calibration.pdf
## Model design and predictor construction: see the header of each model file.

library(tidyverse)
library(arrow)
library(sf)
library(brms)
library(coda)

dir.create("output/posteriors", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

tow_level       <- readRDS("data/derived/tow_level.RDS")
cull_site_level <- readRDS("data/derived/cull_site_level.RDS")
reef_level      <- readRDS("data/derived/reef_level.RDS")
mt              <- readRDS("data/derived/mt.RDS")
edna_raw        <- readRDS("data/derived/edna_raw.RDS")
reef_xwalk      <- readRDS("data/derived/reef_xwalk.RDS")
xwalk           <- reef_xwalk |> filter(has_edna) |> select(edna_reef, reef_label)

doPlot <- TRUE
if (doPlot) pdf("output/figures/cots_calibration.pdf", onefile = TRUE)

# Shared helpers, priors and fit wrappers; sourced into this environment (not
# local) so the per-model files below can see them.
source("R/cots_shared.R")

# ---- eDNA models (reef level) ----

# Data-prep files are not sourced locally: they build tables the model files read.
source("models/cots/edna_calibration_data_prep.R")

# Per-model files are sourced with local = TRUE so their temporaries (e.g. fit
# objects) are discarded when each file finishes.
edna_model_files <- c(
  "edna_conc_manta.R",
  "edna_conc_surf.R",
  "edna_conc_prop.R",
  "edna_prop_manta.R",
  "edna_prop_surf.R"
)
for (f in edna_model_files) {
  cat("\n---- models/cots/", f, " ----\n", sep = "")
  source(file.path("models/cots", f), local = TRUE)
}

# Free the eDNA tables before building the larger cross-method ones.
rm(edna_reef, reef_manta, reef_surf, mt, edna_raw, reef_xwalk, xwalk, manta_cots_per_ha_factor)

# ---- Cross-method models (manta <-> RS Surface <-> RS Deep) ----

# Builds the tow/cull/reef tables every cross-method model file reads.
source("models/cots/cots_crossmethod_data_prep.R")

crossmethod_model_files <- c(
  # Tow level
  "manta_surf_tow_cots.R",
  "manta_deep_tow_cots.R",
  "surf_deep_tow_cots.R",

  # Cull-site level
  "manta_surf_cull_cots.R",
  "manta_deep_cull_cots.R",
  "surf_deep_cull_cots.R",

  # Reef level (no group term)
  "manta_surf_reef_cots.R",
  "manta_deep_reef_cots.R",
  "surf_deep_reef_cots.R"
)
for (f in crossmethod_model_files) {
  cat("\n---- models/cots/", f, " ----\n", sep = "")
  source(file.path("models/cots", f), local = TRUE)
}

if (doPlot) dev.off()

cat("\n02_cots_calibration.R complete (24 cross-method models incl. exploratory tow-level RS Deep + all 5 eDNA files).\n")
cat("Posteriors in output/posteriors/, plots in output/figures/cots_calibration.pdf\n")
