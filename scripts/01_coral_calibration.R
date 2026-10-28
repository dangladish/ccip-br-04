## 01_coral_calibration.R
##
## Driver for the pairwise coral cover calibration models (CCIP D-02 approach):
## manta <-> RS Surface, manta <-> RS Deep, RS Surface <-> RS Deep, both
## directions, at tow, cull-site and reef level. Each pair + level lives in its
## own file under models/coral/ and is sourced in order below.
## Inputs:  data/derived/{tow,cull_site,reef}_level.RDS (from 00_read_and_build_data.R)
## Outputs: output/posteriors/*.RDS, output/figures/coral_calibration.pdf
## Engine:  brms/Stan for tow-level and the large cull-site pair (too slow in
##          JAGS); JAGS for reef-level and small cull-site pairs.

library(tidyverse)
library(rjags)
library(coda)
library(brms)

dir.create("output/posteriors", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

tow_level       <- readRDS("data/derived/tow_level.RDS")
cull_site_level <- readRDS("data/derived/cull_site_level.RDS")
reef_level      <- readRDS("data/derived/reef_level.RDS")

doPlot <- TRUE
if (doPlot) pdf("output/figures/coral_calibration.pdf", onefile = TRUE)

# Shared helpers; sourced into this environment (not local) so the per-model
# files below can see them.
source("R/coral_shared.R")

# Per-model files are sourced with local = TRUE so their temporaries (e.g. fit
# objects) are discarded when each file finishes.
model_files <- c(
  # ---- Tow level (brms/Stan) ----
  "manta_surface_tow.R",
  "manta_deep_tow.R",
  "surface_deep_tow.R",

  # ---- Reef level (JAGS) ----
  "manta_surface_reef.R",
  "manta_deep_reef.R",
  "surface_deep_reef.R",

  # ---- Cull-site level (manta <-> surface: brms; others: JAGS, no reef effect) ----
  "manta_surface_cull.R",
  "manta_deep_cull.R",
  "surface_deep_cull.R"
)

for (f in model_files) {
  cat("\n---- models/coral/", f, " ----\n", sep = "")
  source(file.path("models/coral", f), local = TRUE)
}

if (doPlot) dev.off()

cat("\nAll 18 coral calibration models fit. Posteriors in output/posteriors/, plots in output/figures/coral_calibration.pdf\n")
