# CCIP-BR-04: calibrating RJFMP coral cover and COTS monitoring methods

Analysis code for CCIP-BR-04, which compares monitoring methods used by the
Reef Joint Field Management Program (RJFMP) on the Great Barrier Reef:
manta tow, ReefScan Surface, ReefScan Deep, eDNA, and RHIS (descriptive only).
Bayesian GLMMs (Beta for coral cover; negative binomial and binomial for
Crown-of-Thorns Starfish, COTS) calibrate each method against the others at
manta-tow, cull-site and reef scale, with uncertainty. The work builds on the
CCIP D-02 calibration framework. Methods and results are described in the
final report.

## Repository layout

```
scripts/
  00a_build_cotspotter_detections.R   COTSpotter detections -> tow/cull-site/reef counts
  00_read_and_build_data.R            raw sources -> tow_level / cull_site_level / reef_level
  01_coral_calibration.R              coral cover models (driver)
  02_cots_calibration.R               COTS models (driver)
  03_uncertainty_report.R             uncertainty comparison figures
  04_edna_subreef_necklace.R          sub-reef eDNA necklace figures (24 reefs)
  jobs/                               figure scripts and the Shiny data bundle
R/                                    shared helpers
models/coral/, models/cots/           one file per model pair and spatial level
shiny_app/                            interactive calibration-curve explorer
data/raw/, data/derived/              inputs and derived tables (not in the repository)
output/posteriors/, output/figures/   model outputs (not in the repository)
```

## Running the pipeline

Open `ccip-br-04.Rproj` (scripts use paths relative to the repository root),
place the inputs under `data/raw/`, then run in order:

```r
source("scripts/00a_build_cotspotter_detections.R")
source("scripts/00_read_and_build_data.R")
source("scripts/01_coral_calibration.R")
source("scripts/02_cots_calibration.R")
source("scripts/03_uncertainty_report.R")
source("scripts/04_edna_subreef_necklace.R")   # internet needed on first run (map tiles)
```

The figure scripts then run in any order:

```r
source("scripts/jobs/plot_calibration_figures.R")      # one figure per model + combined PDFs
source("scripts/jobs/plot_report_figures.R")           # combined calibration figures
source("scripts/jobs/plot_cots_calibration_pages.R")   # COTS calibration pages
source("scripts/jobs/plot_reef_spatial_overview.R")    # reef map and example reef (KML tracks; internet for tiles)
source("scripts/jobs/build_shiny_app_data.R")          # Shiny data bundle
```

**Inputs under `data/raw/`:** `manta_tows/` (`manta_tows.parquet`,
`260529-COTS-Manta-Cull-RHIS.xlsx` read from sheet `RHIS`, `manta_scar_counts.csv` (per-tow feeding-scar
counts from the workbook's `Manta` sheet, produced by `scripts/jobs/extract_manta_scar_counts.pl`; run it
only if you start from the full workbook), and the voyage `*.kml` tracks used only for the example-reef map); `reefscan_surface_coral/reefscan_tows.parquet`;
`reefscan_deep/` (ReefScan Deep point-classification CSV); `reefscan_cots/csv_files/*.csv`
(COTSpotter exports); `eDNA/eDNA_BR_02.csv` and `eDNA/RJFMP25_26 Results.csv`;
`COTS Control Program Cull Sites - External/Cull_Sites_Master.shp`.

**Requirements:** R with tidyverse, arrow, sf, readxl, brms (C++ toolchain
for Stan), rjags (JAGS installed), coda, patchwork, maptiles, tidyterra,
ggspatial, magick, shiny and bslib. Fitting dominates run time: `01_` takes
about 15 minutes and `02_` 20 minutes or more.

## Figures in the final report

Written to `output/figures/`. Figure numbers follow the final report.

| Figure | Script | Output |
|---|---|---|
| 1, 2 | `scripts/jobs/plot_reef_spatial_overview.R` | `reef_locations_overview.jpg`, `reef_spatial_overview_21-094.jpg` |
| 3 | `scripts/03_uncertainty_report.R` | `gbr_pooled_combined.png` |
| 4-6 | `scripts/jobs/plot_report_figures.R` | `report_coral_manta_surface_all_levels.png`, `report_cots_crossmethod_cull_all_pairs.png`, `report_edna_reef_pairs.png` |
| B1-B3 | `scripts/03_uncertainty_report.R` | `reef_*_all_reefs*.png` (three files) |
| C1-C6 | `scripts/jobs/plot_cots_calibration_pages.R` | `cots_calibration_predicting_*.png` |
| D1-D24 | `scripts/04_edna_subreef_necklace.R` | `edna_necklace/necklace_<reef>.png` |

## Shiny app

`shiny_app/shiny_data.RDS` is precomputed, so the app runs without the raw
data or posteriors:

```r
shiny::runApp("shiny_app")
```

Re-create the bundle with `scripts/jobs/build_shiny_app_data.R`.

## Data

The raw data and fitted posteriors are not stored in this repository
(`data/` and `output/` are git-ignored placeholders). The input data are
published as a zip file whose contents unpack to `data/raw/`: unzip it into the
repository root and the pipeline above runs as written.

**Data location: TBD.** *(To be added once the data are archived.)*

## Licence and citation

GNU General Public License, version 2 or (at your option) any later version;
see [LICENSE](LICENSE). The raw input data are not covered by this licence.
To cite this repository, see [CITATION.cff](CITATION.cff).
