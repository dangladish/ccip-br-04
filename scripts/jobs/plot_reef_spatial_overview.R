## plot_reef_spatial_overview.R
##
## Two spatial overview figures (JPEG, via save_compact_jpeg()):
##   1. reef_spatial_overview_<TARGET_REEF>.jpg: satellite close-up of one reef
##      showing manta tow tracks, ReefScan Surface and Deep points, RHIS surveys
##      and cull site polygons (when present).
##   2. reef_locations_overview.jpg: map of every reef visited on the RJFMP
##      voyages, coloured by voyage; voyages not used in the calibration models
##      are drawn as an X.
## Inputs: data/derived/{mt,reef_xwalk,rhis}.RDS, raw manta tow KML tracks,
## ReefScan Surface/Deep raw points, cull site shapefile.
## Manta tow geometry is the actual GPS track from the voyage KML exports,
## matched to tows by tow_id (tows curve, so a start-end line understates them).
## Change TARGET_REEF to map a different reef. Run from the repository root.

library(tidyverse)
library(sf)
library(arrow)
library(maptiles)
library(tidyterra)
library(ggspatial)
library(magick)

## Default reef; set TARGET_REEF before sourcing to use a different one.
if (!exists("TARGET_REEF", inherits = FALSE)) TARGET_REEF <- "21-094"

exclude_voyages <- c("FMPCRS 2024-05", "FMPCRS 2025-01", "FMPCRS 2026-03")

fig_dir <- "output/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Save as JPEG: basemap imagery is continuous-tone, so JPEG is far smaller than
## PNG with no visible loss. Written via a temporary PNG, then re-encoded.
save_compact_jpeg <- function(plot, path, width, height, dpi = 300, quality = 90) {
  tmp_png <- sub("\\.jpg$", "_tmp.png", path)
  ggsave(tmp_png, plot, width = width, height = height, dpi = dpi)
  magick::image_read(tmp_png) |>
    magick::image_background("white", flatten = TRUE) |>
    magick::image_write(path, format = "jpeg", quality = quality)
  file.remove(tmp_png)
  cat("Saved:", path, "(", round(file.info(path)$size / 1e6, 2), "MB)\n")
}

# ---- Manta tows, reef crosswalk, RHIS (checkpoints) ----
mt         <- readRDS("data/derived/mt.RDS")
reef_xwalk <- readRDS("data/derived/reef_xwalk.RDS")
rhis       <- readRDS("data/derived/rhis.RDS")

stopifnot(TARGET_REEF %in% reef_xwalk$reef_label)
xwalk_row <- reef_xwalk |> filter(reef_label == TARGET_REEF)

mt_target <- mt |> filter(reef_label == TARGET_REEF, !is.na(start_lat), !is.na(end_lat))

# ---- ReefScan Surface raw points ----
## No saved checkpoint at point grain, so re-read and re-apply the voyage exclusion
## (tow_id -> voyage join; survey-date window where tow_id has no match).
## mt_allvoyages is kept for the GBR-wide map below.
mt_allvoyages <- arrow::read_parquet("data/raw/manta_tows/manta_tows.parquet") |>
  mutate(tow_distance_m = as.numeric(tow_distance_m)) |>
  filter(!is.na(tow_distance_m), tow_distance_m > 0)
voyage_windows <- mt_allvoyages |> summarise(start = min(tow_date), end = max(tow_date), .by = voyage)

in_excluded_window <- function(d) {
  Reduce(`|`, lapply(exclude_voyages, function(v) {
    w <- voyage_windows[voyage_windows$voyage == v, ]
    !is.na(d) & d >= w$start & d <= w$end
  }))
}

rs_surface_raw <- arrow::read_parquet("data/raw/reefscan_surface_coral/reefscan_tows.parquet") |>
  select(tow_id, site_reef_name, image_latitude, image_longitude, survey_date)

rs_surf_voyage  <- rs_surface_raw |> left_join(mt_allvoyages |> distinct(tow_id, voyage), by = "tow_id") |> pull(voyage)
rs_surf_exclude <- (!is.na(rs_surf_voyage) & rs_surf_voyage %in% exclude_voyages) |
  (is.na(rs_surf_voyage) & in_excluded_window(rs_surface_raw$survey_date))
rs_surface_raw  <- rs_surface_raw[!rs_surf_exclude, ]

surf_sites  <- str_split(xwalk_row$rs_surf_site, "; ")[[1]]
surf_target <- rs_surface_raw |> filter(site_reef_name %in% surf_sites)

rm(rs_surf_voyage, rs_surf_exclude, rs_surface_raw)

# ---- ReefScan Deep raw points (shared loader applies the voyage exclusion) ----
source("R/uncertainty_shared.R")
rs_deep_raw <- load_rs_deep_raw_3voyage_excluded()

deep_sites  <- str_split(xwalk_row$rs_deep_sites, "; ")[[1]]
deep_target <- rs_deep_raw |> filter(site_reef_name %in% deep_sites)

rm(rs_deep_raw)

# ---- RHIS surveys ----
rhis_target <- rhis |> filter(reef_label == TARGET_REEF, !is.na(latitude), !is.na(longitude))

cat("Reef", TARGET_REEF, "--", nrow(mt_target), "manta tows,", nrow(surf_target), "RS Surface points,",
    nrow(deep_target), "RS Deep points,", nrow(rhis_target), "RHIS surveys\n")

# ---- Manta tow geometry: GPS tracks from KML exports, matched by tow_id ----
read_voyage_kml_tracks <- function(survey_folder) {
  kml_files <- list.files(file.path("data/raw/manta_tows", survey_folder),
                           pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)
  if (length(kml_files) == 0) return(NULL)
  kml_files |>
    map(~ st_read(.x, quiet = TRUE) |> select(tow_id = Name)) |>
    bind_rows() |>
    distinct(tow_id, .keep_all = TRUE)
}

mt_tracks <- unique(mt_target$survey_folder) |> map(read_voyage_kml_tracks) |> bind_rows()
mt_sf <- mt_tracks |> inner_join(mt_target |> select(tow_id, reef_label, voyage, survey_folder), by = "tow_id")

n_missing <- nrow(mt_target) - nrow(mt_sf)
if (n_missing > 0) {
  warning(n_missing, " of ", nrow(mt_target), " manta tow(s) for reef ", TARGET_REEF,
          " had no matching KML track (older voyage export format?) and are omitted from the figure.")
}
rm(read_voyage_kml_tracks, mt_tracks, n_missing)

# ---- sf geometries ----
## Subsample the dense photo-point clouds for faster rendering.
set.seed(5917)
max_pts <- 6000
surf_sf <- surf_target |> slice_sample(n = min(max_pts, nrow(surf_target))) |>
  st_as_sf(coords = c("image_longitude", "image_latitude"), crs = 4326)
deep_sf <- deep_target |> slice_sample(n = min(max_pts, nrow(deep_target))) |>
  st_as_sf(coords = c("image_longitude", "image_latitude"), crs = 4326)
rhis_sf <- rhis_target |> st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# ---- Satellite basemap (Esri World Imagery) ----
pad <- 0.0008  # ~90 m -- small border around the combined point/line extent
all_lon <- c(mt_target$start_lon, mt_target$end_lon, surf_target$image_longitude, deep_target$image_longitude)
all_lat <- c(mt_target$start_lat, mt_target$end_lat, surf_target$image_latitude, deep_target$image_latitude)
bbox_sf <- st_as_sfc(st_bbox(c(
  xmin = min(all_lon, na.rm = TRUE) - pad, ymin = min(all_lat, na.rm = TRUE) - pad,
  xmax = max(all_lon, na.rm = TRUE) + pad, ymax = max(all_lat, na.rm = TRUE) + pad
), crs = 4326))
bb <- st_bbox(bbox_sf)  # reused for coord_sf() limits and image aspect

basemap <- maptiles::get_tiles(bbox_sf, provider = "Esri.WorldImagery", crop = TRUE, zoom = 17)

# ---- Cull site polygons ----
## Filter to the frame first, then dissolve by name; layer omitted if none intersect.
cull_target <- sf::st_read("data/raw/COTS Control Program Cull Sites - External/Cull_Sites_Master.shp", quiet = TRUE) |>
  sf::st_zm() |>
  sf::st_transform(4326) |>
  sf::st_filter(bbox_sf, .predicate = st_intersects)

has_cull <- nrow(cull_target) > 0
if (has_cull) {
  cull_target <- cull_target |>
    group_by(Name) |>
    summarise(.groups = "drop") |>
    sf::st_make_valid() |>
    rename(cull_site = Name)
}
cat(nrow(cull_target), "cull site polygon(s) found for reef", TARGET_REEF, "\n")

# ---- Reef close-up: plot ----
method_levels <- c("Manta tow", "ReefScan Surface", "ReefScan Deep", if (has_cull) "Cull site boundary")
method_cols   <- c("Manta tow" = "red", "ReefScan Surface" = "gold", "ReefScan Deep" = "magenta",
                    if (has_cull) c("Cull site boundary" = "white"))
is_cull_entry <- method_levels == "Cull site boundary"

reef_spatial_overview <- ggplot() +
  geom_spatraster_rgb(data = basemap, maxcell = 2e6)

if (has_cull) {
  # key_glyph = "path" gives line-style legend keys instead of filled boxes
  reef_spatial_overview <- reef_spatial_overview +
    geom_sf(data = cull_target, aes(color = "Cull site boundary"), fill = NA, linewidth = 0.8,
            linetype = "dashed", key_glyph = "path")
}

reef_spatial_overview <- reef_spatial_overview +
  # key_glyph = "path" so the thick override.aes line keys show for point layers
  geom_sf(data = surf_sf, aes(color = "ReefScan Surface"), size = 0.6, alpha = 0.5, key_glyph = "path") +
  geom_sf(data = deep_sf, aes(color = "ReefScan Deep"), size = 0.6, alpha = 0.5, key_glyph = "path") +
  geom_sf(data = mt_sf, aes(color = "Manta tow"), linewidth = 1.3, key_glyph = "path") +
  geom_sf(data = rhis_sf, aes(fill = "RHIS"), shape = 23, size = 3, color = "black", stroke = 0.8) +
  scale_color_manual(name = NULL, values = method_cols, breaks = method_levels) +
  scale_fill_manual(name = NULL, values = c("RHIS" = "cyan"), breaks = "RHIS") +
  guides(
    # Thick line swatches for all methods; dashed for cull site boundary
    color = guide_legend(order = 1, override.aes = list(
      linetype = ifelse(is_cull_entry, "dashed", "solid"),
      linewidth = ifelse(is_cull_entry, 1, 2.2)
    )),
    fill  = guide_legend(order = 2, override.aes = list(shape = 23, size = 3))
  ) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
  annotation_scale(location = "bl", text_col = "white", line_col = "white", height = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),
                          style = north_arrow_minimal(text_col = "white", line_col = "white")) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(size = 7, color = "grey40")) +
  labs(title = paste0("Reef ", TARGET_REEF, ": sampling footprint by method"),
       caption = "Imagery \u00a9 Esri, Maxar, Earthstar Geographics")

reef_spatial_overview

# ---- Reef close-up: save (height from the bbox aspect ratio plus space for legend/caption) ----
asp <- as.numeric((bb["ymax"] - bb["ymin"]) / ((bb["xmax"] - bb["xmin"]) * cos(mean(bb[c("ymin", "ymax")]) * pi / 180)))
out_w <- 8
out_h <- out_w * asp + 1.3

out_file <- file.path(fig_dir, paste0("reef_spatial_overview_", TARGET_REEF, ".jpg"))
save_compact_jpeg(reef_spatial_overview, out_file, width = out_w, height = out_h)

# ---- GBR-wide reef locations by voyage ----
## Uses mt_allvoyages (all voyages). Each reef is placed at the mean of its tows'
## start/end coordinates.
reef_locations <- mt_allvoyages |>
  rename(reef_label = `Reef label`) |>
  group_by(reef_label, voyage, survey_folder) |>
  summarise(
    centroid_lat = mean(c(start_lat, end_lat), na.rm = TRUE),
    centroid_lon = mean(c(start_lon, end_lon), na.rm = TRUE),
    .groups = "drop"
  ) |>
  # voyage_status drives the point shape in the plot (X = not used in calibration)
  mutate(voyage_status = if_else(voyage %in% exclude_voyages, "Excluded from calibration", "Used in calibration"))
cat(nrow(reef_locations), "reefs across", n_distinct(reef_locations$voyage), "voyages in the GBR-wide overview\n")

# Legend label "voyage (region, n reefs)"; region parsed from survey_folder,
# falling back to the voyage code if the name does not match.
voyage_region <- reef_locations |>
  count(voyage, survey_folder, name = "n_reefs") |>
  mutate(
    region = str_match(survey_folder, "^QPWS_FMPCRS_[0-9]{4}-[0-9]{2}_(.+)_[0-9]{8}$")[, 2],
    voyage_label = if_else(is.na(region), as.character(voyage),
                            paste0(voyage, " (", region, ", ", n_reefs, " reefs)"))
  )
reef_locations <- reef_locations |> left_join(voyage_region |> select(voyage, voyage_label), by = "voyage")

reef_locations_sf <- reef_locations |> st_as_sf(coords = c("centroid_lon", "centroid_lat"), crs = 4326)

gbr_pad <- 0.6  # degrees
gbr_bbox_sf <- st_as_sfc(st_bbox(c(
  xmin = min(reef_locations$centroid_lon) - gbr_pad, ymin = min(reef_locations$centroid_lat) - gbr_pad,
  xmax = max(reef_locations$centroid_lon) + gbr_pad, ymax = max(reef_locations$centroid_lat) + gbr_pad
), crs = 4326))
gbr_bb <- st_bbox(gbr_bbox_sf)

gbr_basemap <- maptiles::get_tiles(gbr_bbox_sf, provider = "Esri.WorldImagery", crop = TRUE, zoom = 7)

# ---- GBR-wide map: plot and save ----
reef_locations_overview <- ggplot() +
  geom_spatraster_rgb(data = gbr_basemap, maxcell = 2e6) +
  # shape marks calibration use; colour marks voyage (an X has no fill)
  geom_sf(data = reef_locations_sf, aes(color = voyage_label, shape = voyage_status), size = 3, stroke = 1.3) +
  scale_color_brewer(name = "Voyage", palette = "Set1") +
  scale_shape_manual(name = NULL, values = c("Used in calibration" = 16, "Excluded from calibration" = 4)) +
  coord_sf(xlim = c(gbr_bb["xmin"], gbr_bb["xmax"]), ylim = c(gbr_bb["ymin"], gbr_bb["ymax"]), expand = FALSE) +
  annotation_scale(location = "bl", text_col = "white", line_col = "white") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm"),
                          style = north_arrow_minimal(text_col = "white", line_col = "white")) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 9),
        plot.caption = element_text(size = 7, color = "grey40")) +
  labs(caption = "Imagery \u00a9 Esri, Maxar, Earthstar Geographics")

reef_locations_overview

gbr_asp   <- as.numeric((gbr_bb["ymax"] - gbr_bb["ymin"]) / ((gbr_bb["xmax"] - gbr_bb["xmin"]) * cos(mean(gbr_bb[c("ymin", "ymax")]) * pi / 180)))
gbr_out_w <- 9
gbr_out_h <- (gbr_out_w - 2.3) * gbr_asp + 0.4  # legend width share; caption headroom

gbr_out_file <- file.path(fig_dir, "reef_locations_overview.jpg")
save_compact_jpeg(reef_locations_overview, gbr_out_file, width = gbr_out_w, height = gbr_out_h)
