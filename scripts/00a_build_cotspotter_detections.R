## 00a_build_cotspotter_detections.R
##
## Ingests COTSpotter detection exports (downloaded manually from
## https://cots.spot-lab.org; the API requires authentication) and builds
## detection-level plus tow / cull-site / reef count tables for BOTH
## ReefScan platforms: Deep (AIMS_RSD) and Surface (AIMS_RSS).
##
## Sourced by 00_read_and_build_data.R after its ingest step; also runnable on its own once
## data/derived/ingest/manta_tows.parquet exists.
##
## Inputs
##   data/raw/reefscan_cots/csv_files/*.csv
##       One file per trip. Columns: trip_number, trip_name, site_name,
##       survey_name, platform, detection_time ("YYYY-MM-DD HH:MM:SS AEST"),
##       latitude, longitude, depth_m, cots_video_sequence_length,
##       ML_model_confidence, certified_by, certified_at.
##       `platform` is AIMS_RSD (deep) or AIMS_RSS (surface).
##   data/derived/ingest/manta_tows.parquet      (tow geometry + reef labels)
##   data/raw/COTS Control Program Cull Sites - External/Cull_Sites_Master.shp
##
## Reef resolution: `site_name` mostly matches the manta `reef_label`
## directly. A few names need overriding (case differences, or a common name
## used instead of a reef code); see `site_name_overrides` below. Detections
## whose reef is not among the manta reefs are dropped with a warning.
##
## Spatial assignment (same method for both platforms)
##   - Cull site: point-in-polygon (st_within) against Cull_Sites_Master.shp;
##     NA if the detection falls in no cull-site polygon.
##   - Manta tow: nearest tow by point-to-line distance, restricted to tows
##     with the same reef_label (tow line = straight segment between the
##     tow's start and end coordinates; no distance cutoff). Approximate: the
##     closest tow line, not necessarily the tow the detection was made in.
##
## Data quality: detections with latitude == longitude == 0 are corrupt
## coordinates; they are flagged (`invalid_coordinates`) and excluded from
## spatial assignment and from all count tables.
##
## Outputs (data/derived/)
##   cotspotter_{deep,surface}_detections.csv      detection level
##   cotspotter_{deep,surface}_tow_counts.csv      COTS count per tow_id
##   cotspotter_{deep,surface}_cull_site_counts.csv  COTS count per cull site
##   cotspotter_{deep,surface}_reef_counts.csv     COTS count per reef_label

library(tidyverse)
library(sf)
library(arrow)

out_dir <- "data/derived"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Nearest-tow helpers ----

# Tow segments as linestrings in EPSG:3577 (GDA94 Australian Albers, metres).
build_manta_lines <- function(mt) {
  mt |>
    filter(!is.na(start_lat), !is.na(end_lat)) |>
    rowwise() |>
    mutate(geometry = st_sfc(st_linestring(
      matrix(c(start_lon, end_lon, start_lat, end_lat), ncol = 2)
    ))) |>
    ungroup() |>
    st_as_sf(crs = 4326) |>
    st_transform(3577) |>
    select(tow_id, reef_label, tow_date, geometry)
}

# Nearest tow line (same reef only) for each point. Points beyond `max_dist_m`
# of every candidate tow, or with NA / (0, 0) coordinates or NA reef, get NA.
# `ref_date` only feeds the reported date_diff_days; it is not used to match.
assign_nearest_tow <- function(lon, lat, reef_label, manta_lines,
                               max_dist_m = Inf, ref_date = NULL) {
  stopifnot(length(lon) == length(lat), length(lon) == length(reef_label))
  n <- length(lon)

  out <- tibble(
    tow_id         = rep(NA_character_, n),
    dist_to_tow_m  = rep(NA_real_, n),
    tow_date       = rep(as.Date(NA), n),
    date_diff_days = rep(NA_real_, n)
  )

  ok <- !is.na(lon) & !is.na(lat) & !(lon == 0 & lat == 0) & !is.na(reef_label)
  if (!any(ok)) return(out)

  idx_ok <- which(ok)
  pts <- st_as_sf(
    tibble(.row = idx_ok, lon = lon[idx_ok], lat = lat[idx_ok], reef_label = reef_label[idx_ok]),
    coords = c("lon", "lat"), crs = 4326, remove = FALSE
  ) |> st_transform(3577)

  # One distance matrix per reef rather than one st_distance() call per point
  for (rl in unique(pts$reef_label)) {
    cand <- manta_lines[manta_lines$reef_label == rl, ]
    if (nrow(cand) == 0) next
    grp <- pts[pts$reef_label == rl, ]

    dmat <- matrix(as.numeric(st_distance(grp, cand)), nrow = nrow(grp))
    j    <- apply(dmat, 1, which.min)
    dmin <- dmat[cbind(seq_len(nrow(grp)), j)]

    keep <- dmin <= max_dist_m
    rows <- grp$.row[keep]

    out$tow_id[rows]        <- cand$tow_id[j[keep]]
    out$dist_to_tow_m[rows] <- dmin[keep]
    out$tow_date[rows]      <- cand$tow_date[j[keep]]
    if (!is.null(ref_date)) {
      out$date_diff_days[rows] <- as.numeric(ref_date[grp$.row[keep]] - cand$tow_date[j[keep]])
    }
  }
  out
}

# ---- Inputs ----

# All voyages are used here; any voyage exclusion happens in
# 00_read_and_build_data.R. Tows without a usable distance are dropped.
manta <- read_parquet("data/derived/ingest/manta_tows.parquet") |>
  rename(reef_label = `Reef label`) |>
  mutate(tow_distance_m = as.numeric(tow_distance_m)) |>
  filter(!is.na(tow_distance_m), tow_distance_m > 0)

cull_sf <- st_read(
  "data/raw/COTS Control Program Cull Sites - External/Cull_Sites_Master.shp",
  quiet = TRUE
) |> st_transform(4326)

manta_lines <- build_manta_lines(manta)

site_name_overrides <- c(
  "21-250A" = "21-250a", "Lavers" = "21-235", "Mystery" = "21-258",
  "Shell" = "19-028", "Stanley" = "19-045", "Stucco" = "19-141"
)

cots_col_types <- cols(
  trip_number = col_character(), trip_name = col_character(), site_name = col_character(),
  survey_name = col_character(), platform = col_character(), detection_time = col_character(),
  latitude = col_double(), longitude = col_double(), depth_m = col_double(),
  cots_video_sequence_length = col_double(), ML_model_confidence = col_double(),
  certified_by = col_character(), certified_at = col_character()
)
files <- list.files("data/raw/reefscan_cots/csv_files", full.names = TRUE, pattern = "\\.csv$")
cots_raw <- map_dfr(files, read_csv, col_types = cots_col_types)

# ---- Process one platform -> detection table + 3 count tables ----

process_cotspotter_platform <- function(cots_raw, platform, label,
                                        site_name_overrides, manta,
                                        cull_sf, manta_lines) {
  dat <- cots_raw |>
    filter(platform == !!platform) |>
    mutate(
      detection_id = row_number(),
      detection_dttm = as.POSIXct(str_remove(detection_time, " AEST$"), tz = "Australia/Brisbane"),
      detection_date = as.Date(detection_dttm),
      reef_label = coalesce(site_name_overrides[site_name], site_name),
      invalid_coordinates = latitude == 0 & longitude == 0
    )

  unresolved <- dat |> filter(!reef_label %in% manta$reef_label)
  if (nrow(unresolved) > 0) {
    warning(
      label, ": dropping ", nrow(unresolved), " detection(s) with reef_label(s) not in ",
      "the manta data (likely a trip predating manta's voyage coverage): ",
      paste(unique(unresolved$reef_label), collapse = ", "),
      call. = FALSE
    )
  }
  dat <- dat |> filter(reef_label %in% manta$reef_label)

  dat_pts <- dat |> st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  dat$cull_site <- ifelse(
    dat$invalid_coordinates, NA_character_,
    st_join(dat_pts, cull_sf["Name"], join = st_within)$Name
  )

  dat <- bind_cols(dat, assign_nearest_tow(
    dat$longitude, dat$latitude, dat$reef_label, manta_lines,
    max_dist_m = Inf, ref_date = dat$detection_date
  ))

  detections_out <- dat |>
    select(
      detection_id, trip_number, site_name, reef_label, cull_site,
      latitude, longitude, depth_m, detection_dttm,
      cots_video_sequence_length, ML_model_confidence, certified_by,
      nearest_tow_id = tow_id, dist_to_nearest_tow_m = dist_to_tow_m,
      nearest_tow_date = tow_date, date_diff_days, invalid_coordinates
    )

  valid <- dat |> filter(!invalid_coordinates)
  col_name <- paste0("cots_reefscan_", label)

  list(
    detections       = detections_out,
    tow_counts       = valid |> filter(!is.na(tow_id)) |> count(tow_id, name = col_name),
    cull_site_counts = valid |> filter(!is.na(cull_site)) |> count(cull_site, name = col_name),
    reef_counts      = valid |> count(reef_label, name = col_name)
  )
}

deep <- process_cotspotter_platform(cots_raw, "AIMS_RSD", "deep",
                                    site_name_overrides, manta, cull_sf, manta_lines)
surf <- process_cotspotter_platform(cots_raw, "AIMS_RSS", "surface",
                                    site_name_overrides, manta, cull_sf, manta_lines)

write_csv(deep$detections,       file.path(out_dir, "cotspotter_deep_detections.csv"))
write_csv(deep$tow_counts,       file.path(out_dir, "cotspotter_deep_tow_counts.csv"))
write_csv(deep$cull_site_counts, file.path(out_dir, "cotspotter_deep_cull_site_counts.csv"))
write_csv(deep$reef_counts,      file.path(out_dir, "cotspotter_deep_reef_counts.csv"))

write_csv(surf$detections,       file.path(out_dir, "cotspotter_surface_detections.csv"))
write_csv(surf$tow_counts,       file.path(out_dir, "cotspotter_surface_tow_counts.csv"))
write_csv(surf$cull_site_counts, file.path(out_dir, "cotspotter_surface_cull_site_counts.csv"))
write_csv(surf$reef_counts,      file.path(out_dir, "cotspotter_surface_reef_counts.csv"))

cat("Deep:    ", nrow(deep$detections), "detections |", nrow(deep$tow_counts), "tows |",
    nrow(deep$cull_site_counts), "cull sites |", nrow(deep$reef_counts), "reefs\n")
cat("Surface: ", nrow(surf$detections), "detections |", nrow(surf$tow_counts), "tows |",
    nrow(surf$cull_site_counts), "cull sites |", nrow(surf$reef_counts), "reefs\n")
