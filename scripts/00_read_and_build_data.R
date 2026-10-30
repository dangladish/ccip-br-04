## 00_read_and_build_data.R
##
## Reads every raw source and builds the three canonical analysis tables used
## by all downstream scripts:
##   tow_level        one row per manta tow segment
##   cull_site_level  one row per cull site
##   reef_level       one row per reef
## Sources: manta tows (EotR workbooks + KML tracks), ReefScan Surface/Deep points,
## eDNA, RHIS, cull-site polygons, and COTSpotter counts (built by
## 00a_build_cotspotter_detections.R, sourced below). Analytic standard errors are attached.
##
## Voyage exclusion: three FMPCRS voyages are dropped from the manta tows,
## ReefScan points and COTSpotter counts (set in `exclude_voyages` below):
##   2024-05 and 2026-03 (precautionary); 2025-01 (RS Surface GPS freeze,
##   i.e. many consecutive images with identical coordinates).
## eDNA and RHIS are not filtered by voyage.
##
## Run from the repository root.

library(tidyverse)
library(sf)
library(arrow)
library(readxl)
library(lubridate)
library(xml2)

dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)

# ---- Coral-cover midpoint scale (manta ordinal categories -> proportions) ----
{
coral_midpoints <- c(
  "0"  = 0.000,
  "1-" = 0.025, "1+" = 0.075,
  "2-" = 0.150, "2+" = 0.250,
  "3-" = 0.350, "3+" = 0.450,
  "4-" = 0.550, "4+" = 0.650,
  "5-" = 0.750, "5+" = 0.875
)
coral_levels <- names(coral_midpoints)
}

# ---- Uncertainty helpers: analytic standard errors ----
# Simple Wald / normal-approximation / Poisson SEs (not bootstrapped, since some
# reefs and cull sites have very few tows). Poisson SEs for COTS counts ignore
# overdispersion, so they are a lower bound.
#   wtd_mean_se()  tow-distance-weighted mean (Cochran 1977 approximation)
#   mean_se()      plain mean
#   binom_se()     proportion (n = classified points, or n = tows for rates)
#   pois_sum_se() / pois_rate_se()  count sum / rate
{
wtd_mean_se <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]; w <- w[ok]; n <- length(x)
  if (n < 2 || sum(w) == 0) return(NA_real_)
  m <- weighted.mean(x, w)
  sqrt(sum(w^2 * (x - m)^2) / sum(w)^2 * n / (n - 1))
}
mean_se <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}
binom_se <- function(p, n) {
  ifelse(is.na(p) | is.na(n) | n == 0, NA_real_, sqrt(p * (1 - p) / n))
}
pois_sum_se  <- function(count) ifelse(is.na(count), NA_real_, sqrt(count))
pois_rate_se <- function(count, n) ifelse(is.na(count) | is.na(n) | n == 0, NA_real_, sqrt(count) / n)
}

# ---- Ingest: original exports -> manta tows + ReefScan Surface points matched to tows ----
# Ported from the ingest scripts in github.com/open-AIMS/digital-manta-tow (same author).
# Inputs (data/raw/): manta_tows/QPWS_FMPCRS_*/ EotR report workbooks and HardCoral KML tracks
# (Whitsundays 2024-05: old-format tender KMLs); reefscan_surface_coral/*.csv[.gz], the ReefCloud
# point-classification export (only the columns read below are needed).
# Outputs (data/derived/ingest/): manta_tows.parquet (one row per tow, all voyages) and
# reefscan_surface_tows.parquet (one row per point, with matched tow_id).
# Image -> tow, dist_cutoff_m = 25: (1) candidate tows = same reef code; (2) tow whose UTC window
# (AEST start - 10 h, + distance/speed) contains the image time, nearest KML track if several;
# (3) a temporal match > 25 m from its track goes to the nearest same-reef tow, else dropped;
# (4) no temporal match: nearest same-reef tow within 25 m, else dropped.
# Cached; rebuilt if an output is missing or options(br04.rebuild_ingest = TRUE).
ingest_dir     <- "data/derived/ingest"
ingest_manta   <- file.path(ingest_dir, "manta_tows.parquet")
ingest_rss     <- file.path(ingest_dir, "reefscan_surface_tows.parquet")
rebuild_ingest <- isTRUE(getOption("br04.rebuild_ingest", FALSE)) || !all(file.exists(c(ingest_manta, ingest_rss)))
dir.create(ingest_dir, showWarnings = FALSE, recursive = TRUE)

if (rebuild_ingest) local({
  manta_dir     <- "data/raw/manta_tows"
  reefscan_dir  <- "data/raw/reefscan_surface_coral"
  dist_cutoff_m <- 25

  ## ---- A. manta tows: EotR xlsx exports -> one row per tow ----
  all_xlsx <- list.files(manta_dir, pattern = "Eotr_CotsSurveillanceReport.*\\.xlsx$", recursive = TRUE, full.names = TRUE)
  # tow-level files have no surveyor breakdown; observer-level files add Surveyor id/level
  tow_files <- Filter(\(x) !grepl("Obs.lev", x, ignore.case = TRUE), all_xlsx)
  obs_files <- Filter(\(x)  grepl("Obs.lev", x, ignore.case = TRUE), all_xlsx)
  stopifnot(length(tow_files) > 0)
  message("INGEST: reading ", length(tow_files), " tow-level and ", length(obs_files), " observer-level EotR workbooks")

  # all columns read as text (types differ across file versions); the id column's capitalisation varies
  manta_raw <- lapply(tow_files, \(f) {
    read_excel(f, col_types = "text") |>
      rename_with(\(nm) sub("(?i)cots surveillance id", "COTS Surveillance id", nm, perl = TRUE)) |>
      mutate(survey_folder = basename(dirname(f)), source_file = basename(f))
  }) |>
    bind_rows() |>
    filter(!is.na(`Tow time`))

  # `Tow time` is a "dd/mm/yyyy HH:MM AM/PM" string or an Excel serial number. The values are AEST
  # but labelled UTC here (kept: the matcher subtracts 10 h to get true UTC).
  manta_tows <- manta_raw |>
    mutate(
      tow_id = trimws(`COTS Surveillance id`),
      tow_datetime = case_when(
        grepl("^[0-9.]+$", `Tow time`) ~ as.POSIXct(as.numeric(`Tow time`) * 86400, origin = "1899-12-30", tz = "UTC"),
        TRUE ~ suppressWarnings(parse_date_time(`Tow time`, orders = c("dmy I:M p", "dmy HM"), tz = "UTC"))
      ),
      tow_date = as_date(tow_datetime),
      across(c(`Start lat`, `Start long`, `End lat`, `End long`, `Distance (m)`, `Average speed (km/h)`), as.numeric)
    ) |>
    rename(
      voyage = Voyage, vessel = Vessel,
      start_lat = `Start lat`, start_lon = `Start long`, end_lat = `End lat`, end_lon = `End long`,
      tow_distance_m = `Distance (m)`, avg_speed_kmh = `Average speed (km/h)`
    ) |>
    select(
      tow_id, tow_date, tow_datetime, voyage, vessel, survey_folder, source_file,
      start_lat, start_lon, end_lat, end_lon, tow_distance_m, avg_speed_kmh,
      everything(), -`COTS Surveillance id`, -`Tow time`
    )

  surveyor_info <- lapply(obs_files, \(f) {
    read_excel(f, col_types = "text") |>
      rename_with(\(nm) sub("(?i)cots surveillance id", "tow_id", nm, perl = TRUE)) |>
      filter(!is.na(`Tow time`)) |>
      mutate(tow_id = trimws(tow_id)) |>
      select(tow_id, surveyor_id = `Surveyor id`, surveyor_level = `Surveyor level`)
  }) |>
    bind_rows() |>
    distinct(tow_id, .keep_all = TRUE)
  manta_tows <- manta_tows |> left_join(surveyor_info, by = "tow_id")

  stopifnot(!anyNA(manta_tows$tow_id), !anyNA(manta_tows$tow_date))
  # Written and re-read before the matching: parquet
  # stores times at microsecond precision, and the image-in-tow-window test is sensitive to that
  # at window edges (17 of ~161,000 images switch temporal/spatial label otherwise).
  arrow::write_parquet(manta_tows, ingest_manta)
  manta_tows <- arrow::read_parquet(ingest_manta)
  message("INGEST: ", nrow(manta_tows), " tows; ", sum(manta_tows$tow_id %in% surveyor_info$tow_id), " with surveyor info")

  ## ---- B. KML track of each tow ----
  # New format (every voyage but Whitsundays 2024-05): tow id in the <name> tag of the
  # Eotr_CotsSurveillance_HardCoral* KML. Whitsundays 2024-05 (old format): two tender KMLs
  # whose placemarks carry only a "Survey time", matched to tows by time window.
  old_format_folders <- "QPWS_FMPCRS_2024-05_Whitsundays_20241015"

  kml_tow_ids <- function(kml_path) {
    kml <- read_xml(kml_path); ns <- xml_ns(kml)
    pms <- xml_find_all(kml, ".//d1:Placemark[.//d1:LineString]", ns)
    ids <- xml_text(xml_find_first(pms, "d1:name", ns))
    trimws(ids[grepl("^[0-9]+$", trimws(ids))])
  }
  kml_datetimes <- function(kml_path) {
    kml <- read_xml(kml_path); ns <- xml_ns(kml)
    pms <- xml_find_all(kml, ".//d1:Placemark[.//d1:LineString]", ns)
    times <- xml_text(xml_find_first(pms, ".//d1:Data[@name='Survey time']/d1:value", ns))
    suppressWarnings(parse_date_time(times, orders = c("dmy I:M p", "dmy HM")))
  }
  kml_centroids <- function(kml_path) {
    kml <- read_xml(kml_path); ns <- xml_ns(kml)
    pms <- xml_find_all(kml, ".//d1:Placemark[.//d1:LineString]", ns)
    lapply(pms, \(pm) {
      coord_txt <- xml_text(xml_find_first(pm, ".//d1:LineString/d1:coordinates", ns))
      pts <- strsplit(strsplit(trimws(coord_txt), "\\s+")[[1]], ",")
      c(cen_lon = mean(as.numeric(sapply(pts, `[`, 1)), na.rm = TRUE),
        cen_lat = mean(as.numeric(sapply(pts, `[`, 2)), na.rm = TRUE))
    }) |>
      do.call(rbind, args = _) |>
      as.data.frame()
  }

  tows <- manta_tows |>
    mutate(
      kml_format = if_else(survey_folder %in% old_format_folders, "old", "new"),
      kml_path   = NA_character_,
      kml_match  = NA_character_
    )

  for (folder in unique(tows$survey_folder[tows$kml_format == "new"])) {
    kmls <- list.files(file.path(manta_dir, folder), pattern = "HardCoral.*\\.kml$",
                       recursive = FALSE, full.names = TRUE, ignore.case = TRUE)
    if (length(kmls) == 0) { warning("No HardCoral KML found for: ", folder); next }
    missing_ids <- setdiff(tows$tow_id[tows$survey_folder == folder], kml_tow_ids(kmls[[1]]))
    if (length(missing_ids) > 0) warning(length(missing_ids), " tow IDs not found in KML for ", folder)
    kml_found <- kmls[[1]]
    tows <- tows |>
      mutate(
        kml_path  = if_else(survey_folder == folder, kml_found, kml_path),
        kml_match = if_else(survey_folder == folder, "tow_id",  kml_match)
      )
  }

  # Whitsundays 2024-05: KML survey_time is AEST like tow_datetime. A placemark belongs to a tow if
  # its time is inside [tow_datetime, tow_datetime + distance/speed]; if several do, the one whose
  # centroid is nearest the tow start wins.
  whits_folder <- old_format_folders
  whits_dir    <- file.path(manta_dir, whits_folder, "EotR_App_Exports")
  whits_kmls   <- list(
    tender1 = list.files(file.path(whits_dir, "Tender 1"), pattern = "Hard_Coral.*\\.kml$", full.names = TRUE, ignore.case = TRUE)[[1]],
    tender2 = list.files(file.path(whits_dir, "Tender 2"), pattern = "Hard_Coral.*\\.kml$", full.names = TRUE, ignore.case = TRUE)[[1]]
  )
  whits_rows       <- which(tows$survey_folder == whits_folder)
  whits_tow_starts <- manta_tows$tow_datetime[whits_rows]
  whits_speeds     <- manta_tows$avg_speed_kmh[whits_rows]
  # speed <= 0 would give an infinite end time; NA skips those tows
  whits_tow_ends <- whits_tow_starts +
    ifelse(!is.na(whits_speeds) & whits_speeds > 0,
           manta_tows$tow_distance_m[whits_rows] / (whits_speeds * 1000 / 3600), NA_real_)

  t1_dts  <- kml_datetimes(whits_kmls$tender1);  t2_dts  <- kml_datetimes(whits_kmls$tender2)
  t1_cens <- kml_centroids(whits_kmls$tender1);  t2_cens <- kml_centroids(whits_kmls$tender2)

  for (i in seq_along(whits_rows)) {
    row_i   <- whits_rows[i]
    t_start <- whits_tow_starts[i]
    t_end   <- whits_tow_ends[i]
    if (is.na(t_start) || is.na(t_end)) next
    t1_hit <- which(t1_dts >= t_start & t1_dts <= t_end)
    t2_hit <- which(t2_dts >= t_start & t2_dts <= t_end)
    n_hits <- length(t1_hit) + length(t2_hit)
    if (n_hits == 0) next
    if (n_hits == 1) {
      tows$kml_path[row_i]  <- if (length(t1_hit) == 1) whits_kmls$tender1 else whits_kmls$tender2
      tows$kml_match[row_i] <- "datetime"
      next
    }
    candidates <- bind_rows(
      if (length(t1_hit) > 0) data.frame(kml = whits_kmls$tender1, cen_lon = t1_cens$cen_lon[t1_hit], cen_lat = t1_cens$cen_lat[t1_hit]),
      if (length(t2_hit) > 0) data.frame(kml = whits_kmls$tender2, cen_lon = t2_cens$cen_lon[t2_hit], cen_lat = t2_cens$cen_lat[t2_hit])
    ) |>
      mutate(dist = sqrt((cen_lon - tows$start_lon[row_i])^2 + (cen_lat - tows$start_lat[row_i])^2))
    best <- candidates[which.min(candidates$dist), ]
    tows$kml_path[row_i]  <- best$kml
    tows$kml_match[row_i] <- "datetime_proximity"
  }
  message("INGEST: tows assigned a KML: ", sum(!is.na(tows$kml_path)), " / ", nrow(tows))

  ## ---- C. ReefScan Surface points (ReefCloud export) ----
  # Only the columns used downstream are read; the matcher is image-level.
  rs_csv <- list.files(reefscan_dir, pattern = "\\.csv(\\.gz)?$", full.names = TRUE)
  stopifnot("expected exactly one ReefScan Surface CSV" = length(rs_csv) == 1)
  rs_cols <- c("image_id", "image_name", "site_reef_name", "site_code", "survey_start_date (UTC)",
               "image_latitude", "image_longitude", "point_machine_group_code", "point_human_group_code")
  message("INGEST: reading ", basename(rs_csv))
  reefscan <- arrow::read_csv_arrow(
    rs_csv, col_select = dplyr::all_of(rs_cols),
    col_types = arrow::schema(image_id = arrow::float64(), image_latitude = arrow::float64(),
                              image_longitude = arrow::float64(), `survey_start_date (UTC)` = arrow::string()),
    parse_options = arrow::CsvParseOptions$create(newlines_in_values = TRUE)
  )
  # the export column is UTC ("YYYY-MM-DD HH:MM:SS+00:00"), so the date is its first 10 characters
  stopifnot(all(grepl("\\+00:00$", reefscan$`survey_start_date (UTC)`[!is.na(reefscan$`survey_start_date (UTC)`)][1:1000])))
  reefscan <- reefscan |>
    mutate(survey_date = as.Date(substr(`survey_start_date (UTC)`, 1, 10)), `survey_start_date (UTC)` = NULL)
  message("INGEST: ", nrow(reefscan), " ReefScan Surface points, ", n_distinct(reefscan$image_id), " images; ",
          sum(is.na(reefscan$survey_date)), " with no survey date")

  ## ---- D. match ReefScan images to manta tows ----
  # site_codes with no ##-### reef code
  reef_code_overrides <- c("Hardy Reef" = "19-135")
  extract_reef_code <- function(site_code) {
    overridden <- reef_code_overrides[site_code]
    m          <- regexpr("\\d{2}-\\d{3,4}[A-Za-z]?", site_code)
    extracted  <- ifelse(m != -1, tolower(substr(site_code, m, m + attr(m, "match.length") - 1)), NA_character_)
    ifelse(!is.na(overridden), overridden, extracted)
  }

  parse_kml_new <- function(kml_path) {
    kml <- read_xml(kml_path); ns <- xml_ns(kml)
    pms <- xml_find_all(kml, ".//d1:Placemark[.//d1:LineString]", ns)
    results <- lapply(pms, \(pm) {
      tow_id <- xml_text(xml_find_first(pm, "d1:name", ns))
      if (is.na(tow_id) || !grepl("^[0-9]+$", trimws(tow_id))) return(NULL)
      coord_txt <- xml_text(xml_find_first(pm, ".//d1:LineString/d1:coordinates", ns))
      pts <- strsplit(strsplit(trimws(coord_txt), "\\s+")[[1]], ",")
      lon <- as.numeric(sapply(pts, `[`, 1)); lat <- as.numeric(sapply(pts, `[`, 2))
      n <- length(lon)
      if (n < 2) return(NULL)
      data.frame(tow_id = trimws(tow_id), x1 = lon[-n], y1 = lat[-n], x2 = lon[-1], y2 = lat[-1])
    })
    bind_rows(results[!sapply(results, is.null)])
  }
  parse_kml_old <- function(kml_path) {
    kml <- read_xml(kml_path); ns <- xml_ns(kml)
    pms <- xml_find_all(kml, ".//d1:Placemark[.//d1:LineString]", ns)
    results <- lapply(pms, \(pm) {
      st <- xml_text(xml_find_first(pm, ".//d1:Data[@name='Survey time']/d1:value", ns))
      coord_txt <- xml_text(xml_find_first(pm, ".//d1:LineString/d1:coordinates", ns))
      pts <- strsplit(strsplit(trimws(coord_txt), "\\s+")[[1]], ",")
      lon <- as.numeric(sapply(pts, `[`, 1)); lat <- as.numeric(sapply(pts, `[`, 2))
      n <- length(lon)
      if (n < 2 || is.na(st)) return(NULL)
      data.frame(
        # UTC label to match tow_datetime's storage convention (AEST values labelled UTC)
        dt = suppressWarnings(parse_date_time(st, orders = c("dmy I:M p", "dmy HM"), tz = "UTC")),
        x1 = lon[-n], y1 = lat[-n], x2 = lon[-1], y2 = lat[-1]
      )
    })
    bind_rows(results[!sapply(results, is.null)])
  }
  # centroid-to-point distance (m): tiebreak when several KML tracks fall in one tow window
  track_centroid_dist_m <- function(x1, y1, x2, y2, target_lon, target_lat) {
    cx <- mean(c(x1, x2)); cy <- mean(c(y1, y2))
    cos_lat <- cos(target_lat * pi / 180)
    sqrt(((cx - target_lon) * 111320 * cos_lat)^2 + ((cy - target_lat) * 110540)^2)
  }

  tows_matched <- tows |> filter(!is.na(kml_path), !is.na(kml_match))

  new_segs <- tows_matched |>
    filter(kml_format == "new") |>
    pull(kml_path) |>
    unique() |>
    lapply(\(p) parse_kml_new(p) |> filter(tow_id %in% tows_matched$tow_id)) |>
    bind_rows()

  old_tows <- tows_matched |>
    filter(kml_format == "old", !is.na(avg_speed_kmh), avg_speed_kmh > 0) |>
    mutate(
      tow_duration_s = tow_distance_m / (avg_speed_kmh * 1000 / 3600),
      tow_end_aest   = tow_datetime + tow_duration_s
    ) |>
    select(tow_id, tow_datetime, tow_end_aest, start_lat, start_lon, kml_path)

  old_segs <- old_tows |>
    pull(kml_path) |>
    unique() |>
    lapply(\(p) {
      raw <- parse_kml_old(p)
      if (is.null(raw) || nrow(raw) == 0) return(NULL)
      tows_this <- old_tows |> filter(kml_path == p)
      lapply(seq_len(nrow(tows_this)), \(i) {
        t_start <- tows_this$tow_datetime[i]
        t_end   <- tows_this$tow_end_aest[i]
        if (is.na(t_start) || is.na(t_end)) return(NULL)
        in_window <- which(raw$dt >= t_start & raw$dt <= t_end)
        if (length(in_window) == 0) return(NULL)
        unique_dts <- unique(raw$dt[in_window])
        if (length(unique_dts) > 1) {
          centroid_dists <- vapply(unique_dts, \(dt_) {
            idx <- which(raw$dt == dt_)[1]
            track_centroid_dist_m(raw$x1[idx], raw$y1[idx], raw$x2[idx], raw$y2[idx],
                                  tows_this$start_lon[i], tows_this$start_lat[i])
          }, numeric(1))
          in_window <- which(raw$dt == unique_dts[which.min(centroid_dists)])
        }
        raw[in_window, ] |>
          mutate(tow_id = tows_this$tow_id[i]) |>
          select(tow_id, x1, y1, x2, y2)
      }) |>
        bind_rows()
    }) |>
    bind_rows()

  kml_segments <- bind_rows(new_segs, old_segs) |> distinct(tow_id, x1, y1, x2, y2)
  message(sprintf("INGEST: KML segments: %d across %d tows", nrow(kml_segments), n_distinct(kml_segments$tow_id)))

  # Planar metre approximation at the mean latitude of all images (reef scale, tens of km)
  m_per_deg_lat <- 111000
  m_per_deg_lon <- 111000 * cos(mean(range(reefscan$image_latitude, na.rm = TRUE)) * pi / 180)

  kml_seg_sc <- kml_segments |>
    mutate(sx1 = x1 * m_per_deg_lon, sy1 = y1 * m_per_deg_lat, sx2 = x2 * m_per_deg_lon, sy2 = y2 * m_per_deg_lat) |>
    group_by(tow_id) |>
    mutate(
      seg_len_m   = sqrt((sx2 - sx1)^2 + (sy2 - sy1)^2),
      # distance along the tow from its start to the beginning of each segment
      seg_start_m = cumsum(lag(seg_len_m, default = 0))
    ) |>
    ungroup()

  # image time (UTC) is embedded in image_name as YYYYMMDD_HHMMSS
  reefscan <- reefscan |>
    mutate(
      image_datetime_utc = as.POSIXct(
        ifelse(grepl("\\d{8}_\\d{6}", image_name), sub(".*(\\d{8}_\\d{6}).*", "\\1", image_name), NA_character_),
        format = "%Y%m%d_%H%M%S", tz = "UTC"
      )
    )

  # tow_datetime is AEST (UTC+10) with a UTC label: subtract 10 h; end time from speed and distance
  tow_windows <- tows |>
    filter(!is.na(tow_datetime), !is.na(avg_speed_kmh), avg_speed_kmh > 0, !is.na(tow_distance_m)) |>
    mutate(
      tow_duration_s = tow_distance_m / (avg_speed_kmh * 1000 / 3600),
      tow_start_utc  = tow_datetime - 10 * 3600,
      tow_end_utc    = tow_start_utc + tow_duration_s,
      reef_code      = tolower(`Reef label`)
    ) |>
    select(tow_id, reef_code, tow_start_utc, tow_end_utc)

  # all points of an image share one GPS fix: match once per image, join back to the points
  rs_images <- reefscan |>
    distinct(image_id, .keep_all = TRUE) |>
    filter(!is.na(image_latitude), !is.na(image_longitude)) |>
    mutate(reef_code = extract_reef_code(site_code), rx = image_longitude * m_per_deg_lon, ry = image_latitude * m_per_deg_lat)
  message(sprintf("INGEST: matching %d images (%d with a reef code)", nrow(rs_images), sum(!is.na(rs_images$reef_code))))

  # per (image, candidate tow): minimum distance to any segment of the tow's KML track, and the
  # along-track distance to the foot of that perpendicular
  kml_min_dist <- function(img_tow_df) {
    img_tow_df |>
      left_join(kml_seg_sc, by = "tow_id", relationship = "many-to-many") |>
      mutate(
        dx = sx2 - sx1, dy = sy2 - sy1, len_sq = dx^2 + dy^2,
        t = ifelse(len_sq == 0, 0, pmax(0, pmin(1, ((rx - sx1) * dx + (ry - sy1) * dy) / len_sq))),
        dist_m = sqrt((rx - (sx1 + t * dx))^2 + (ry - (sy1 + t * dy))^2),
        cumul_at_foot = seg_start_m + t * seg_len_m
      ) |>
      group_by(image_id, tow_id) |>
      slice_min(dist_m, n = 1, with_ties = FALSE) |>
      ungroup() |>
      select(image_id, tow_id, dist_to_tow_m = dist_m, cumul_dist_tow = cumul_at_foot)
  }

  # Stage 1: every image x every tow at the same reef
  candidates <- rs_images |>
    select(image_id, reef_code, image_datetime_utc, rx, ry) |>
    inner_join(tow_windows, by = "reef_code", relationship = "many-to-many") |>
    mutate(in_window = image_datetime_utc >= tow_start_utc & image_datetime_utc <= tow_end_utc)

  temporal_ids <- candidates |> filter(in_window) |> pull(image_id) |> unique()
  single_t_ids <- candidates |> filter(in_window) |> count(image_id) |> filter(n == 1) |> pull(image_id)
  multi_t_ids  <- candidates |> filter(in_window) |> count(image_id) |> filter(n > 1)  |> pull(image_id)
  no_temp_ids  <- setdiff(candidates |> pull(image_id) |> unique(), temporal_ids)

  # Stage 2: temporal assignment (several windows -> KML-nearest tow)
  assign_temporal <- candidates |>
    filter(image_id %in% single_t_ids, in_window) |>
    select(image_id, tow_id, rx, ry) |>
    kml_min_dist() |>
    mutate(match_method = "temporal")
  assign_temporal <- bind_rows(
    assign_temporal,
    candidates |>
      filter(image_id %in% multi_t_ids, in_window) |>
      select(image_id, tow_id, rx, ry) |>
      kml_min_dist() |>
      group_by(image_id) |>
      slice_min(dist_to_tow_m, n = 1, with_ties = FALSE) |>
      ungroup() |>
      mutate(match_method = "temporal")
  )

  # Stage 3: temporal matches beyond the cutoff -> nearest same-reef tow, dropped if still beyond
  temporal_ok     <- assign_temporal |> filter(dist_to_tow_m <= dist_cutoff_m)
  need_rescue_ids <- assign_temporal |> filter(dist_to_tow_m > dist_cutoff_m) |> pull(image_id)
  assign_rescued <- candidates |>
    filter(image_id %in% need_rescue_ids) |>
    distinct(image_id, tow_id, rx, ry) |>
    kml_min_dist() |>
    group_by(image_id) |>
    slice_min(dist_to_tow_m, n = 1, with_ties = FALSE) |>
    ungroup() |>
    filter(dist_to_tow_m <= dist_cutoff_m) |>
    mutate(match_method = "spatial")

  # Stage 4: no temporal match -> nearest same-reef tow within the cutoff
  assign_spatial <- candidates |>
    filter(image_id %in% no_temp_ids) |>
    distinct(image_id, tow_id, rx, ry) |>
    kml_min_dist() |>
    group_by(image_id) |>
    slice_min(dist_to_tow_m, n = 1, with_ties = FALSE) |>
    ungroup() |>
    filter(dist_to_tow_m <= dist_cutoff_m) |>
    mutate(match_method = "spatial")

  assignments <- bind_rows(temporal_ok, assign_rescued, assign_spatial) |> distinct(image_id, .keep_all = TRUE)
  message(sprintf("INGEST: matched %d / %d images (%d temporal, %d spatial); %d without a reef/tow candidate",
                  nrow(assignments), nrow(rs_images), sum(assignments$match_method == "temporal"),
                  sum(assignments$match_method == "spatial"), nrow(rs_images) - n_distinct(candidates$image_id)))

  reefscan_tows <- reefscan |>
    select(-image_datetime_utc, -image_name) |>
    left_join(assignments |> select(image_id, tow_id, dist_to_tow_m, cumul_dist_tow, match_method), by = "image_id")

  # (manta_tows.parquet was already written, before the matching)
  arrow::write_parquet(reefscan_tows, ingest_rss)
  message("INGEST: written ", ingest_manta, " and ", ingest_rss)
})
if (!rebuild_ingest) message("INGEST: using cached ", ingest_dir, " (options(br04.rebuild_ingest = TRUE) to rebuild)")

# ---- Read manta tow data ----
{
mt <- arrow::read_parquet(ingest_manta) |>
  rename(
    reef_label    = `Reef label`,
    reef_name     = `Reef name`,
    cull_site     = `Cull site`,
    dist_cull_m   = `Distance to cull site (m)`,
    cots_observed = `COTS observed`,
    feeding_scars = `Feeding scars`,
    hard_coral    = `Hard coral`,
    soft_coral    = `Soft coral`,
    recently_dead = `Recently dead coral`,
    bleach_sev    = `Bleaching severity`,
    bleach_ext    = `Bleaching extent (%)`
  ) |>
  mutate(
    cots_observed  = as.integer(cots_observed),
    tow_distance_m = as.numeric(tow_distance_m),
    avg_speed_kmh  = as.numeric(avg_speed_kmh),

    hard_coral    = factor(hard_coral, levels = coral_levels, ordered = TRUE),
    soft_coral    = factor(soft_coral, levels = coral_levels, ordered = TRUE),
    recently_dead = factor(recently_dead, levels = coral_levels, ordered = TRUE),

    perc_hard = coral_midpoints[as.character(hard_coral)],
    perc_soft = coral_midpoints[as.character(soft_coral)],

    reef_id = factor(reef_label),
    voyage  = factor(voyage),
    vessel  = factor(vessel),

    # Presence indicator; the parquet holds only the a/p/c scar category
    # (absent / present = 1-10 / common = >10 scars per tow).
    scar_present = as.integer(feeding_scars %in% c("p", "c"))
  ) |>
  filter(!is.na(tow_distance_m), tow_distance_m > 0)

# ---- Feeding-scar count per tow ----
# ScarsCount is recorded for tows from 14 June 2025 (FMPCRS 2025-03 onward); extracted from the
# workbook by scripts/jobs/extract_manta_scar_counts.pl. FMPCRS 2025-02 has no count, so its
# scarred tows are imputed from the category: present = 2, common = 11.
scar_impute <- c(p = 2, c = 11)
scar_counts <- readr::read_csv("data/raw/manta_tows/manta_scar_counts.csv", col_types = "ccd", progress = FALSE) |>
  transmute(tow_id, scar_code_wb = scar_code, scars_count_obs = scars_count)
mt <- mt |>
  left_join(scar_counts, by = "tow_id") |>
  mutate(
    scar_count_imputed = is.na(scars_count_obs) & scar_present == 1,
    scar_count = case_when(
      !is.na(scars_count_obs) ~ scars_count_obs,
      scar_present == 0       ~ 0,
      TRUE                    ~ unname(scar_impute[feeding_scars])
    )
  )
stopifnot(
  all(mt$scar_code_wb == mt$feeding_scars, na.rm = TRUE),
  all(mt$scars_count_obs[mt$feeding_scars == "a"] == 0, na.rm = TRUE),
  all(mt$scars_count_obs[mt$feeding_scars == "p"] %in% c(NA, 1:10)),
  all(mt$scars_count_obs[mt$feeding_scars == "c"] >= 11, na.rm = TRUE),
  !anyNA(mt$scar_count)
)
cat("scar count: recorded for", sum(!is.na(mt$scars_count_obs)), "of", nrow(mt), "tows; imputed for",
    sum(mt$scar_count_imputed), "scarred tows (present = 2, common = 11)
")

cat("manta (pre-exclusion):", nrow(mt), "tows,", n_distinct(mt$reef_label), "reefs\n")

## ---- Voyage exclusion (see header) ----
## mt_allvoyages, voyage_windows and exclude_voyages are reused below to filter
## ReefScan points and COTSpotter counts.
exclude_voyages <- c("FMPCRS 2024-05", "FMPCRS 2025-01", "FMPCRS 2026-03")
mt_allvoyages   <- mt
voyage_windows  <- mt_allvoyages |> summarise(start = min(tow_date), end = max(tow_date), .by = voyage)
mt <- mt |> filter(!voyage %in% exclude_voyages)

cat("manta (post 3-voyage exclusion):", nrow(mt), "tows,", n_distinct(mt$reef_label), "reefs",
    "-- excluded", sum(mt_allvoyages$voyage %in% exclude_voyages), "tows from", paste(exclude_voyages, collapse = ", "), "\n")
}

# ---- Read ReefScan Surface / Deep point-level data ----
{
rs_surface_raw <- arrow::read_parquet(ingest_rss) |>
  select(
    tow_id, site_reef_name, site_code,
    survey_date, image_latitude, image_longitude,
    point_machine_group_code, point_human_group_code
  )

# No tow_id in this source -- known limitation; RS Deep is matched to tows spatially below.
rs_deep_raw <- read_csv(
  "data/raw/reefscan_deep/reefcloud-point-classification-GRB-ReefScan-Deep-2026-08-13-filtered.csv",
  show_col_types = FALSE, progress = FALSE
) |>
  select(
    site_reef_name, site_code, `survey_start_date (UTC)`,
    image_latitude, image_longitude,
    point_machine_group_code, point_human_group_code
  )

cat("RS Surface points (pre-exclusion):", nrow(rs_surface_raw), " RS Deep points (pre-exclusion):", nrow(rs_deep_raw), "\n")

## ---- Drop excluded voyages from ReefScan points ----
## Cull-site/reef aggregation of ReefScan points is independent of manta tows,
## so the imagery itself must be filtered. RS Surface: voyage via its tow_id in
## the full manta table, falling back to the survey-date window of each voyage.
## RS Deep has no tow_id, so only the date window is used.
in_excluded_window <- function(d) {
  Reduce(`|`, lapply(exclude_voyages, function(v) {
    w <- voyage_windows[voyage_windows$voyage == v, ]
    !is.na(d) & d >= w$start & d <= w$end
  }))
}

rs_surf_voyage <- rs_surface_raw |>
  left_join(mt_allvoyages |> distinct(tow_id, voyage), by = "tow_id") |>
  pull(voyage)
rs_surf_exclude <- (!is.na(rs_surf_voyage) & rs_surf_voyage %in% exclude_voyages) |
  (is.na(rs_surf_voyage) & in_excluded_window(rs_surface_raw$survey_date))
rs_deep_exclude <- in_excluded_window(as.Date(rs_deep_raw$`survey_start_date (UTC)`))

n_surf_pre <- nrow(rs_surface_raw); n_deep_pre <- nrow(rs_deep_raw)
rs_surface_raw <- rs_surface_raw[!rs_surf_exclude, ]
rs_deep_raw    <- rs_deep_raw[!rs_deep_exclude, ]
rm(rs_surf_voyage, rs_surf_exclude, rs_deep_exclude, in_excluded_window, mt_allvoyages, voyage_windows)

cat("RS Surface points (post 3-voyage exclusion):", nrow(rs_surface_raw), "of", n_surf_pre,
    "-- excluded", n_surf_pre - nrow(rs_surface_raw),
    paste0("(", round(100 * (n_surf_pre - nrow(rs_surface_raw)) / n_surf_pre, 1), "%)"), "\n")
cat("RS Deep points (post 3-voyage exclusion):", nrow(rs_deep_raw), "of", n_deep_pre,
    "-- excluded", n_deep_pre - nrow(rs_deep_raw),
    paste0("(", round(100 * (n_deep_pre - nrow(rs_deep_raw)) / n_deep_pre, 1), "%)"), "\n")
rm(n_surf_pre, n_deep_pre)
}

# ---- Read eDNA data (reef-level summary + raw samples) ----
{
edna_reef_raw <- read_csv("data/raw/eDNA/RJFMP25_26 Results.csv", show_col_types = FALSE)

edna_raw <- read_csv("data/raw/eDNA/eDNA_BR_02.csv", show_col_types = FALSE) |>
  rename(
    sample_id           = `...1`,
    date                = Date,
    year                = Year,
    edna_reef           = Reef,
    edna_reef_name      = ReefName,
    site_name           = Site_name,
    site_no             = Site_No,
    lat                 = Lat,
    long                = Long,
    voyage              = Voyage,
    project             = Project,
    collection_org      = Collection.organisation,
    collection_vessel   = Collection.vessel,
    reef_rep            = Reef_rep,
    site_rep            = Site_rep,
    pcr1_binary         = PCR1_binary,
    pcr2_binary         = PCR2_binary,
    lod_sample_positive = LOD_sample_positive,
    pcr1_conc           = PCR1_conc,
    pcr2_conc           = PCR2_conc,
    conc_mean           = Conc_mean
  ) |>
  mutate(
    date      = lubridate::dmy(date),
    voyage    = factor(voyage),
    project   = factor(project),
    edna_reef = as.character(edna_reef)
  )

cat("eDNA reef summary:", nrow(edna_reef_raw), "reefs. eDNA raw samples:", nrow(edna_raw), "\n")
}

# ---- Read RHIS data ----
# One row per survey; ReefLabel matches manta's reef_label directly.
{
rhis_raw <- read_excel(
  "data/raw/manta_tows/260529-COTS-Manta-Cull-RHIS.xlsx",
  sheet = "RHIS"
)

# Blank predator/scar counts are zero (workbook README), so they are zero-filled and every
# survey counts toward COTS means and rates; has_predator_survey records non-blank entries.
# Hard coral cover = live coral (% of benthos) x summed live share of the 7 hard growth forms;
# soft = live x soft share. Cover is 0 where live coral is 0; surveys without growth-form data
# (79; 18 with live coral > 0) stay NA.
rhis_hard_cols <- paste0("Coral", c("Branching", "Bushy", "PlateTable", "VaseFoliose",
                                    "Encrusting", "Mushroom", "Massive"), "LiveProportion")
rhis_scar_cols <- paste0("Predator", c("Soft", "Branching", "Bushy", "PlateTable", "VaseFoliose",
                                       "Encrusting", "Mushroom", "Massive"), "CrownOfThornsStarfishScarCount")
# RHIS is restricted to the two seasons the voyages cover (2024/25 and 2025/26, July-June).
rhis_window <- as.Date(c("2024-07-01", "2026-06-30"))
rhis <- rhis_raw |>
  mutate(
    hard_live_share = rowSums(across(all_of(rhis_hard_cols))),
    scar_total      = rowSums(across(all_of(rhis_scar_cols)), na.rm = TRUE)
  ) |>
  transmute(
    reef_label = ReefLabel,
    cull_site  = CullSite,
    survey_id  = SurveyId,
    survey_date = as.Date(SurveyTime),
    latitude   = Latitude,
    longitude  = Longitude,
    hc_prop_survey = if_else(BenthosLiveCoralPercentage == 0, 0,
      BenthosLiveCoralPercentage / 100 * hard_live_share / 100),
    sc_prop_survey = if_else(BenthosLiveCoralPercentage == 0, 0,
      BenthosLiveCoralPercentage / 100 * CoralSoftLiveProportion / 100),
    has_predator_survey = !is.na(PredatorCrownOfThornsStarfishAdultIndividualCount),
    cots_count_survey = coalesce(PredatorCrownOfThornsStarfishAdultIndividualCount, 0) +
      coalesce(PredatorCrownOfThornsStarfishJuvenileIndividualCount, 0),
    scar_count_survey = scar_total
  ) |>
  filter(survey_date >= rhis_window[1], survey_date <= rhis_window[2])

cat("RHIS:", nrow(rhis), "surveys,", n_distinct(rhis$reef_label), "reefs,",
    sum(rhis$has_predator_survey), "with predator counts\n")
}

# ---- Read cull-site polygons (duplicate-Name polygons dissolved) ----
{
name_re <- "^(.*)_([0-9]{2}-[0-9]{3}[a-z]?)_([^_]+)$"

cull_raw <- sf::st_read(
  "data/raw/COTS Control Program Cull Sites - External/Cull_Sites_Master.shp", quiet = TRUE
) |>
  sf::st_zm()

cull_raw <- cull_raw |>
  mutate(
    .m         = str_match(Name, name_re),
    cs_prefix  = .m[, 2],
    cs_reef    = .m[, 3],
    cs_sitenum = .m[, 4]
  ) |>
  select(-.m)

cull_polygons <- cull_raw |>
  group_by(Name) |>
  summarise(
    GlobalID   = first(GlobalID),
    cs_prefix  = first(cs_prefix),
    cs_reef    = first(cs_reef),
    cs_sitenum = first(cs_sitenum),
    n_polys    = n(),
    .groups    = "drop"
  ) |>
  sf::st_make_valid()

stopifnot(nrow(cull_polygons) == n_distinct(cull_raw$Name))

cull_polygons <- cull_polygons |>
  rename(cull_site = Name) |>
  sf::st_transform(4326) |>
  sf::st_make_valid()

rm(cull_raw, name_re)
cat("Cull-site polygons:", nrow(cull_polygons), "distinct sites\n")
}

# ---- Reef ID crosswalk (manta reef_label <-> RS Surface / RS Deep / eDNA) ----
{
mt_canon <- mt |>
  distinct(reef_label, reef_name) |>
  mutate(reef_name_clean = str_remove(reef_name, "\\s*\\([^)]+\\)\\s*$") |> str_trim()) |>
  arrange(reef_label)

strip_reef_suffix <- function(x) str_remove(x, "\\s+Reef$") |> str_to_lower() |> str_trim()

# --- RS Surface: try site_code, then site_reef_name, then reef_name_clean ---
rs_surf_xwalk <- rs_surface_raw |>
  distinct(site_reef_name, site_code) |>
  mutate(key_code = strip_reef_suffix(site_code), key_name = strip_reef_suffix(site_reef_name)) |>
  left_join(
    mt_canon |> transmute(reef_label_by_code = reef_label, key_code = str_to_lower(reef_label)),
    by = "key_code"
  ) |>
  left_join(
    mt_canon |> transmute(reef_label_by_name = reef_label, key_name = str_to_lower(reef_label)),
    by = "key_name"
  ) |>
  left_join(
    mt_canon |> transmute(reef_label_by_clean_name = reef_label, key_name = strip_reef_suffix(reef_name_clean)),
    by = "key_name"
  ) |>
  mutate(reef_label = coalesce(reef_label_by_code, reef_label_by_name, reef_label_by_clean_name)) |>
  select(rs_surf_site_reef_name = site_reef_name, reef_label)

# --- RS Deep: fix the "20-263" -> "21-263" typo, strip " RSD [AB]"/" Reef" ---
rs_deep_xwalk <- rs_deep_raw |>
  distinct(site_reef_name, site_code) |>
  mutate(
    site_reef_name_fixed = if_else(
      str_starts(site_reef_name, "20-263"),
      str_replace(site_reef_name, "^20-263", "21-263"),
      site_reef_name
    ),
    join_key = str_remove(site_reef_name_fixed, "\\s+RSD\\s*[AB]?\\s*$") |>
      str_remove("\\s+Reef$") |> str_to_lower() |> str_trim()
  ) |>
  left_join(mt_canon |> transmute(reef_label, join_key = str_to_lower(reef_label)), by = "join_key") |>
  select(rs_deep_site_reef_name = site_reef_name, reef_label)

# --- eDNA: strip trailing " Reef"/plural "s", fall back to raw name match ---
edna_step1 <- edna_reef_raw |>
  distinct(Reef) |>
  mutate(edna_candidate = str_remove(Reef, "\\s+Reef$") |> str_remove("s$") |> str_trim()) |>
  left_join(mt_canon |> transmute(reef_label, reef_name_clean, join_key = reef_label), by = c("edna_candidate" = "join_key"))

edna_xwalk <- edna_step1 |>
  left_join(mt_canon |> transmute(reef_label2 = reef_label, join_key2 = reef_name_clean), by = c("Reef" = "join_key2")) |>
  mutate(reef_label = coalesce(reef_label, reef_label2)) |>
  select(edna_reef = Reef, reef_label)

reef_xwalk <- mt_canon |>
  select(reef_label, reef_name_clean) |>
  left_join(
    rs_surf_xwalk |> filter(!is.na(reef_label)) |> group_by(reef_label) |>
      summarise(rs_surf_site = paste(rs_surf_site_reef_name, collapse = "; "), .groups = "drop"),
    by = "reef_label"
  ) |>
  left_join(
    rs_deep_xwalk |> filter(!is.na(reef_label)) |> group_by(reef_label) |>
      summarise(rs_deep_sites = paste(rs_deep_site_reef_name, collapse = "; "), .groups = "drop"),
    by = "reef_label"
  ) |>
  left_join(edna_xwalk |> filter(!is.na(reef_label)) |> select(edna_reef, reef_label), by = "reef_label") |>
  mutate(
    has_manta   = TRUE,
    has_rs_surf = !is.na(rs_surf_site),
    has_rs_deep = !is.na(rs_deep_sites),
    has_edna    = !is.na(edna_reef),
    has_rhis    = reef_label %in% unique(rhis$reef_label)
  )

# Saved: the eDNA data prep needs the row-level RS Surface site_reef_name ->
# reef_label mapping (used with the 280 m eDNA-site-assignment cutoff).
saveRDS(rs_surf_xwalk |> filter(!is.na(reef_label)), "data/derived/rs_surf_site_xwalk.RDS")

rm(mt_canon, rs_surf_xwalk, rs_deep_xwalk, edna_step1, edna_xwalk, strip_reef_suffix)
cat("Reef crosswalk:", nrow(reef_xwalk), "reefs -- has_rs_surf:", sum(reef_xwalk$has_rs_surf),
    " has_rs_deep:", sum(reef_xwalk$has_rs_deep), " has_edna:", sum(reef_xwalk$has_edna),
    " has_rhis:", sum(reef_xwalk$has_rhis), "\n")
}

# ---- Manta tow linestrings (for nearest-tow matching) ----
{
manta_lines <- mt |>
  filter(!is.na(start_lat), !is.na(end_lat)) |>
  rowwise() |>
  mutate(geometry = st_sfc(st_linestring(matrix(c(start_lon, end_lon, start_lat, end_lat), ncol = 2)))) |>
  ungroup() |>
  st_as_sf(crs = 4326) |>
  st_transform(3577) |>
  select(tow_id, reef_label, tow_date, geometry)

cat("Manta tow lines built:", nrow(manta_lines), "\n")
}

# ---- Point-in-polygon cull-site assignment ----
assign_cull_site_inline <- function(lon, lat, cull_polygons) {
  n  <- length(lon)
  ok <- !is.na(lon) & !is.na(lat) & !(lon == 0 & lat == 0)
  out <- tibble(cull_site = rep(NA_character_, n), ambiguous_flag = rep(NA, n))
  if (!any(ok)) return(out)
  pts <- st_as_sf(tibble(.row = which(ok), lon = lon[ok], lat = lat[ok]), coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  hits <- st_within(pts, cull_polygons)
  n_sites <- vapply(hits, function(idx) if (length(idx) == 0) 0L else n_distinct(cull_polygons$cull_site[idx]), integer(1))
  idx1  <- vapply(hits, function(h) if (length(h) >= 1) h[[1]] else NA_integer_, integer(1))
  site1 <- cull_polygons$cull_site[idx1]
  out$cull_site[ok]      <- if_else(n_sites == 1L, site1, NA_character_)
  out$ambiguous_flag[ok] <- n_sites > 1L
  out
}

# ---- Nearest manta tow (used for RS Deep tow-level matching) ----
assign_nearest_tow_inline <- function(lon, lat, reef_label, manta_lines, max_dist_m = Inf) {
  n <- length(lon)
  out <- tibble(tow_id = rep(NA_character_, n), dist_to_tow_m = rep(NA_real_, n))
  ok <- !is.na(lon) & !is.na(lat) & !(lon == 0 & lat == 0) & !is.na(reef_label)
  if (!any(ok)) return(out)
  idx_ok <- which(ok)
  pts <- st_as_sf(
    tibble(.row = idx_ok, lon = lon[idx_ok], lat = lat[idx_ok], reef_label = reef_label[idx_ok]),
    coords = c("lon", "lat"), crs = 4326, remove = FALSE
  ) |> st_transform(3577)
  for (rl in unique(pts$reef_label)) {
    cand <- manta_lines[manta_lines$reef_label == rl, ]
    if (nrow(cand) == 0) next
    grp  <- pts[pts$reef_label == rl, ]
    dmat <- matrix(as.numeric(st_distance(grp, cand)), nrow = nrow(grp))
    j    <- apply(dmat, 1, which.min)
    dmin <- dmat[cbind(seq_len(nrow(grp)), j)]
    keep <- dmin <= max_dist_m
    rows <- grp$.row[keep]
    out$tow_id[rows]        <- cand$tow_id[j[keep]]
    out$dist_to_tow_m[rows] <- dmin[keep]
  }
  out
}

# ---- Point -> group HC/SC % aggregator ----
# Unclassified ("IN") points are excluded from the denominator.
aggregate_hc_sc_inline <- function(points, group_col) {
  points |>
    filter(!is.na(.data[[group_col]])) |>
    group_by(.data[[group_col]]) |>
    summarise(
      n_points     = n(),
      n_classified = sum(point_machine_group_code != "IN", na.rm = TRUE),
      pct_HC       = sum(point_machine_group_code == "HC", na.rm = TRUE) / n_classified,
      pct_SC       = sum(point_machine_group_code == "SC", na.rm = TRUE) / n_classified,
      .groups      = "drop"
    )
}

# Where coral cover exists from the same imagery, a missing detection count means zero.
zero_fill_when_observed_inline <- function(x, observed) if_else(observed & is.na(x), 0, x)

# ---- COTSpotter counts (built by 00a_build_cotspotter_detections.R from the ingested tows) ----
{
source("scripts/00a_build_cotspotter_detections.R", local = new.env())
rs_surface_cots      <- read_csv("data/derived/cotspotter_surface_tow_counts.csv", show_col_types = FALSE) |> mutate(tow_id = as.character(tow_id))
rsd_cots_tow         <- read_csv("data/derived/cotspotter_deep_tow_counts.csv", show_col_types = FALSE) |> mutate(tow_id = as.character(tow_id))
rsd_cots_cull        <- read_csv("data/derived/cotspotter_deep_cull_site_counts.csv", show_col_types = FALSE)
rsd_cots_reef        <- read_csv("data/derived/cotspotter_deep_reef_counts.csv", show_col_types = FALSE)
rs_surface_cots_cull <- read_csv("data/derived/cotspotter_surface_cull_site_counts.csv", show_col_types = FALSE)
rs_surface_cots_reef <- read_csv("data/derived/cotspotter_surface_reef_counts.csv", show_col_types = FALSE)

## ---- Cull-site/reef counts excluding dropped voyages ----
## The pre-built cull-site/reef counts carry no voyage or date, so re-derive
## them from the detection-level files (trip_number matches the manta voyage
## label), using the same aggregation as 00a_build_cotspotter_detections.R.
rebuild_cotspotter_cull_reef <- function(detections_path, col_name) {
  det <- read_csv(detections_path, show_col_types = FALSE) |>
    filter(!invalid_coordinates, !trip_number %in% exclude_voyages)
  list(
    cull = det |> filter(!is.na(cull_site)) |> count(cull_site, name = col_name),
    reef = det |> count(reef_label, name = col_name)
  )
}

surf_cotspotter      <- rebuild_cotspotter_cull_reef("data/derived/cotspotter_surface_detections.csv", "cots_reefscan_surface")
deep_cotspotter      <- rebuild_cotspotter_cull_reef("data/derived/cotspotter_deep_detections.csv", "cots_reefscan_deep")
rs_surface_cots_cull <- surf_cotspotter$cull
rs_surface_cots_reef <- surf_cotspotter$reef
rsd_cots_cull        <- deep_cotspotter$cull
rsd_cots_reef        <- deep_cotspotter$reef
rm(surf_cotspotter, deep_cotspotter, rebuild_cotspotter_cull_reef)
}

# ---- Tow level ----
# RS Surface joined on tow_id; RS Deep via nearest-tow match (25 m cutoff);
# cull_site from manta's own column, with a spatial fallback for gaps.
{
rs_surf_tow <- aggregate_hc_sc_inline(rs_surface_raw, "tow_id") |>
  rename(pct_HC_surf = pct_HC, pct_SC_surf = pct_SC, n_pts_surf = n_points) |>
  select(-n_classified)

tow_level <- mt |>
  mutate(tow_mid_lon = (start_lon + end_lon) / 2, tow_mid_lat = (start_lat + end_lat) / 2) |>
  left_join(rs_surf_tow, by = "tow_id")

needs_fallback <- is.na(tow_level$cull_site)
if (any(needs_fallback)) {
  fallback <- assign_cull_site_inline(tow_level$tow_mid_lon[needs_fallback], tow_level$tow_mid_lat[needs_fallback], cull_polygons)
  tow_level$cull_site[needs_fallback] <- fallback$cull_site
}

deep_site_to_reef <- reef_xwalk |>
  filter(has_rs_deep) |>
  separate_rows(rs_deep_sites, sep = "; ") |>
  select(rs_deep_site = rs_deep_sites, reef_label)

rs_deep_labeled <- rs_deep_raw |> inner_join(deep_site_to_reef, by = c("site_reef_name" = "rs_deep_site"))

deep_images <- rs_deep_labeled |> distinct(site_reef_name, image_latitude, image_longitude, reef_label)
deep_match  <- assign_nearest_tow_inline(deep_images$image_longitude, deep_images$image_latitude, deep_images$reef_label, manta_lines, max_dist_m = 25)
deep_images$tow_id <- deep_match$tow_id

rs_deep_tow <- rs_deep_labeled |>
  left_join(deep_images |> select(site_reef_name, image_latitude, image_longitude, tow_id), by = c("site_reef_name", "image_latitude", "image_longitude")) |>
  aggregate_hc_sc_inline("tow_id") |>
  rename(pct_HC_deep = pct_HC, pct_SC_deep = pct_SC, n_pts_deep = n_points) |>
  select(-n_classified)

tow_level <- tow_level |>
  left_join(rs_deep_tow, by = "tow_id") |>
  left_join(rs_surface_cots, by = "tow_id") |>
  left_join(rsd_cots_tow, by = "tow_id")

tow_level$cots_reefscan_surface <- zero_fill_when_observed_inline(tow_level$cots_reefscan_surface, !is.na(tow_level$pct_HC_surf))
tow_level$cots_reefscan_deep    <- zero_fill_when_observed_inline(tow_level$cots_reefscan_deep, !is.na(tow_level$pct_HC_deep))

tow_level <- tow_level |>
  select(
    tow_id, tow_date, tow_datetime, voyage, vessel, survey_folder, source_file,
    reef_label, reef_name, cull_site, dist_cull_m,
    cots_observed, feeding_scars, scar_present, scar_count, scar_count_imputed,
    hard_coral, soft_coral, recently_dead, bleach_sev, bleach_ext,
    surveyor_id, surveyor_level, tow_distance_m,
    tow_mid_lon, tow_mid_lat, perc_hard, perc_soft,
    pct_HC_surf, pct_SC_surf, n_pts_surf,
    pct_HC_deep, pct_SC_deep, n_pts_deep,
    cots_reefscan_surface, cots_reefscan_deep
  ) |>
  ## has_manta_tow: feeding scar present OR >= 1 COTS observed (response in manta-as-response COTS models).
  mutate(has_manta_tow = as.integer(scar_present == 1 | cots_observed > 0), .after = scar_present)

rm(rs_surf_tow, needs_fallback, deep_site_to_reef, rs_deep_labeled, deep_images, deep_match, rs_deep_tow)
cat("tow_level:", nrow(tow_level), "rows,", n_distinct(tow_level$reef_label), "reefs\n")
}

# ---- Cull-site level ----
# manta / RS Surface from tow_level's cull_site; RS Deep and RS Surface points
# assigned to cull sites directly.
{
manta_surf_cull <- tow_level |>
  filter(!is.na(cull_site)) |>
  group_by(cull_site) |>
  summarise(
    reef_label                       = first(reef_label),
    n_tows                           = n(),
    total_tow_distance_m             = sum(tow_distance_m, na.rm = TRUE),
    perc_hard_se                     = wtd_mean_se(perc_hard, tow_distance_m),
    perc_hard_simple_se              = mean_se(perc_hard),
    perc_hard_simple                 = mean(perc_hard, na.rm = TRUE),
    perc_hard                        = weighted.mean(perc_hard, tow_distance_m, na.rm = TRUE),
    perc_soft_se                     = wtd_mean_se(perc_soft, tow_distance_m),
    perc_soft_simple_se              = mean_se(perc_soft),
    perc_soft_simple                 = mean(perc_soft, na.rm = TRUE),
    perc_soft                        = weighted.mean(perc_soft, tow_distance_m, na.rm = TRUE),
    cots_observed                    = sum(cots_observed, na.rm = TRUE),
    cots_observed_se                 = pois_sum_se(cots_observed),
    feeding_scars                    = sum(scar_present, na.rm = TRUE),
    scar_count                       = sum(scar_count, na.rm = TRUE),
    n_scar_tows_imputed              = sum(scar_count_imputed, na.rm = TRUE),
    manta_combined                   = sum(has_manta_tow, na.rm = TRUE),
    n_tows_with_surf                 = sum(!is.na(pct_HC_surf)),
    pct_HC_surf_se                   = mean_se(pct_HC_surf),
    pct_HC_surf                      = mean(pct_HC_surf, na.rm = TRUE),
    pct_SC_surf_se                   = mean_se(pct_SC_surf),
    pct_SC_surf                      = mean(pct_SC_surf, na.rm = TRUE),
    n_tows_reefscan_surface_cots_obs = sum(!is.na(cots_reefscan_surface)),
    .groups = "drop"
  ) |>
  mutate(
    scar_rate    = feeding_scars / n_tows,
    scar_rate_se = binom_se(scar_rate, n_tows),
    manta_combined_rate    = manta_combined / n_tows,
    manta_combined_rate_se = binom_se(manta_combined_rate, n_tows),
    cots_rate    = cots_observed / n_tows,
    cots_rate_se = pois_rate_se(cots_observed, n_tows)
  )

deep_assigned_cull <- assign_cull_site_inline(rs_deep_raw$image_longitude, rs_deep_raw$image_latitude, cull_polygons)
rs_deep_cull <- rs_deep_raw |>
  mutate(cull_site = deep_assigned_cull$cull_site) |>
  aggregate_hc_sc_inline("cull_site") |>
  rename(pct_HC_deep = pct_HC, pct_SC_deep = pct_SC, n_pts_deep = n_points, n_classified_deep = n_classified) |>
  mutate(pct_HC_deep_se = binom_se(pct_HC_deep, n_classified_deep))

surf_assigned_cull <- assign_cull_site_inline(rs_surface_raw$image_longitude, rs_surface_raw$image_latitude, cull_polygons)
rs_surf_cull_pts <- rs_surface_raw |>
  mutate(cull_site = surf_assigned_cull$cull_site) |>
  aggregate_hc_sc_inline("cull_site") |>
  select(cull_site, n_pts_surf = n_points)

# RHIS: restricted to manta reefs; its own CullSite column is authoritative,
# with a spatial fallback for gaps.
rhis_scope <- rhis |> filter(reef_label %in% reef_xwalk$reef_label)
needs_fallback_rhis <- is.na(rhis_scope$cull_site)
rhis_cull_site <- rhis_scope$cull_site
if (any(needs_fallback_rhis)) {
  fallback_rhis <- assign_cull_site_inline(rhis_scope$longitude[needs_fallback_rhis], rhis_scope$latitude[needs_fallback_rhis], cull_polygons)
  rhis_cull_site[needs_fallback_rhis] <- fallback_rhis$cull_site
}

rhis_cull <- rhis_scope |>
  mutate(cull_site = rhis_cull_site) |>
  filter(!is.na(cull_site)) |>
  group_by(cull_site) |>
  summarise(
    n_surveys_rhis          = n(),
    pct_HC_rhis_se          = mean_se(hc_prop_survey),
    pct_HC_rhis             = mean(hc_prop_survey, na.rm = TRUE),
    pct_SC_rhis_se          = mean_se(sc_prop_survey),
    pct_SC_rhis             = mean(sc_prop_survey, na.rm = TRUE),
    n_surveys_predator_rhis = sum(has_predator_survey),
    cots_observed_rhis      = sum(cots_count_survey, na.rm = TRUE),
    scar_count_rhis         = sum(scar_count_survey, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    cots_observed_rhis_se = pois_sum_se(cots_observed_rhis),
    cots_rate_rhis         = cots_observed_rhis / n_surveys_rhis,
    cots_rate_rhis_se     = pois_rate_se(cots_observed_rhis, n_surveys_rhis),
    scar_count_rhis_se    = pois_sum_se(scar_count_rhis),
    scar_rate_rhis        = scar_count_rhis / n_surveys_rhis,
    scar_rate_rhis_se     = pois_rate_se(scar_count_rhis, n_surveys_rhis)
  )

cull_site_level <- manta_surf_cull |>
  full_join(rs_deep_cull, by = "cull_site") |>
  full_join(rs_surf_cull_pts, by = "cull_site") |>
  full_join(rhis_cull, by = "cull_site") |>
  arrange(cull_site) |>
  left_join(rsd_cots_cull, by = "cull_site") |>
  left_join(rs_surface_cots_cull, by = "cull_site")

cull_site_level$cots_reefscan_deep    <- zero_fill_when_observed_inline(cull_site_level$cots_reefscan_deep, !is.na(cull_site_level$n_pts_deep) & cull_site_level$n_pts_deep > 0)
cull_site_level$cots_reefscan_surface <- zero_fill_when_observed_inline(cull_site_level$cots_reefscan_surface, !is.na(cull_site_level$n_tows_with_surf) & cull_site_level$n_tows_with_surf > 0)
cull_site_level$cots_reefscan_deep_se    <- pois_sum_se(cull_site_level$cots_reefscan_deep)
cull_site_level$cots_reefscan_surface_se <- pois_sum_se(cull_site_level$cots_reefscan_surface)

rm(manta_surf_cull, deep_assigned_cull, rs_deep_cull, surf_assigned_cull, rs_surf_cull_pts,
   rhis_scope, needs_fallback_rhis, rhis_cull_site, rhis_cull)
if (exists("fallback_rhis")) rm(fallback_rhis)
cat("cull_site_level:", nrow(cull_site_level), "rows\n")
}

# ---- Reef level ----
# manta / RS Surface from tow_level means; RS points via the crosswalk; eDNA and RHIS joined here.
{
manta_reef <- tow_level |>
  group_by(reef_label) |>
  summarise(
    n_tows                            = n(),
    mean_HC_manta                     = mean(perc_hard, na.rm = TRUE),
    mean_HC_manta_se                  = mean_se(perc_hard),
    cots_observed                     = sum(cots_observed, na.rm = TRUE),
    cots_observed_se                  = pois_sum_se(cots_observed),
    feeding_scars                     = sum(scar_present, na.rm = TRUE),
    scar_count                        = sum(scar_count, na.rm = TRUE),
    n_scar_tows_imputed               = sum(scar_count_imputed, na.rm = TRUE),
    manta_combined                    = sum(has_manta_tow, na.rm = TRUE),
    n_tows_reefscan_surface_cots_obs  = sum(!is.na(cots_reefscan_surface)),
    .groups = "drop"
  ) |>
  mutate(
    scar_rate    = feeding_scars / n_tows,
    scar_rate_se = binom_se(scar_rate, n_tows),
    manta_combined_rate    = manta_combined / n_tows,
    manta_combined_rate_se = binom_se(manta_combined_rate, n_tows),
    cots_rate    = cots_observed / n_tows,
    cots_rate_se = pois_rate_se(cots_observed, n_tows)
  )

surf_reef <- tow_level |>
  filter(!is.na(pct_HC_surf)) |>
  group_by(reef_label) |>
  summarise(
    n_tows_with_surf = n(),
    mean_HC_surf     = mean(pct_HC_surf, na.rm = TRUE),
    mean_HC_surf_se  = mean_se(pct_HC_surf),
    .groups = "drop"
  )

deep_site_to_reef_r <- reef_xwalk |>
  filter(has_rs_deep) |>
  separate_rows(rs_deep_sites, sep = "; ") |>
  select(rs_deep_site = rs_deep_sites, reef_label)

deep_reef <- rs_deep_raw |>
  filter(!is.na(site_reef_name)) |>
  inner_join(deep_site_to_reef_r, by = c("site_reef_name" = "rs_deep_site")) |>
  aggregate_hc_sc_inline("reef_label") |>
  transmute(
    reef_label,
    n_pts_deep_reef = n_points,
    mean_HC_deep    = pct_HC,
    mean_HC_deep_se = binom_se(pct_HC, n_classified)
  )

surf_site_to_reef_r <- reef_xwalk |>
  filter(has_rs_surf) |>
  separate_rows(rs_surf_site, sep = "; ") |>
  select(rs_surf_site, reef_label)

surf_reef_pts <- rs_surface_raw |>
  filter(!is.na(site_reef_name)) |>
  group_by(site_reef_name) |>
  summarise(n_pts_surf_reef = n(), .groups = "drop") |>
  inner_join(surf_site_to_reef_r, by = c("site_reef_name" = "rs_surf_site")) |>
  group_by(reef_label) |>
  summarise(n_pts_surf_reef = sum(n_pts_surf_reef), .groups = "drop")

edna_reef <- edna_reef_raw |>
  left_join(reef_xwalk |> filter(has_edna) |> select(edna_reef, reef_label), by = c("Reef" = "edna_reef")) |>
  filter(!is.na(reef_label)) |>
  transmute(
    reef_label,
    edna_prop = ProportionPositive, edna_prop_lower_hpd = lower.HPD, edna_prop_upper_hpd = upper.HPD,
    edna_conc = Concentration, edna_conc_se = SE,
    edna_conc_asymp_lcl = asymp.LCL, edna_conc_asymp_ucl = asymp.UCL
  )

rhis_reef <- rhis |>
  filter(reef_label %in% reef_xwalk$reef_label) |>
  group_by(reef_label) |>
  summarise(
    n_surveys_rhis          = n(),
    pct_HC_rhis_se          = mean_se(hc_prop_survey),
    pct_HC_rhis             = mean(hc_prop_survey, na.rm = TRUE),
    pct_SC_rhis_se          = mean_se(sc_prop_survey),
    pct_SC_rhis             = mean(sc_prop_survey, na.rm = TRUE),
    n_surveys_predator_rhis = sum(has_predator_survey),
    cots_observed_rhis      = sum(cots_count_survey, na.rm = TRUE),
    scar_count_rhis         = sum(scar_count_survey, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    cots_observed_rhis_se = pois_sum_se(cots_observed_rhis),
    cots_rate_rhis         = cots_observed_rhis / n_surveys_rhis,
    cots_rate_rhis_se     = pois_rate_se(cots_observed_rhis, n_surveys_rhis),
    scar_count_rhis_se    = pois_sum_se(scar_count_rhis),
    scar_rate_rhis        = scar_count_rhis / n_surveys_rhis,
    scar_rate_rhis_se     = pois_rate_se(scar_count_rhis, n_surveys_rhis)
  )

reef_level <- manta_reef |>
  left_join(surf_reef, by = "reef_label") |>
  left_join(deep_reef, by = "reef_label") |>
  left_join(surf_reef_pts, by = "reef_label") |>
  left_join(edna_reef, by = "reef_label") |>
  left_join(rhis_reef, by = "reef_label") |>
  mutate(has_surf = !is.na(mean_HC_surf), has_deep = !is.na(mean_HC_deep), has_edna = !is.na(edna_prop), has_rhis = !is.na(pct_HC_rhis)) |>
  relocate(reef_label, n_tows, n_tows_with_surf, mean_HC_manta, mean_HC_surf, mean_HC_deep) |>
  left_join(rsd_cots_reef, by = "reef_label") |>
  left_join(rs_surface_cots_reef, by = "reef_label")

reef_level$cots_reefscan_deep    <- zero_fill_when_observed_inline(reef_level$cots_reefscan_deep, !is.na(reef_level$n_pts_deep_reef) & reef_level$n_pts_deep_reef > 0)
reef_level$cots_reefscan_surface <- zero_fill_when_observed_inline(reef_level$cots_reefscan_surface, !is.na(reef_level$n_tows_with_surf) & reef_level$n_tows_with_surf > 0)
reef_level$cots_reefscan_deep_se    <- pois_sum_se(reef_level$cots_reefscan_deep)
reef_level$cots_reefscan_surface_se <- pois_sum_se(reef_level$cots_reefscan_surface)

rm(manta_reef, surf_reef, deep_site_to_reef_r, deep_reef, surf_site_to_reef_r, surf_reef_pts, edna_reef, rhis_reef)
cat("reef_level:", nrow(reef_level), "rows -- has_surf:", sum(reef_level$has_surf),
    " has_deep:", sum(reef_level$has_deep), " has_edna:", sum(reef_level$has_edna),
    " has_rhis:", sum(reef_level$has_rhis), "\n")
}

# ---- Save canonical tables ----
{
saveRDS(tow_level, "data/derived/tow_level.RDS")
saveRDS(cull_site_level, "data/derived/cull_site_level.RDS")
saveRDS(reef_level, "data/derived/reef_level.RDS")

# Checkpoint intermediates reused by 02_cots_calibration.R.
saveRDS(mt, "data/derived/mt.RDS")
saveRDS(edna_raw, "data/derived/edna_raw.RDS")
saveRDS(rhis, "data/derived/rhis.RDS")
saveRDS(reef_xwalk, "data/derived/reef_xwalk.RDS")
}
