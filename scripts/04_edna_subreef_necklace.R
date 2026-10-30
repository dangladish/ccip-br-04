########################################################################
## 04_edna_subreef_necklace.R
##
## Sub-reef eDNA maps (descriptive only; fits no model). For each reef with
## both eDNA sites and manta tows, builds one continuous chain ("necklace")
## of ~10-ha beads through the eDNA sites and along the manta tow track, and
## maps eDNA concentration, eDNA proportion positive, manta density and
## ReefScan Surface density bead by bead over satellite imagery.
##
## Inputs:  data/derived/{tow_level,edna_raw,reef_xwalk}.RDS (Esri World
##          Imagery tiles are downloaded on first run, then cached).
## Outputs: output/figures/edna_necklace/necklace_<reef>.png
##          data/derived/edna_subreef_beads.RDS (bead table, no geometry)
##          data/derived/edna_subreef_summary.csv (one row per reef)
##          output/figures/edna_necklace_summary.png (only if make_summary_figure)
##
## Per reef:
##   1. eDNA sites -> smooth path through them (open path or closed loop),
##      extended/rotated so each site sits mid-bead.
##   2. Manta tows, in time order -> smoothed ribbons kept clear of the chain.
##   3. Ribbons joined to the chain by gentle Bezier links or side branches.
##   4. Each piece is a constant-width (200 m) strip cut into ~500 m beads
##      by Voronoi dividers, so beads never overlap.
##   5. Tows and eDNA sites are assigned to beads; per-bead summaries computed.
########################################################################

library(tidyverse)
library(sf)
library(maptiles)
library(tidyterra)
library(patchwork)
library(magick)

## ---- Inputs ----
tow_level  <- readRDS("data/derived/tow_level.RDS")
edna_raw   <- readRDS("data/derived/edna_raw.RDS")
reef_xwalk <- readRDS("data/derived/reef_xwalk.RDS")

fig_dir   <- "output/figures/edna_necklace"
tile_dir  <- "data/derived/_tile_cache"
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tile_dir, showWarnings = FALSE, recursive = TRUE)

## Map eDNA reef codes to manta reef_label
edna_mapped <- edna_raw |>
  inner_join(
    reef_xwalk |> filter(!is.na(edna_reef)) |> select(edna_reef, reef_label),
    by = "edna_reef"
  )

target_area_m2 <- 100000   # 10 ha, cull-site size
tile_zoom      <- 15
fig_w_max_cm   <- 15.9     # A4 text width
fig_h_max_cm   <- 11.5     # so two figures fit on an A4 page
make_summary_figure <- FALSE   # TRUE also builds the one-page summary figure

path_cum <- \(P) c(0, cumsum(sqrt(rowSums((P[-1, , drop = FALSE] - P[-nrow(P), , drop = FALSE])^2))))
interp_pt <- \(P, cum, u) cbind(approx(cum, P[, 1], u)$y, approx(cum, P[, 2], u)$y)
turn_max <- \(Q) { dd <- diff(Q); aa <- atan2(dd[, 2], dd[, 1]); tt <- abs(diff(aa)); max(pmin(tt, 2 * pi - tt)) * 180 / pi }

## ---- eDNA chain path ----
## Smooth path through the sites; open paths are extended past the end sites (closed paths rotated) so sites sit mid-bead.
chain_path <- function(sites, L = 500, ds = 25, loop_ratio = 1.6, n_iter = 5, ext_range = c(120, 380)) {
  xy <- st_coordinates(sites)[, 1:2, drop = FALSE]; n <- nrow(xy)
  d <- as.matrix(dist(xy))
  perms <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), \(i) lapply(perms(v[-i]), \(p) c(v[i], p))))
  allp <- perms(seq_len(n))
  open_len <- vapply(allp, \(p) sum(d[cbind(p[-n], p[-1])]), 1)
  closed_len <- vapply(allp, \(p) sum(d[cbind(p, c(p[-1], p[1]))]), 1)
  closed <- n >= 3 && min(closed_len) < loop_ratio * min(open_len)
  C <- xy[if (closed) allp[[which.min(closed_len)]] else allp[[which.min(open_len)]], , drop = FALSE]
  subdiv <- function(Pm, closed) {
    m <- nrow(Pm)
    get <- \(i) if (closed) Pm[((i - 1) %% m) + 1, ] else {
      if (i < 1) 2 * Pm[1, ] - Pm[min(2 - i, m), ] else if (i > m) 2 * Pm[m, ] - Pm[max(2 * m - i, 1), ] else Pm[i, ] }
    nseg <- if (closed) m else m - 1
    out <- vector("list", 2 * nseg + (if (closed) 0 else 1)); k <- 0
    for (i in seq_len(nseg)) {
      k <- k + 1; out[[k]] <- Pm[i, ]
      k <- k + 1; out[[k]] <- (-1 / 16) * (get(i - 1) + get(i + 2)) + (9 / 16) * (get(i) + get(i + 1))
    }
    if (!closed) { k <- k + 1; out[[k]] <- Pm[m, ] }
    do.call(rbind, out)
  }
  Pm <- C
  for (it in seq_len(n_iter)) Pm <- subdiv(Pm, closed)
  if (closed) Pm <- rbind(Pm, Pm[1, ])
  cum <- path_cum(Pm); ell0 <- max(cum)
  g <- seq(0, ell0, by = ds); if (tail(g, 1) < ell0 - 1e-6) g <- c(g, ell0)
  P <- interp_pt(Pm, cum, g); m <- nrow(P)
  if (closed) {
    P <- P[-m, , drop = FALSE]; m <- nrow(P)
    sm <- \(v, w = 3) { idx <- outer(seq_len(m), -w:w, \(i, j) ((i + j - 1) %% m) + 1); rowMeans(matrix(v[idx], m)) }
    P <- cbind(sm(P[, 1]), sm(P[, 2])); P <- rbind(P, P[1, ])
  }
  cum <- path_cum(P); ell0 <- max(cum)
  u0 <- vapply(seq_len(n), \(i) cum[which.min((P[, 1] - xy[i, 1])^2 + (P[, 2] - xy[i, 2])^2)], 1)
  if (closed) {
    nb <- max(6, round(ell0 / L)); spacing <- ell0 / nb
    phis <- seq(0, spacing - 1, by = 10)
    marg <- vapply(phis, \(ph) { cuts <- (ph + (0:(nb - 1)) * spacing) %% ell0
      min(vapply(u0, \(u) { dd <- abs(u - cuts); min(pmin(dd, ell0 - dd)) }, 1)) }, 1)
    ph <- phis[which.max(marg)]
    return(list(P = P, ctr = (ph + (seq_len(nb) - 0.5) * spacing) %% ell0, closed = TRUE, ell = ell0))
  }
  e1 <- P[1, ] - P[3, ]; e1 <- e1 / sqrt(sum(e1^2)); e2 <- P[m, ] - P[m - 2, ]; e2 <- e2 / sqrt(sum(e2^2))
  grid <- expand.grid(es = seq(ext_range[1], ext_range[2], by = 20), ee = seq(ext_range[1], ext_range[2], by = 20))
  sc <- vapply(seq_len(nrow(grid)), \(i) {
    len <- ell0 + grid$es[i] + grid$ee[i]; nb <- max(1, round(len / L)); sp <- len / nb
    cuts <- (0:nb) * sp; mg <- min(vapply(u0 + grid$es[i], \(u) min(abs(u - cuts)), 1))
    mg - 0.5 * abs(sp - L) }, 1)
  bi <- which.max(sc); es <- grid$es[bi]; ee <- grid$ee[bi]
  ext_pts <- \(p0, e, len) { dd <- seq(0, len, length.out = max(2, ceiling(len / ds) + 1)); cbind(p0[1] + e[1] * dd, p0[2] + e[2] * dd) }
  head_ext <- ext_pts(P[1, ], e1, es); head_ext <- head_ext[nrow(head_ext):2, , drop = FALSE]
  tail_ext <- ext_pts(P[m, ], e2, ee)[-1, , drop = FALSE]
  Pext <- rbind(head_ext, P, tail_ext)
  ell <- max(path_cum(Pext)); nb <- max(1, round(ell / L)); spacing <- ell / nb
  list(P = Pext, ctr = (seq_len(nb) - 0.5) * spacing, closed = FALSE, ell = ell)
}

## ---- Tow ribbons ----
tow_paths <- function(tows, avoid_lines, dt, ds = 25, break_m = 800, min_len = 1000, sep_m = 230, max_turn = 10) {
  crs <- st_crs(tows)
  xy <- st_coordinates(tows)[, 1:2, drop = FALSE]
  o <- order(dt, na.last = TRUE); xy <- xy[o, , drop = FALSE]
  d <- sqrt(rowSums((xy[-1, , drop = FALSE] - xy[-nrow(xy), , drop = FALSE])^2))
  segs <- split(seq_len(nrow(xy)), c(1L, 1L + cumsum(d > break_m)))
  segs <- segs[lengths(segs) >= 3]
  build_path <- function(idx) {
    P <- xy[idx, , drop = FALSE]; k <- nrow(P); hk <- min(8, floor((k - 1) / 2))
    if (k >= 5 && hk >= 2) {
      w <- dnorm(-hk:hk, sd = hk / 2.3); w <- w / sum(w)
      pad <- \(v) c(2 * v[1] - v[(hk + 1):2], v, 2 * v[k] - v[(k - 1):(k - hk)])
      Ps <- cbind(as.numeric(stats::filter(pad(P[, 1]), w))[(hk + 1):(k + hk)], as.numeric(stats::filter(pad(P[, 2]), w))[(hk + 1):(k + hk)])
      Ps[1, ] <- P[1, ]; Ps[k, ] <- P[k, ]
    } else Ps <- P
    cum <- path_cum(Ps); ell <- max(cum)
    if (ell < min_len) return(NULL)
    g <- seq(0, ell, by = ds); if (tail(g, 1) < ell - 1e-6) g <- c(g, ell)
    Q <- interp_pt(Ps, cum, g); m <- nrow(Q)
    if (m >= 5) for (it in 1:400) {
      if (turn_max(Q) <= max_turn) break
      for (j in 1:5) { Q2 <- Q; Q2[2:(m - 1), ] <- (Q[1:(m - 2), ] + 2 * Q[2:(m - 1), ] + Q[3:m, ]) / 4; Q <- Q2 }
    }
    Q
  }
  paths <- Filter(Negate(is.null), lapply(segs, build_path))
  paths <- paths[order(-vapply(paths, \(Q) max(path_cum(Q)), 1))]
  avoid <- avoid_lines                            # sfc of lines already placed
  out <- list()
  for (Q in paths) {
    pts <- st_as_sf(data.frame(x = Q[, 1], y = Q[, 2]), coords = c("x", "y"), crs = crs)
    far <- as.numeric(st_distance(st_geometry(pts), avoid)) >= sep_m
    r <- rle(far); ends <- cumsum(r$lengths); begs <- ends - r$lengths + 1
    for (q in which(r$values)) {
      sub <- Q[begs[q]:ends[q], , drop = FALSE]
      if (nrow(sub) >= 3 && max(path_cum(sub)) >= min_len) {
        out[[length(out) + 1]] <- sub
        avoid <- st_union(avoid, st_sfc(st_linestring(sub), crs = crs))
      }
    }
  }
  out
}

## ---- Join ribbons to the chain ----
## Nearest ribbon first: (1) Bezier end join, (2) hairpin end join with a gentle-curve limit,
## (3) branch to the side of the necklace (the only option for a closed loop).
## A ribbon no stage can attach stays a separate piece.
link_bezier <- function(p0, tn, p3, tr, k, ds) {
  dd <- sqrt(sum((p3 - p0)^2)); p1 <- p0 + tn * k; p2 <- p3 + tr * k
  tt <- seq(0, 1, length.out = 400)
  B <- t(vapply(tt, \(t) (1 - t)^3 * p0 + 3 * (1 - t)^2 * t * p1 + 3 * (1 - t) * t^2 * p2 + t^3 * p3, c(0, 0)))
  cum <- path_cum(B); g <- seq(0, max(cum), length.out = max(4, round(max(cum) / ds) + 1))
  interp_pt(B, cum, g)
}
clear_of <- function(B, rest, min_d = 190, skip_m = 80) {
  if (is.null(rest) || nrow(rest) == 0) return(TRUE)
  cb <- path_cum(B); inner <- B[cb > skip_m & cb < max(cb) - skip_m, , drop = FALSE]
  if (nrow(inner) == 0) return(TRUE)
  dx <- outer(inner[, 1], rest[, 1], "-"); dy <- outer(inner[, 2], rest[, 2], "-")
  sqrt(min(dx^2 + dy^2)) >= min_d
}
join_pieces <- function(base, ribbons, L = 500, ds = 25, max_link = 2500, max_angle = 135,
                        hw = 100, max_turn_link = 14, near_m = 260) {
  unit <- \(v) v / sqrt(sum(v^2))
  tang_out <- \(P, end) if (end == "start") unit(P[1, ] - P[min(4, nrow(P)), ]) else unit(P[nrow(P), ] - P[max(nrow(P) - 3, 1), ])
  ang <- \(a, b) acos(max(-1, min(1, sum(a * b)))) * 180 / pi
  eq_ctr <- \(P, nmin = 1) { len <- max(path_cum(P)); nb <- max(nmin, round(len / L)); (seq_len(nb) - 0.5) * len / nb }
  concat <- function(sec) {
    Pall <- NULL; ctr_all <- numeric(0); type_all <- character(0)
    for (s_i in sec) {
      Pi <- s_i$P; if (!is.null(Pall)) Pi <- Pi[-1, , drop = FALSE]
      off <- if (is.null(Pall)) 0 else max(path_cum(Pall))
      Pall <- rbind(Pall, Pi); ctr_all <- c(ctr_all, off + s_i$ctr); type_all <- c(type_all, rep(s_i$type, length(s_i$ctr)))
    }
    list(P = Pall, ctr = ctr_all, type = type_all)
  }
  sec <- list(list(P = base$P, ctr = base$ctr, type = "eDNA"))
  branches <- list(); rem <- ribbons
  while (length(rem) > 0) {
    cur <- concat(sec); Pn <- cur$P; cumn <- path_cum(Pn); total <- max(cumn)
    best <- NULL
    # stages 1-2: end joins (open necklace only)
    if (!base$closed) {
      for (stage in 1:2) {
        for (ri in seq_along(rem)) {
          R <- rem[[ri]]
          for (re in c("start", "end")) {
            pr <- if (re == "start") R[1, ] else R[nrow(R), ]; tr_ <- tang_out(R, re)
            for (ne in c("start", "end")) {
              pn <- if (ne == "start") Pn[1, ] else Pn[nrow(Pn), ]; tn <- tang_out(Pn, ne)
              dd <- sqrt(sum((pr - pn)^2)); if (dd > max_link || dd < 1) next
              dirv <- unit(pr - pn); a1 <- ang(tn, dirv); a2 <- ang(tr_, -dirv)
              if (stage == 1) {
                if (a1 > max_angle || a2 > max_angle) next
                B <- link_bezier(pn, tn, pr, tr_, 0.4 * dd, ds)
                score <- dd * (1 + (a1 + a2) / 180)
              } else {
                if (a1 > 175 || a2 > 175) next
                B <- NULL
                for (km in seq(0.3, 1.6, by = 0.05)) {
                  Bk <- link_bezier(pn, tn, pr, tr_, km * dd, ds)
                  if (turn_max(Bk) <= max_turn_link) { B <- Bk; break }
                }
                if (is.null(B)) next
                score <- max(path_cum(B))
                rest <- if (ne == "start") Pn[cumn > near_m, , drop = FALSE] else Pn[cumn < total - near_m, , drop = FALSE]
                if (!clear_of(B, rest)) next
              }
              if (is.null(best) || score < best$score) best <- list(kind = "end", score = score, ri = ri, re = re, ne = ne, B = B)
            }
          }
        }
        if (!is.null(best)) break
      }
    }
    # stage 3: branch join
    if (is.null(best)) {
      for (ri in seq_along(rem)) {
        R <- rem[[ri]]
        for (re in c("start", "end")) {
          pr <- if (re == "start") R[1, ] else R[nrow(R), ]; tr_ <- tang_out(R, re)
          k <- which.min((Pn[, 1] - pr[1])^2 + (Pn[, 2] - pr[2])^2)
          i0 <- max(1, k - 3); i1 <- min(nrow(Pn), k + 3); tg <- unit(Pn[i1, ] - Pn[i0, ]); nm <- c(-tg[2], tg[1])
          if (sum(nm * (pr - Pn[k, ])) < 0) nm <- -nm
          e <- Pn[k, ] + (hw + 2) * nm; dd <- sqrt(sum((pr - e)^2)); if (dd > max_link) next
          if (dd < 100) {
            B <- rbind(e, pr)
          } else {
            B <- NULL
            for (km in seq(0.3, 1.6, by = 0.05)) {
              Bk <- link_bezier(e, nm, pr, tr_, km * dd, ds)
              if (turn_max(Bk) <= max_turn_link) { B <- Bk; break }
            }
            if (is.null(B)) next
          }
          near_q <- sqrt((Pn[, 1] - Pn[k, 1])^2 + (Pn[, 2] - Pn[k, 2])^2) < 450
          if (!clear_of(B, Pn[!near_q, , drop = FALSE])) next
          score <- sum(sqrt(rowSums((B[-1, , drop = FALSE] - B[-nrow(B), , drop = FALSE])^2)))
          if (is.null(best) || score < best$score) best <- list(kind = "branch", score = score, ri = ri, re = re, B = B)
        }
      }
    }
    if (is.null(best)) break
    R <- rem[[best$ri]]
    if (best$re == "end") R <- R[nrow(R):1, , drop = FALSE]               # R now starts at the joined end
    B <- best$B
    if (best$kind == "end") {
      link <- list(P = B, ctr = eq_ctr(B), type = "link"); rib <- list(P = R, ctr = eq_ctr(R, 2), type = "tow")
      if (best$ne == "end") sec <- c(sec, list(link, rib)) else {
        link$P <- link$P[nrow(link$P):1, , drop = FALSE]; link$ctr <- max(path_cum(link$P)) - rev(link$ctr)
        rib$P <- rib$P[nrow(rib$P):1, , drop = FALSE]; rib$ctr <- max(path_cum(rib$P)) - rev(rib$ctr)
        sec <- c(list(rib, link), sec)
      }
    } else {
      # branch: one piece = link (from the necklace side) + ribbon
      Pb <- rbind(B, R[-1, , drop = FALSE]); len_link <- max(path_cum(B)); ctr <- eq_ctr(Pb, 2)
      branches[[length(branches) + 1]] <- list(P = Pb, ctr = ctr, type = ifelse(ctr < len_link, "link", "tow"))
    }
    rem <- rem[-best$ri]
  }
  cur <- concat(sec)
  list(path = if (base$closed) NULL else list(P = cur$P, ctr = cur$ctr, type = cur$type),
       loop = if (base$closed) base else NULL, branches = branches, ribbons_sep = rem)
}

## ---- Beads ----
beads_from_path <- function(P, ctr, closed, w, crs) {
  hw <- w / 2
  line <- st_sfc(st_linestring(P), crs = crs)
  strip <- st_buffer(line, hw, endCapStyle = "FLAT", joinStyle = "ROUND", nQuadSegs = 10)
  cum <- path_cum(P); seeds <- interp_pt(P, cum, ctr)
  if (nrow(seeds) == 1) {
    cells <- strip
  } else {
    vor <- st_voronoi(st_sfc(st_multipoint(seeds), crs = crs), envelope = st_buffer(strip, 300)) |> st_collection_extract("POLYGON")
    vor <- st_sfc(vor, crs = crs)
    seed_pts <- st_sfc(lapply(seq_len(nrow(seeds)), \(i) st_point(seeds[i, ])), crs = crs)
    idx <- vapply(seq_len(nrow(seeds)), \(i) st_within(seed_pts[i], vor)[[1]][1], 1L)
    cells <- st_sfc(lapply(seq_len(nrow(seeds)), \(i) {
      pc <- suppressWarnings(st_intersection(vor[idx[i]], strip))
      pc <- st_cast(st_cast(pc, "MULTIPOLYGON"), "POLYGON")
      pc[which.max(as.numeric(st_area(pc)))][[1]] }), crs = crs)
  }
  st_sf(geometry = cells, area_ha = as.numeric(st_area(cells)) / 10000)
}

## ---- One reef's necklace ----
make_necklace <- function(sites, tows, dt, L = 500, area = target_area_m2, ds = 25) {
  crs <- st_crs(sites); w <- area / L
  base <- chain_path(sites, L, ds)
  chain_line <- st_sfc(st_linestring(base$P), crs = crs)
  rib <- tryCatch(tow_paths(tows, chain_line, dt, ds), error = \(e) list())
  jn <- join_pieces(base, rib, L, ds)
  pieces <- list()
  if (!is.null(jn$loop)) pieces[[1]] <- list(P = jn$loop$P, ctr = jn$loop$ctr, closed = TRUE, type = rep("eDNA", length(jn$loop$ctr)))
  if (!is.null(jn$path)) pieces[[length(pieces) + 1]] <- list(P = jn$path$P, ctr = jn$path$ctr, closed = FALSE, type = jn$path$type)
  for (b in jn$branches) pieces[[length(pieces) + 1]] <- list(P = b$P, ctr = b$ctr, closed = FALSE, type = b$type)
  for (R in jn$ribbons_sep) {
    len <- max(path_cum(R)); nb <- max(2, round(len / L))
    pieces[[length(pieces) + 1]] <- list(P = R, ctr = (seq_len(nb) - 0.5) * len / nb, closed = FALSE, type = rep("tow", nb))
  }
  out <- lapply(pieces, \(pc) { b <- beads_from_path(pc$P, pc$ctr, pc$closed, w, crs); b$source <- pc$type; b })
  cells <- do.call(rbind, out)
  # beads that touch count as one piece
  comp <- length(st_cast(st_union(st_buffer(st_geometry(cells), 1)), "POLYGON"))
  list(cells = cells, closed = base$closed, length_m = base$ell, n_pieces = comp)
}

## ---- Figure layout ----
## Arrange the four maps 1x4, 2x2 or 4x1, whichever gives the largest maps within the size box.
plan_layout <- function(asp, W = fig_w_max_cm, Hc = fig_h_max_cm, oh_w = 0.5, oh_h = 1.6, oh_top = 1.4) {
  o <- data.frame(ncol = c(1, 2, 4), nrow = c(4, 2, 1))
  o$mw_max <- (W - o$ncol * oh_w) / o$ncol
  o$mh_max <- (Hc - o$nrow * oh_h - oh_top) / o$nrow
  o$mw <- pmin(o$mw_max, pmax(o$mh_max, 0.1) * asp)
  o$mh <- o$mw / asp
  best <- o[which.max(o$mw * o$mh), ]
  list(ncol = best$ncol, nrow = best$nrow, mw = best$mw, mh = best$mh,
       w_cm = best$ncol * (best$mw + oh_w), h_cm = best$nrow * (best$mh + oh_h) + oh_top)
}

## Render, trim white margins, add a thin border and write an 8-bit RGB PNG.
save_fit <- function(plot, file, w_cm, h_cm, dpi = 250, pad_px = 12) {
  tmp <- tempfile(fileext = ".png"); on.exit(unlink(tmp))
  ggsave(tmp, plot, width = w_cm, height = h_cm, units = "cm", dpi = dpi, bg = "white")
  im <- image_border(image_trim(image_read(tmp), fuzz = 2), "white", sprintf("%dx%d", pad_px, pad_px))
  image_write(image_convert(im, type = "TrueColor", depth = 8), file, format = "png", depth = 8, quality = 95)
}

## ---- Per-reef map ----
make_necklace_reef <- function(reef_lbl, target_area_m2 = 100000, zoom = 15) {

  site_pts <- edna_mapped |>
    filter(reef_label == reef_lbl) |>
    group_by(site_no) |>
    summarise(lat = mean(lat), long = mean(long),
              mean_prop = mean(lod_sample_positive, na.rm = TRUE),
              mean_conc = mean(conc_mean, na.rm = TRUE),
              .groups = "drop")
  if (nrow(site_pts) < 2 || anyNA(site_pts$lat))
    return(list(ok = FALSE, reason = "insufficient eDNA sites"))

  tow_r <- tow_level |>
    filter(reef_label == reef_lbl, !is.na(tow_mid_lon), !is.na(tow_mid_lat))
  if (nrow(tow_r) < 3)
    return(list(ok = FALSE, reason = "insufficient manta tows"))

  # Local UTM zone so geometry is planar (m)
  zone    <- floor((mean(c(tow_r$tow_mid_lon, site_pts$long)) + 180) / 6) + 1
  utm_crs <- if (mean(tow_r$tow_mid_lat) < 0) 32700 + zone else 32600 + zone

  sites_utm <- st_as_sf(site_pts, coords = c("long", "lat"), crs = 4326) |>
    st_transform(utm_crs)
  tows_utm <- st_as_sf(tow_r, coords = c("tow_mid_lon", "tow_mid_lat"),
                       crs = 4326, remove = FALSE) |>
    st_transform(utm_crs)

  chain <- make_necklace(st_geometry(sites_utm), st_geometry(tows_utm), tow_r$tow_datetime)
  beads <- chain$cells |>
    mutate(bead_idx = row_number(),
           cull_site_artificial = paste0(reef_lbl, "_", source, "_", bead_idx))

  # Assign tows and eDNA sites to beads (1 m tolerance; first bead wins)
  first_hit <- \(m) apply(m, 1, \(x) if (any(x)) which(x)[1] else NA_integer_)
  tow_r$bead_idx       <- first_hit(st_is_within_distance(tows_utm,  beads, dist = 1, sparse = FALSE))
  site_pts$bead_idx    <- first_hit(st_is_within_distance(sites_utm, beads, dist = 1, sparse = FALSE))

  bead_tows <- tow_r |>
    filter(!is.na(bead_idx)) |>
    group_by(bead_idx) |>
    summarise(
      n_tows              = n(),
      manta_cots_rate     = sum(cots_observed, na.rm = TRUE) / n_tows,
      scars_per_tow       = mean(scar_count, na.rm = TRUE),
      manta_combined_rate = manta_cots_rate + scars_per_tow,
      n_tows_with_surf    = sum(!is.na(cots_reefscan_surface)),
      rs_surf_rate        = if_else(n_tows_with_surf > 0,
                                    sum(cots_reefscan_surface, na.rm = TRUE) / n_tows_with_surf,
                                    NA_real_),
      .groups = "drop"
    )

  # Sites sharing a bead are averaged
  bead_edna <- site_pts |>
    filter(!is.na(bead_idx)) |>
    group_by(bead_idx) |>
    summarise(n_edna_sites = n(), mean_prop = mean(mean_prop),
              mean_conc = mean(mean_conc), .groups = "drop")

  beads_full <- st_transform(beads, 4326) |>
    left_join(bead_tows, by = "bead_idx") |>
    left_join(bead_edna, by = "bead_idx") |>
    mutate(has_edna = !is.na(n_edna_sites))

  # Map extent: beads plus a thin margin
  bb <- st_bbox(beads_full)
  x0 <- bb[["xmin"]]; x1 <- bb[["xmax"]]; y0 <- bb[["ymin"]]; y1 <- bb[["ymax"]]
  pad_x <- max(0.03 * (x1 - x0), 0.0015); pad_y <- max(0.03 * (y1 - y0), 0.0015)
  x0 <- x0 - pad_x; x1 <- x1 + pad_x; y0 <- y0 - pad_y; y1 <- y1 + pad_y
  tile_bbox <- st_bbox(c(xmin = x0, xmax = x1, ymin = y0, ymax = y1), crs = 4326)
  span_x  <- x1 - x0
  asp_wh  <- (x1 - x0) * cos(mean(c(y0, y1)) * pi / 180) / (y1 - y0)   # width / height on the ground

  # Satellite tiles (cached): finer zoom for small extents, coarser for wide ones
  zoom_use <- if (span_x > 0.12) 14 else if (span_x < 0.025) 16 else zoom
  fetch <- \(z) tryCatch(get_tiles(st_as_sfc(tile_bbox), provider = "Esri.WorldImagery",
                                   zoom = z, crop = TRUE, cachedir = tile_dir), error = function(e) NULL)
  tiles <- fetch(zoom_use)
  if (is.null(tiles) && zoom_use != zoom) tiles <- fetch(zoom)
  if (is.null(tiles)) return(list(ok = FALSE, reason = "tile fetch failed"))

  lay <- plan_layout(asp_wh)
  narrow <- lay$mw < 5
  wt <- \(x) if (narrow) stringr::str_wrap(x, 14) else x
  mk_panel <- function(fill_var, legend_title, title_txt, trans = "identity", scale_bar = FALSE) {
    bar_w <- max(1.5, min(3, 0.8 * lay$mw))
    ggplot() +
      geom_spatraster_rgb(data = tiles) +
      geom_sf(data = beads_full, aes(fill = .data[[fill_var]]),
              color = "white", linewidth = 0.15, alpha = 0.75) +
      geom_sf(data = filter(beads_full, has_edna), fill = NA,
              color = "red", linewidth = 0.5) +
      coord_sf(xlim = c(tile_bbox["xmin"], tile_bbox["xmax"]),
               ylim = c(tile_bbox["ymin"], tile_bbox["ymax"]), expand = FALSE) +
      scale_fill_viridis_c(option = "plasma", na.value = NA, trans = trans,
                           guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                                  barwidth = unit(bar_w, "cm"), barheight = unit(0.2, "cm"))) +
      (if (scale_bar && requireNamespace("ggspatial", quietly = TRUE))
         ggspatial::annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.5, height = unit(0.1, "cm"),
                                     bar_cols = c("white", "black"), text_col = "white", line_col = "white") else NULL) +
      labs(title = wt(title_txt), fill = legend_title, x = NULL, y = NULL) +
      theme_void(base_size = 7) +
      theme(plot.title = element_text(size = 6.5, face = "bold", margin = margin(0, 0, 2, 0)),
            legend.position = "bottom", legend.title = element_text(size = 5.5),
            legend.text = element_text(size = 5), legend.margin = margin(1, 0, 0, 0),
            plot.margin = margin(1, 2, 1, 2),
            panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.3))
  }

  panels <- list(mk_panel("mean_conc", "eDNA conc.", "eDNA concentration", trans = "sqrt", scale_bar = TRUE),
                 mk_panel("mean_prop", "proportion positive", "eDNA proportion positive"),
                 mk_panel("manta_combined_rate", "COTS + scars per tow", "Manta density"),
                 mk_panel("rs_surf_rate", "COTS per tow", "ReefScan Surface density"))
  plot <- wrap_plots(panels, ncol = lay$ncol) +
    plot_annotation(
      title = paste0("Reef ", reef_lbl),
      caption = "Red outline: bead contains an eDNA site; unfilled: no data. Imagery: Esri World Imagery.",
      theme = theme(plot.title = element_text(size = 9, face = "bold"),
                    plot.caption = element_text(size = 5.5, colour = "grey30"),
                    plot.margin = margin(2, 2, 2, 2))
    )

  list(ok = TRUE, plot = plot,
       beads = st_drop_geometry(beads_full) |> mutate(reef_label = reef_lbl, .before = 1),
       n_beads = nrow(beads_full), n_sites = nrow(site_pts),
       n_sites_in_bead = sum(!is.na(site_pts$bead_idx)),
       n_tows = nrow(tow_r), n_tows_in_bead = sum(!is.na(tow_r$bead_idx)),
       median_area_ha = median(beads_full$area_ha),
       chain_closed = chain$closed, chain_km = chain$length_m / 1000,
       n_site_cells = sum(beads_full$source == "eDNA"), n_tow_cells = sum(beads_full$source == "tow"),
       n_link_cells = sum(beads_full$source == "link"), n_pieces = chain$n_pieces,
       fig_w = lay$w_cm, fig_h = lay$h_cm)
}

## ---- Run all reefs with eDNA and manta tows ----
reefs <- sort(intersect(unique(edna_mapped$reef_label), unique(tow_level$reef_label)))

# Each figure is saved as soon as it is built and the plot dropped, to limit memory use.
results <- map(set_names(reefs), \(rl) {
  res <- tryCatch(
    make_necklace_reef(rl, target_area_m2, tile_zoom),
    error = \(e) list(ok = FALSE, reason = paste("error:", conditionMessage(e)))
  )
  if (isTRUE(res$ok)) {
    save_fit(res$plot, file.path(fig_dir, paste0("necklace_", rl, ".png")), res$fig_w, res$fig_h)
    res$plot <- NULL
  }
  res
})

## ---- Tables ----
edna_subreef_beads <- results |> keep(\(r) isTRUE(r$ok)) |> map("beads") |> list_rbind()

edna_subreef_summary <- imap(results, \(res, rl) {
  if (!isTRUE(res$ok)) return(tibble(reef_label = rl, status = res$reason))
  tibble(reef_label = rl, status = "ok", n_beads = res$n_beads,
         n_edna_chain_cells = res$n_site_cells, n_tow_ribbon_cells = res$n_tow_cells,
         n_link_cells = res$n_link_cells, n_pieces = res$n_pieces,
         median_area_ha = round(res$median_area_ha, 2),
         chain_type = if (res$chain_closed) "loop" else "open", chain_km = round(res$chain_km, 1),
         n_edna_sites = res$n_sites, n_edna_sites_in_bead = res$n_sites_in_bead,
         n_tows = res$n_tows, n_tows_in_bead = res$n_tows_in_bead)
}) |> list_rbind()

saveRDS(edna_subreef_beads, "data/derived/edna_subreef_beads.RDS")
write_csv(edna_subreef_summary, "data/derived/edna_subreef_summary.csv")

## ---- Optional one-page summary figure ----
if (isTRUE(make_summary_figure)) {
## A: share of each reef's tows inside a bead.
## B: bead-level eDNA split by whether the imaging method detected COTS in the bead
##    (beads with an eDNA site only; descriptive, no test).
bead_eDNA <- edna_subreef_beads |>
  filter(has_edna) |>
  select(reef_label, mean_prop, mean_conc, manta = manta_combined_rate, `RS Surface` = rs_surf_rate) |>
  pivot_longer(c(manta, `RS Surface`), names_to = "method", values_to = "density") |>
  filter(!is.na(density)) |>
  mutate(method = recode(method, manta = "Manta (COTS + scars)"),
         detected = factor(if_else(density > 0, "detected", "not detected"),
                           levels = c("not detected", "detected"))) |>
  pivot_longer(c(mean_prop, mean_conc), names_to = "metric", values_to = "edna") |>
  mutate(metric = recode(metric, mean_prop = "eDNA proportion positive",
                         mean_conc = "eDNA concentration"))

n_lab <- bead_eDNA |>
  filter(metric == "eDNA proportion positive") |>
  count(method, detected) |>
  mutate(xlab = paste0(detected, "
(n = ", n, ")"))
bead_eDNA <- bead_eDNA |>
  left_join(select(n_lab, method, detected, xlab), by = c("method", "detected")) |>
  mutate(xlab = fct_reorder(xlab, as.integer(detected)))

p_cov <- edna_subreef_summary |>
  mutate(pct_tows = 100 * n_tows_in_bead / n_tows,
         reef_label = fct_reorder(reef_label, pct_tows)) |>
  ggplot(aes(x = pct_tows, y = reef_label)) +
  geom_segment(aes(x = 0, xend = pct_tows, yend = reef_label), colour = "grey70") +
  geom_point() +
  geom_text(aes(x = 108, label = n_beads), size = 2.4, hjust = 1) +
  annotate("text", x = 108, y = nrow(edna_subreef_summary) + 1.2, label = "beads",
           size = 2.4, hjust = 1, vjust = 0) +
  scale_x_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  coord_cartesian(clip = "off") +
  labs(title = "A. Tows captured by beads", x = "% of manta tows inside a bead", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(plot.margin = margin(14, 4, 4, 4), plot.title = element_text(size = 8.5, face = "bold"))

p_det <- ggplot(bead_eDNA, aes(x = xlab, y = edna)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6, size = 1.1) +
  # median tick
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "errorbar",
               width = 0.5, colour = "red", linewidth = 0.6) +
  facet_grid(metric ~ method, scales = "free", switch = "y") +
  labs(title = "B. eDNA by COTS detection in the same bead",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
        plot.title = element_text(size = 8.5, face = "bold"))

p_summary <- (p_cov | p_det) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(
    title = sprintf("Sub-reef eDNA: %d beads of ~10 ha (median %.1f ha) across %d reefs",
                    nrow(edna_subreef_beads), median(edna_subreef_beads$area_ha),
                    nrow(edna_subreef_summary)),
    subtitle = sprintf("All %d eDNA sites sit in a bead; %d of %d tows are in a bead. Panel B shows only beads containing an eDNA site.",
                       sum(edna_subreef_summary$n_edna_sites_in_bead),
                       sum(edna_subreef_summary$n_tows_in_bead), sum(edna_subreef_summary$n_tows)),
    caption = "Red tick = median. Beads are clustered within reefs and COTS detections are rare per bead; descriptive only.",
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7.5), plot.caption = element_text(size = 6.5))
  )

ggsave("output/figures/edna_necklace_summary.png", p_summary,
       width = fig_w_max_cm, height = 15, units = "cm", dpi = 300, bg = "white")
}

edna_subreef_summary
