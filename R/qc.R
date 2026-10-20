# qc.R -- quality-control flags for human review.
#
# Records are never withheld from the update files because of a flag: every
# record that can be computed gets its value, and problems are flagged.
# Flags (one row per record and flag: site, group, ID, flag, detail):
#   FORMAT_DISAGREEMENT  N/E and polar both submitted and their positions differ by more
#                        than qc$format_disagreement_max_m (plot lat/long values are earlier
#                        uploads, not user input, and are not compared here)
#   OUTSIDE_TA           computed point outside the station's Target Area polygon (only when
#                        the TA is available); detail = distance to the TA boundary
#   LARGE_SHIFT          new position differs from the previous value by more than
#                        qc$large_shift_m; previous value = the calculated coordinates
#                        returned by the API, else the --previous merged files
#   NOT_COMPUTED         record with N/E or polar input whose position could not be
#                        computed (e.g. unknown reference plot, magnetic record without a
#                        date); detail = reason. These are the only records absent from
#                        the update files.
#   GAI_METHOD_UNKNOWN   GAI record with coordinates whose GAI_METHOD is empty or not in
#                        GAI_METHOD_PAIRS; the first complete pair (DHP, then CEPT) is used
#   GAI_PAIR_MISMATCH    the pair matching GAI_METHOD is empty and the other pair is used

QC_DEFAULTS <- list(format_disagreement_max_m = 1.0, large_shift_m = 20)
QC_FLAGS <- c("FORMAT_DISAGREEMENT", "OUTSIDE_TA", "LARGE_SHIFT", "NOT_COMPUTED", "GAI_METHOD_UNKNOWN", "GAI_PAIR_MISMATCH")

#' QC thresholds from settings, completed with QC_DEFAULTS.
qc_settings <- function(settings) {
    q <- QC_DEFAULTS
    for (k in names(settings$qc)) q[[k]] <- settings$qc[[k]]
    q
}

.flag_df <- function(site, group, id, flag, detail) {
    if (!length(id)) return(NULL)
    data.frame(site = site, group = group, ID = id, flag = flag, detail = detail, stringsAsFactors = FALSE)
}

#' Distances (m) between point pairs near a station, in the AEQD plane of the tower.
#'
#' Accurate to well below 1 mm for points within a few km of the tower.
local_distance_m <- function(lat1, lon1, lat2, lon2, tower_lat, tower_lon) {
    a <- latlon_to_ne(tower_lat, tower_lon, lat1, lon1)
    b <- latlon_to_ne(tower_lat, tower_lon, lat2, lon2)
    sqrt((a$east - b$east)^2 + (a$north - b$north)^2)
}

#' Previous values from a directory of merged files (GRP_<X>.tsv with ID and calculated columns).
#'
#' @param dir directory or NULL.
#' @return named list group -> data.frame(ID, lat, lon), or NULL.
read_previous_values <- function(dir) {
    if (is.null(dir)) return(NULL)
    if (!dir.exists(dir)) stop("--previous directory not found: ", dir, call. = FALSE)
    out <- list()
    for (g in c("GRP_PLOT", CONVERTED_GROUPS)) {
        p <- file.path(dir, paste0(g, ".tsv"))
        if (!file.exists(p)) next
        x <- read_api_tsv(p, paste("previous", g))
        sp <- GROUP_SPECS[[g]]
        if (!nrow(x) || !all(c("ID", sp$out_lat, sp$out_lon) %in% names(x))) next
        x <- x[!duplicated(x$ID), , drop = FALSE]
        out[[g]] <- data.frame(ID = x$ID, lat = as_num(x[[sp$out_lat]]), lon = as_num(x[[sp$out_lon]]))
    }
    out
}

#' QC flags of one processed group.
#'
#' @param site station code; @param group group name.
#' @param full full table of the group (process_group()/process_plots()).
#' @param qc thresholds (qc_settings()).
#' @param tower list(lat, lon) for local distances.
#' @param ta sf polygon(s) of the Target Area, or NULL.
#' @param previous data.frame(ID, lat, lon) from --previous, or NULL.
#' @return data.frame of flags (may have 0 rows).
qc_group_flags <- function(site, group, full, qc, tower, ta = NULL, previous = NULL) {
    if (is.null(full) || !nrow(full)) return(NULL)
    spec <- GROUP_SPECS[[group]]
    fl <- list()
    ok <- full$status == "ok" & full$input_format %in% c("ne", "polar")
    if (group == "GRP_PLOT") ok <- ok & full$plot_type != "TOWER"
    u <- full[ok, , drop = FALSE]

    # FORMAT_DISAGREEMENT
    d <- u$ne_polar_disagreement_m
    k <- !is.na(d) & d > qc$format_disagreement_max_m
    fl$fd <- .flag_df(site, group, u$ID[k], "FORMAT_DISAGREEMENT",
                      sprintf("disagreement_m=%.2f; format_used=%s", d[k], u$input_format[k]))

    # NOT_COMPUTED
    nc <- grepl("^skipped", full$status)
    fl$nc <- .flag_df(site, group, full$ID[nc], "NOT_COMPUTED", full$status[nc])

    # Instrument-pair flags set during processing (GAI): GAI_METHOD_UNKNOWN, GAI_PAIR_MISMATCH
    if (!is.null(full$pair_flag)) {
        k <- !is.na(full$pair_flag)
        fl$pf <- .flag_df(site, group, full$ID[k], full$pair_flag[k], full$pair_detail[k])
    }

    if (nrow(u)) {
        # LARGE_SHIFT: previous = API calculated values, else --previous files
        if (group == "GRP_PLOT") {
            plat <- u$in_lat; plon <- u$in_lon
        } else {
            plat <- u[[paste0("db_", spec$out_lat)]]; plon <- u[[paste0("db_", spec$out_lon)]]
            if (is.null(plat)) plat <- plon <- rep(NA_real_, nrow(u))
        }
        psrc <- ifelse(!is.na(plat) & !is.na(plon), "api", NA)
        if (!is.null(previous)) {
            m <- match(u$ID, previous$ID)
            fill <- is.na(psrc) & !is.na(m)
            plat[fill] <- previous$lat[m[fill]]; plon[fill] <- previous$lon[m[fill]]; psrc[fill] <- "previous_dir"
        }
        has <- !is.na(psrc)
        if (any(has)) {
            sh <- rep(NA_real_, nrow(u))
            sh[has] <- local_distance_m(u$calc_lat[has], u$calc_lon[has], plat[has], plon[has], tower$lat, tower$lon)
            k <- !is.na(sh) & sh > qc$large_shift_m
            fl$ls <- .flag_df(site, group, u$ID[k], "LARGE_SHIFT",
                              sprintf("shift_m=%.2f; D_applied=%.3f; previous=%s", sh[k], u$magdec_applied[k], psrc[k]))
        }

        # OUTSIDE_TA
        if (!is.null(ta) && nrow(ta)) {
            crs <- sf::st_crs(aeqd_crs(tower$lat, tower$lon))
            poly <- sf::st_union(sf::st_make_valid(sf::st_transform(sf::st_geometry(ta), crs)))
            pts <- sf::st_transform(sf::st_as_sf(data.frame(x = u$calc_lon, y = u$calc_lat), coords = c("x", "y"), crs = 4326), crs)
            inside <- lengths(sf::st_intersects(pts, poly)) > 0
            if (any(!inside)) {
                dist <- as.numeric(sf::st_distance(pts[!inside, ], poly))
                fl$ta <- .flag_df(site, group, u$ID[!inside], "OUTSIDE_TA", sprintf("distance_to_TA_m=%.1f", dist))
            }
        }
    }
    out <- do.call(rbind, fl)
    if (is.null(out)) .flag_df(site, group, character(0), character(0), character(0)) else out
}

#' QC flags of a station (all processed groups).
#'
#' @param results named list group -> list(full, ...).
#' @param previous named list group -> data.frame(ID, lat, lon) (read_previous_values()).
qc_station_flags <- function(site, results, qc, tower, ta = NULL, previous = NULL) {
    fl <- lapply(names(results), function(g)
        qc_group_flags(site, g, results[[g]]$full, qc, tower, ta, previous[[g]]))
    out <- do.call(rbind, fl)
    if (is.null(out) || !nrow(out)) return(data.frame(site = character(0), group = character(0), ID = character(0),
                                                      flag = character(0), detail = character(0)))
    out
}
