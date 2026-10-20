# process_plots.R -- conversion of the sampling plots and plot areas.
#
# Inputs : plot table (kplot, character columns), tower positions, eco label.
# Outputs: list(update, full, points, areas, counts)
#   update  ID, PLOT_LOCATION_LAT, PLOT_LOCATION_LONG for plots submitted as N/E or polar
#   full    every plot with format, reference, declination, N/E, polar, lat/long
#   points  data.frame(plot_id, type, lat, lon) of resolved plots (also used as
#           references by the BADM groups)
#   areas   sf data.frame(ID, type, geometry) of SP-I and CP areas
# Hierarchy: TOWER (GRP_LOCATION valid at the date) -> CP, SP-I ; SP-I -> SP-II.

#' Process the sampling-plot table of a station.
#'
#' @param plot_df plot table from fetch_plot_table().
#' @param towers output of tower_positions(); if empty, the plot TOWER row is used.
#' @param eco eco label (forest, crop, grassland, mire, urban, water or NA).
#' @param site station code.
#' @return list(update, full, points, areas, counts, towers).
process_plots <- function(plot_df, towers, eco, site = NA) {
    spec <- GROUP_SPECS$GRP_PLOT
    g <- "GRP_PLOT"
    if (is.null(plot_df) || !nrow(plot_df) || !"PLOT_ID" %in% names(plot_df)) {
        log_warn(g, ": plot table empty or without PLOT_ID")
        return(NULL)
    }
    type <- plot_type(plot_df$PLOT_ID)
    if (any(is.na(type)))
        log_info(g, ": ", sum(is.na(type)), " row(s) with an unrecognised PLOT_ID ignored")
    df <- plot_df[!is.na(type), , drop = FALSE]
    type <- type[!is.na(type)]

    rec <- .rec_template(nrow(df))
    rec$id <- col_or_na(df, "ID")
    rec$date <- parse_badm_date(coalesce_cols(df, spec$date_cols))
    rec$northref <- coalesce_cols(df, spec$northref_cols)
    ne <- .pick_pair(df, spec$ne); pol <- .pick_pair(df, spec$polar); ll <- .pick_pair(df, list(spec$latlong))
    rec$in_east <- ne$a; rec$in_north <- ne$b; rec$in_distance <- pol$a; rec$in_azimuth <- pol$b
    rec$in_lat <- ll$a; rec$in_lon <- ll$b
    plot_id <- trimws(df$PLOT_ID)

    # Tower: GRP_LOCATION by date; the plot TOWER row is only a fallback / check
    tw_rows <- which(type == "TOWER")
    if (!nrow(towers)) {
        if (length(tw_rows) && any(!is.na(rec$in_lat[tw_rows]))) {
            log_warn(g, ": GRP_LOCATION has no position; using the TOWER row of the plot table")
            towers <- data.frame(location_id = paste0("PLOT", rec$id[tw_rows]), date = rec$date[tw_rows],
                                 lat = rec$in_lat[tw_rows], lon = rec$in_lon[tw_rows])
            towers <- towers[!is.na(towers$lat), , drop = FALSE]
            towers <- towers[order(!is.na(towers$date), towers$date), , drop = FALSE]
        } else {
            stop("no tower position (GRP_LOCATION empty and no TOWER row in the plot table)")
        }
    }

    # TOWER rows, CP and SP-I: reference = tower at the record date
    lvl1 <- type %in% c("TOWER", "CP", "SP-I")
    tw <- tower_at(towers, rec$date[lvl1])
    rec$ref_id[lvl1] <- tw$ref_id; rec$ref_lat[lvl1] <- tw$ref_lat; rec$ref_lon[lvl1] <- tw$ref_lon
    rec$ref_rule[lvl1] <- "tower"
    nd <- sum(is.na(rec$date[lvl1]) & type[lvl1] != "TOWER" &
              (!is.na(rec$in_east[lvl1]) | !is.na(rec$in_distance[lvl1])))
    if (nd) log_warn(g, ": ", nd, " plot(s) without a date use the latest tower position")
    r1 <- compute_records(rec[lvl1, , drop = FALSE], g)

    # Check: plot TOWER row vs GRP_LOCATION (> 1 m -> WARN)
    tr <- which(type[lvl1] == "TOWER")
    if (length(tr)) {
        d <- sqrt(r1$east[tr]^2 + r1$north[tr]^2)
        if (any(d > 1, na.rm = TRUE))
            log_warn(g, ": plot TOWER row differs from GRP_LOCATION by up to ", round(max(d, na.rm = TRUE), 2), " m")
    }

    # SP-II: reference = resolved parent SP-I
    lvl2 <- type == "SP-II"
    spi <- data.frame(plot_id = plot_id[lvl1][type[lvl1] == "SP-I"], lat = r1$lat[type[lvl1] == "SP-I"],
                      lon = r1$lon[type[lvl1] == "SP-I"], stringsAsFactors = FALSE)
    r2 <- NULL
    if (any(lvl2)) {
        parent <- sp2_parent(plot_id[lvl2], spi$plot_id)
        m <- match(parent, spi$plot_id)
        rec$ref_id[lvl2] <- ifelse(is.na(parent), NA, parent)
        rec$ref_lat[lvl2] <- spi$lat[m]; rec$ref_lon[lvl2] <- spi$lon[m]
        rec$ref_rule[lvl2] <- "parent SP-I"
        if (any(is.na(m))) log_warn(g, ": ", sum(is.na(m)), " SP-II plot(s) without a resolved parent SP-I")
        r2 <- compute_records(rec[lvl2, , drop = FALSE], g)
    }
    res <- rbind(r1, r2)
    pid <- c(plot_id[lvl1], plot_id[lvl2])
    ptype <- c(type[lvl1], type[lvl2])
    src_rows <- c(which(lvl1), which(lvl2))

    full <- .full_table(res, g, data.frame(PLOT_ID = pid, plot_type = ptype, stringsAsFactors = FALSE))
    full <- full[, c("ID", "PLOT_ID", "plot_type", setdiff(names(full), c("ID", "PLOT_ID", "plot_type")))]
    upd <- res$status == "ok" & res$input_format %in% c("ne", "polar") & ptype != "TOWER"
    update <- data.frame(ID = res$id[upd], PLOT_LOCATION_LAT = res$lat[upd], PLOT_LOCATION_LONG = res$lon[upd],
                         stringsAsFactors = FALSE)
    if (!any(res$input_format %in% c("ne", "polar") & ptype != "TOWER"))
        log_info(g, ": no plot update file: all points submitted in lat/long")

    ok <- !is.na(res$lat) & ptype != "TOWER"
    points <- data.frame(plot_id = pid[ok], type = ptype[ok], lat = res$lat[ok], lon = res$lon[ok],
                         stringsAsFactors = FALSE)
    latest <- tower_at(towers, as.Date(NA))
    points <- rbind(data.frame(plot_id = "TOWER", type = "TOWER", lat = latest$ref_lat, lon = latest$ref_lon), points)

    areas <- plot_area_polygons(points, df[src_rows[ok], , drop = FALSE], res[ok, , drop = FALSE], eco, g)
    counts <- .counts(site, g, res, nrow(update))
    log_info(g, ": ", nrow(res), " plots (", sum(ptype == "CP"), " CP, ", sum(ptype == "SP-I"), " SP-I, ",
             sum(ptype == "SP-II"), " SP-II), ", nrow(update), " with calculated coordinates")
    list(update = update, full = full, points = points, areas = areas, counts = counts, towers = towers)
}

#' SP-I and CP area polygons.
#'
#' SP-I: circle of 15 m (forest) or 10 m (other). CP: rectangle when
#' PLOT_RECT_X/Y/DIR are complete (DIR corrected for declination if the plot's
#' north reference is magnetic); otherwise forest circle of 2000 m^2, crop circle
#' with the half-diagonal of a 5 x 10 m plot, mire 0.6 m north-aligned square;
#' no CP area for grassland, urban, water or unknown eco types.
#'
#' @param points resolved points (first row is the tower).
#' @param df source plot rows matching points[-1, ].
#' @param res processed records matching points[-1, ].
#' @return sf data.frame(ID, type, geometry) in WGS84, or NULL.
plot_area_polygons <- function(points, df, res, eco, g = "GRP_PLOT") {
    p <- points[-1, , drop = FALSE]
    if (!nrow(p)) return(NULL)
    geoms <- list(); ids <- character(0); types <- character(0)
    spi_r <- if (identical(eco, "forest")) PLOT_AREA_RULES$spi_radius[["forest"]] else PLOT_AREA_RULES$spi_radius[["other"]]
    for (i in which(p$type == "SP-I")) {
        geoms[[length(geoms) + 1]] <- polygon_from_ne(p$lat[i], p$lon[i], circle_ne(spi_r))
        ids <- c(ids, p$plot_id[i]); types <- c(types, "SP-I area")
    }
    cp <- which(p$type == "CP")
    if (length(cp)) {
        rx <- as_num(col_or_na(df, "PLOT_RECT_X")); ry <- as_num(col_or_na(df, "PLOT_RECT_Y"))
        rd <- as_num(col_or_na(df, "PLOT_RECT_DIR"))
        no_area <- 0
        for (i in cp) {
            ring <- NULL
            if (!is.na(rx[i]) && !is.na(ry[i]) && !is.na(rd[i])) {
                dir <- rd[i]
                if (is_magnetic(res$northref[i])) {
                    if (is.na(res$date[i])) {
                        log_warn(g, ": rectangular ", p$plot_id[i], " has a magnetic direction but no date; area skipped")
                        next
                    }
                    dir <- correct_magnetic_azimuth(dir, declination(p$lat[i], p$lon[i], res$date[i]))
                }
                ring <- rect_corners_ne(min(rx[i], ry[i]), max(rx[i], ry[i]), dir)
            } else if (identical(eco, "forest")) {
                ring <- circle_ne(PLOT_AREA_RULES$cp_forest_r)
            } else if (identical(eco, "crop")) {
                ring <- circle_ne(PLOT_AREA_RULES$cp_crop_r)
            } else if (identical(eco, "mire")) {
                ring <- square_ne(PLOT_AREA_RULES$cp_mire_side)
            }
            if (is.null(ring)) { no_area <- no_area + 1; next }
            geoms[[length(geoms) + 1]] <- polygon_from_ne(p$lat[i], p$lon[i], ring)
            ids <- c(ids, p$plot_id[i]); types <- c(types, "CP area")
        }
        if (no_area) log_info(g, ": no CP area defined for eco type '", eco, "' (", no_area, " CP)")
    }
    if (!length(geoms)) return(NULL)
    sf::st_sf(ID = ids, type = types, geometry = sf::st_sfc(geoms, crs = 4326))
}
