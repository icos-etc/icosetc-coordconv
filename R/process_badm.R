# process_badm.R -- conversion of the BADM groups with spatial information.
#
# Inputs : BADM tables (character columns), tower positions, resolved plot points.
# Outputs: per group a list(update, full, counts):
#   update  ID [+ extra] + <X>_LOCATION_CALC_LONG/LAT for records submitted as N/E or polar
#   full    every record with input format, reference, declination, N/E, polar, lat/long
#   counts  one summary row
# compute_records() is the shared engine also used by process_plots.R.

#' Standard empty record table columns.
.rec_template <- function(n) {
    data.frame(id = rep(NA_character_, n), date = as.Date(rep(NA, n)), ref_id = NA_character_,
               ref_rule = NA_character_, ref_lat = NA_real_, ref_lon = NA_real_, northref = NA_character_,
               in_east = NA_real_, in_north = NA_real_, in_distance = NA_real_, in_azimuth = NA_real_,
               in_lat = NA_real_, in_lon = NA_real_, stringsAsFactors = FALSE)
}

#' Declination, correction and positions for a set of records.
#'
#' Records with a magnetic north reference and N/E or polar input get the
#' declination at their reference point and date. Magnetic records without a
#' date, and N/E or polar records without a reference point, are skipped.
#'
#' @param rec data.frame from .rec_template() (filled).
#' @param what label used in log messages (e.g. "GRP_TREE").
#' @return rec with positions (see compute_positions()) plus magdec_applied and status.
compute_records <- function(rec, what) {
    if (!nrow(rec)) return(cbind(rec, input_format = character(0), status = character(0)))
    has_pol <- !is.na(rec$in_distance) & !is.na(rec$in_azimuth)
    has_ne  <- !is.na(rec$in_east) & !is.na(rec$in_north)
    mag <- is_magnetic(rec$northref) & (has_pol | has_ne)
    no_date_mag <- mag & is.na(rec$date)
    rec$magdec <- 0
    need <- mag & !no_date_mag & !is.na(rec$ref_lat) & !is.na(rec$ref_lon)
    if (any(need)) {
        rec$magdec[need] <- withCallingHandlers(
            declination(rec$ref_lat[need], rec$ref_lon[need], rec$date[need]),
            warning = function(w) { log_warn(what, ": ", conditionMessage(w)); invokeRestart("muffleWarning") })
    }
    rec <- compute_positions(rec)
    fmt <- rec$input_format
    status <- rep("ok", nrow(rec))
    status[is.na(fmt)] <- "no coordinates"
    no_ref <- fmt %in% c("ne", "polar") & (is.na(rec$ref_lat) | is.na(rec$ref_lon))
    status[no_ref] <- "skipped: no reference point"
    status[no_date_mag & fmt %in% c("ne", "polar")] <- "skipped: magnetic record without date"
    bad <- status != "ok" & status != "no coordinates"
    for (cl in c("east", "north", "distance", "azimuth", "lat", "lon", "format_disagreement_m", "ne_polar_disagreement_m"))
        rec[[cl]][bad] <- NA
    rec$magdec_applied <- ifelse(bad, NA, ifelse(fmt %in% c("ne", "polar") & mag, rec$magdec, 0))
    rec$status <- status
    if (sum(no_date_mag & fmt %in% c("ne", "polar")))
        log_warn(what, ": ", sum(no_date_mag & fmt %in% c("ne", "polar")),
                 " magnetic record(s) without a date skipped (declination cannot be computed)")
    if (sum(rec$n_formats > 1, na.rm = TRUE))
        log_info(what, ": ", sum(rec$n_formats > 1, na.rm = TRUE),
                 " record(s) with more than one coordinate format (precedence polar > N/E > lat/long)")
    rec
}

# First complete (a, b) numeric pair over a list of column pairs.
.pick_pair <- function(df, pairs) {
    n <- nrow(df)
    a <- b <- rep(NA_real_, n)
    used <- rep(NA_integer_, n)                      # index of the pair used
    for (k in seq_along(pairs)) {
        p <- pairs[[k]]
        x <- as_num(col_or_na(df, p[1])); y <- as_num(col_or_na(df, p[2]))
        take <- is.na(a) & is.na(b) & !is.na(x) & !is.na(y)
        a[take] <- x[take]; b[take] <- y[take]; used[take] <- k
    }
    list(a = a, b = b, used = used)
}

#' Instrument-aware choice of the GAI position pair.
#'
#' The pair is chosen from the method column (GAI_METHOD_PAIRS); within the pair the
#' N/E and polar values are passed on and compute_positions() applies polar > N/E.
#' Known method but empty pair while the other pair is filled -> the other pair is used
#' (flag GAI_PAIR_MISMATCH). Unknown/empty method with coordinates -> first complete pair
#' in spec order (flag GAI_METHOD_UNKNOWN).
#'
#' @return list(ne, pol (each list(a, b)), method, pair, same_pair, flag, detail).
.select_instrument_pairs <- function(df, spec) {
    n <- nrow(df)
    P <- spec$instrument_pairs
    method <- trimws(col_or_na(df, spec$method_col))
    inst <- unname(GAI_METHOD_PAIRS[toupper(method)])
    inst[is.na(inst)] <- "unknown"
    val <- lapply(P, function(p) list(ne = .pick_pair(df, list(p$ne)), pol = .pick_pair(df, list(p$polar))))
    has <- vapply(val, function(v) !is.na(v$ne$a) | !is.na(v$pol$a), logical(n))
    if (n == 1) has <- matrix(has, nrow = 1, dimnames = list(NULL, names(P)))
    other <- setNames(rev(names(P)), names(P))
    use <- rep(NA_character_, n); flag <- detail <- rep(NA_character_, n)
    known <- inst %in% names(P)
    own <- known & has[cbind(seq_len(n), match(inst, names(P)))] %in% TRUE
    use[own] <- inst[own]
    alt <- known & !own & has[cbind(seq_len(n), match(other[inst], names(P)))] %in% TRUE
    use[alt] <- other[inst[alt]]
    flag[alt] <- "GAI_PAIR_MISMATCH"
    detail[alt] <- sprintf("method=%s; pair_used=%s", method[alt], use[alt])
    ne_a <- ne_b <- pol_a <- pol_b <- rep(NA_real_, n)
    for (k in names(P)) {
        r <- which(use == k)
        ne_a[r] <- val[[k]]$ne$a[r]; ne_b[r] <- val[[k]]$ne$b[r]
        pol_a[r] <- val[[k]]$pol$a[r]; pol_b[r] <- val[[k]]$pol$b[r]
    }
    same_pair <- rep(TRUE, n)
    fb <- !known & apply(has, 1, any)
    if (any(fb)) {
        ne <- .pick_pair(df, spec$ne); pol <- .pick_pair(df, spec$polar)
        ne_a[fb] <- ne$a[fb]; ne_b[fb] <- ne$b[fb]; pol_a[fb] <- pol$a[fb]; pol_b[fb] <- pol$b[fb]
        same_pair[fb] <- is.na(ne$used[fb]) | is.na(pol$used[fb]) | ne$used[fb] == pol$used[fb]
        u <- ifelse(!is.na(pol$used[fb]), pol$used[fb], ne$used[fb])
        use[fb] <- names(P)[u]
        flag[fb] <- "GAI_METHOD_UNKNOWN"
        detail[fb] <- sprintf("method=%s; fallback pair_used=%s", ifelse(is.na(method[fb]) | !nzchar(method[fb]), "<empty>", method[fb]), use[fb])
    }
    list(ne = list(a = ne_a, b = ne_b), pol = list(a = pol_a, b = pol_b), method = method, pair = use,
         same_pair = same_pair, flag = flag, detail = detail)
}

#' Build the full output table from processed records.
.full_table <- function(rec, group, extra = NULL) {
    out <- data.frame(
        ID = rec$id, group = group, date_used = format(rec$date), input_format = rec$input_format,
        n_formats = rec$n_formats, ref_id = rec$ref_id, ref_rule = rec$ref_rule,
        ref_lat = rec$ref_lat, ref_lon = rec$ref_lon, north_ref = rec$northref,
        magdec_applied = rec$magdec_applied,
        in_east_dist = rec$in_east, in_north_dist = rec$in_north,
        in_distance_polar = rec$in_distance, in_angle_polar = rec$in_azimuth,
        in_lat = rec$in_lat, in_lon = rec$in_lon,
        east_dist = rec$east, north_dist = rec$north, distance_polar = rec$distance, angle_polar = rec$azimuth,
        calc_lat = rec$lat, calc_lon = rec$lon, format_disagreement_m = rec$format_disagreement_m,
        ne_polar_disagreement_m = rec$ne_polar_disagreement_m,
        status = rec$status, stringsAsFactors = FALSE)
    if (!is.null(extra)) out <- cbind(out, extra)
    out
}

#' Summary counts of a processed group.
.counts <- function(site, group, rec, n_update) {
    data.frame(site = site, group = group, n_records = nrow(rec),
               n_ne = sum(rec$input_format %in% "ne"), n_polar = sum(rec$input_format %in% "polar"),
               n_latlong = sum(rec$input_format %in% "latlong"),
               n_magnetic = sum(rec$magdec_applied != 0, na.rm = TRUE),
               n_update = n_update, n_skipped = sum(grepl("^skipped", rec$status)),
               status = "ok", message = "", stringsAsFactors = FALSE)
}

#' Process one BADM group.
#'
#' @param df group table (character columns).
#' @param spec entry of GROUP_SPECS.
#' @param towers output of tower_positions().
#' @param plot_points data.frame(plot_id, lat, lon) of resolved plots.
#' @param site station code (for the summary).
#' @return list(update, full, counts) or NULL if the table is empty.
process_group <- function(df, spec, towers, plot_points, site = NA) {
    g <- spec$group
    if (is.null(df) || !nrow(df)) return(NULL)
    id <- col_or_na(df, "ID")
    noid <- is.na(id) | !nzchar(trimws(id))
    if (any(noid)) {
        log_warn(g, ": ", sum(noid), " record(s) without ID dropped")
        df <- df[!noid, , drop = FALSE]; id <- id[!noid]
    }
    n <- nrow(df)
    rec <- .rec_template(n)
    rec$id <- id
    rec$date <- parse_badm_date(coalesce_cols(df, spec$date_cols))
    rec$northref <- coalesce_cols(df, spec$northref_cols)
    sel <- NULL
    if (!is.null(spec$method_col)) {
        sel <- .select_instrument_pairs(df, spec)
        ne <- sel$ne; pol <- sel$pol
    } else {
        ne <- .pick_pair(df, spec$ne); pol <- .pick_pair(df, spec$polar)
    }
    rec$in_east <- ne$a; rec$in_north <- ne$b
    rec$in_distance <- pol$a; rec$in_azimuth <- pol$b
    if (length(spec$latlong)) {
        ll <- .pick_pair(df, list(spec$latlong)); rec$in_lat <- ll$a; rec$in_lon <- ll$b
    }
    has_xy <- !is.na(rec$in_east) | !is.na(rec$in_distance)
    plot_codes <- if (is.null(spec$plot_col)) NULL else col_or_na(df, spec$plot_col)
    ref <- resolve_record_reference(plot_codes, rec$date, plot_points, towers, spec$reference)
    rec$ref_id <- ref$ref_id; rec$ref_lat <- ref$ref_lat; rec$ref_lon <- ref$ref_lon; rec$ref_rule <- ref$ref_rule

    nf <- sum(ref$ref_rule == "tower_fallback" & has_xy)
    if (nf) log_info(g, ": ", nf, " record(s) without a plot reference (empty/outside/EC) referenced to the tower")
    nd <- sum(is.na(rec$date) & has_xy & ref$ref_rule %in% c("tower", "tower_fallback"))
    if (nd) log_warn(g, ": ", nd, " record(s) without a date use the latest tower position")
    unk <- ref$ref_rule == "unknown_plot" & has_xy
    if (any(unk)) log_warn(g, ": ", sum(unk), " record(s) skipped, unknown reference plot(s): ",
                           paste(head(unique(ref$ref_id[unk]), 5), collapse = ", "))

    rec <- compute_records(rec, g)
    rec$status[ref$ref_rule == "unknown_plot" & grepl("^skipped: no reference", rec$status)] <- "skipped: unknown reference plot"

    extra <- NULL
    if (!is.null(sel)) {
        # formats are compared only within one instrument pair
        rec$ne_polar_disagreement_m[!sel$same_pair] <- NA
        extra <- data.frame(gai_method = sel$method, instrument_pair = sel$pair,
                            pair_flag = sel$flag, pair_detail = sel$detail, stringsAsFactors = FALSE)
        for (fl in c("GAI_METHOD_UNKNOWN", "GAI_PAIR_MISMATCH"))
            if (sum(sel$flag %in% fl)) log_info(g, ": ", sum(sel$flag %in% fl), " record(s) with ", fl)
    }
    if (length(spec$db_cols)) {
        dbx <- as.data.frame(lapply(setNames(spec$db_cols, paste0("db_", spec$db_cols)),
                                    function(cl) as_num(col_or_na(df, cl))), check.names = FALSE)
        dlat <- dbx[[paste0("db_", spec$out_lat)]]; dlon <- dbx[[paste0("db_", spec$out_lon)]]
        ok <- !is.na(dlat) & !is.na(dlon) & !is.na(rec$lat)
        dbx$diff_vs_db_m <- NA_real_
        if (any(ok)) dbx$diff_vs_db_m[ok] <- distance_m(rec$lat[ok], rec$lon[ok], dlat[ok], dlon[ok])
        extra <- if (is.null(extra)) dbx else cbind(extra, dbx)
    }
    full <- .full_table(rec, g, extra)

    upd_rows <- rec$status == "ok" & rec$input_format %in% c("ne", "polar")
    update <- data.frame(ID = rec$id[upd_rows], stringsAsFactors = FALSE)
    if ("EC_SA_NORTH_MAGDEC" %in% spec$extra_out) update$EC_SA_NORTH_MAGDEC <- rec$magdec_applied[upd_rows]
    update[[spec$out_lon]] <- rec$lon[upd_rows]
    update[[spec$out_lat]] <- rec$lat[upd_rows]
    list(update = update, full = full, counts = .counts(site, g, rec, nrow(update)))
}

#' Process all BADM groups of a station.
#'
#' Each group is processed independently: an error is logged and recorded in the
#' counts, and the other groups continue.
#' @param badm named list of group tables.
#' @return named list of process_group() results.
process_badm <- function(badm, towers, plot_points, site = NA) {
    out <- list()
    for (g in CONVERTED_GROUPS) {
        df <- badm[[g]]
        if (is.null(df) || !nrow(df)) { log_info(g, ": no records"); next }
        res <- tryCatch(process_group(df, GROUP_SPECS[[g]], towers, plot_points, site),
                        error = function(e) {
                            log_error(g, ": ", conditionMessage(e))
                            cnt <- .counts(site, g, .rec_template(0), 0)
                            cnt$status <- "error"; cnt$message <- conditionMessage(e)
                            list(update = NULL, full = NULL, counts = cnt)
                        })
        if (!is.null(res)) {
            if (identical(res$counts$status, "ok"))
                log_info(g, ": ", res$counts$n_records, " records, ", res$counts$n_update,
                         " with calculated coordinates (", res$counts$n_magnetic, " magnetic, ",
                         res$counts$n_skipped, " skipped)")
            out[[g]] <- res
        }
    }
    out
}
