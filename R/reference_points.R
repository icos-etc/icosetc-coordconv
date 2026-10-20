# reference_points.R -- tower positions over time, plot hierarchy, record references.
#
# Inputs : GRP_LOCATION records, resolved plot points, record dates/plot codes.
# Outputs: reference point (id, lat, lon) for every record.
#
# Rules
#   Tower : the GRP_LOCATION record valid at the record date, i.e. the latest
#           LOCATION_DATE <= record date; if none, the earliest. A record without
#           a date gets the latest position (flagged). Undated LOCATION records
#           count as the earliest ones.
#   Plots : CP and SP-I are referenced to the tower; SP-II_<i>-<j> to SP-I_<i> (real IDs use a hyphen; an underscore is also accepted).
#   BADM  : tower groups (EC, BM, STO, INSTMAN) -> tower at the record date.
#           plot groups (GAI, TREE, LITTERPNT, WTDPNT, D_SNOW, SPPS):
#             <X>_PLOT NA/empty/"outside..." -> tower at the record date;
#             <X>_PLOT names a known plot    -> that plot's resolved position;
#             <X>_PLOT names an unknown plot -> no reference (record skipped).

#' Tower positions from GRP_LOCATION.
#'
#' @param location_df GRP_LOCATION table (character columns).
#' @return data.frame(location_id, date, lat, lon) ordered by date (undated first).
tower_positions <- function(location_df) {
    if (is.null(location_df) || !nrow(location_df) ||
        !all(c("LOCATION_LAT", "LOCATION_LONG") %in% names(location_df)))
        return(data.frame(location_id = character(0), date = as.Date(character(0)),
                          lat = numeric(0), lon = numeric(0)))
    out <- data.frame(location_id = col_or_na(location_df, "ID"),
                      date = parse_badm_date(col_or_na(location_df, "LOCATION_DATE")),
                      lat = as_num(location_df$LOCATION_LAT), lon = as_num(location_df$LOCATION_LONG),
                      stringsAsFactors = FALSE)
    out <- out[!is.na(out$lat) & !is.na(out$lon), , drop = FALSE]
    out[order(!is.na(out$date), out$date), , drop = FALSE]
}

#' Tower position valid at each date.
#'
#' @param towers output of tower_positions().
#' @param date Date vector.
#' @return data.frame(ref_id, ref_lat, ref_lon, no_date) with one row per date.
tower_at <- function(towers, date) {
    n <- length(date)
    out <- data.frame(ref_id = rep(NA_character_, n), ref_lat = NA_real_, ref_lon = NA_real_,
                      no_date = is.na(date), stringsAsFactors = FALSE)
    if (!nrow(towers) || !n) return(out)
    nt <- nrow(towers)
    idx <- rep(nt, n)                                   # no date -> latest
    dated <- !is.na(towers$date)
    has <- !is.na(date)
    if (any(has)) {
        if (any(dated)) {
            td <- which(dated)
            k <- findInterval(as.numeric(date[has]), as.numeric(towers$date[td]))
            idx[has] <- ifelse(k == 0, 1L, td[pmax(k, 1)])   # none <= date -> earliest
        } else {
            idx[has] <- 1L
        }
    }
    out$ref_id  <- paste0("TOWER:", towers$location_id[idx])
    out$ref_lat <- towers$lat[idx]
    out$ref_lon <- towers$lon[idx]
    out
}

#' Plot type from a PLOT_ID: "TOWER", "CP", "SP-I", "SP-II" or NA.
plot_type <- function(plot_id) {
    id <- toupper(trimws(plot_id))
    ifelse(grepl("^TOWER", id), "TOWER",
    ifelse(grepl("^CP_", id), "CP",
    ifelse(grepl("^SP-II_", id), "SP-II",
    ifelse(grepl("^SP-I_", id), "SP-I", NA_character_))))
}

#' Parent SP-I of an SP-II plot ("SP-II_03-02" or "SP-II_03_02" -> "SP-I_03"),
#' matched numerically.
#'
#' @param id SP-II plot IDs.
#' @param spi_ids available SP-I IDs.
#' @return parent SP-I ID (NA if not found).
sp2_parent <- function(id, spi_ids) {
    num_of <- function(x, re) suppressWarnings(as.integer(sub(re, "\\1", toupper(trimws(x)))))
    want <- num_of(id, "^SP-II_0*([0-9]+)[_-].*$")
    have <- num_of(spi_ids, "^SP-I_0*([0-9]+).*$")
    spi_ids[match(want, have)]
}

#' Normalise a plot code for matching (upper case, trimmed).
norm_plot_code <- function(x) toupper(trimws(as.character(x)))

#' Reference point of BADM records.
#'
#' @param plot_codes <X>_PLOT values (NULL/NA for tower groups).
#' @param dates record dates.
#' @param plot_points data.frame(plot_id, lat, lon) of resolved plots.
#' @param towers output of tower_positions().
#' @param reference "tower" or "plot".
#' @return data.frame(ref_id, ref_lat, ref_lon, ref_rule) where ref_rule is
#'   "tower", "tower_fallback" (no/outside plot), "plot" or "unknown_plot".
resolve_record_reference <- function(plot_codes, dates, plot_points, towers, reference) {
    n <- length(dates)
    tw <- tower_at(towers, dates)
    out <- data.frame(ref_id = tw$ref_id, ref_lat = tw$ref_lat, ref_lon = tw$ref_lon,
                      ref_rule = rep("tower", n), stringsAsFactors = FALSE)
    if (reference == "tower") return(out)
    code <- norm_plot_code(if (is.null(plot_codes)) rep(NA, n) else plot_codes)
    to_tower <- is.na(code) | !nzchar(code) | grepl("OUTSIDE", code) | code == "EC" | grepl("^TOWER", code)
    out$ref_rule[to_tower] <- "tower_fallback"
    m <- match(code, norm_plot_code(plot_points$plot_id))
    use_plot <- !to_tower & !is.na(m)
    out$ref_id[use_plot]  <- plot_points$plot_id[m[use_plot]]
    out$ref_lat[use_plot] <- plot_points$lat[m[use_plot]]
    out$ref_lon[use_plot] <- plot_points$lon[m[use_plot]]
    out$ref_rule[use_plot] <- "plot"
    unknown <- !to_tower & is.na(m)
    out$ref_id[unknown] <- code[unknown]
    out$ref_lat[unknown] <- NA
    out$ref_lon[unknown] <- NA
    out$ref_rule[unknown] <- "unknown_plot"
    out
}
