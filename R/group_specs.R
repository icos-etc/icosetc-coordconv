# group_specs.R -- column specification of the sampling plots and of every BADM group.
#
# One entry per group. Fields:
#   group          BADM variable group (API name) or "GRP_PLOT"
#   date_cols      record date, first non-empty column wins
#   plot_col       column naming the reference plot (plot-referenced groups)
#   ne             list of (east, north) column pairs, first complete pair wins
#   polar          list of (distance, azimuth) column pairs, first complete pair wins
#   latlong        (lat, long) columns of coordinates submitted as lat/long
#   northref_cols  north-reference columns ("Magnetic" -> declination applied),
#                  first non-empty wins; absent columns mean true north. As of
#                  2026-10 the API provides NORTHREF only for GRP_PLOT, EC, BM and
#                  INSTMAN; the other groups are therefore true north (their
#                  NORTHREF names are kept so a future column is picked up).
#   reference      "tower", "plot", "hierarchy" (GRP_PLOT) or "provider" (GRP_LOCATION)
#   out_lat/out_lon columns of the calculated coordinates in the update file
#   extra_out      extra update-file columns (EC: declination applied)
#   db_cols        existing DB values echoed in the full file for comparison
#
# Reference rule for "plot" groups (see reference_points.R): <X>_PLOT empty,
# "outside...", "EC" or "TOWER" -> tower valid at the record date; a known plot
# ID -> that plot's resolved position; an unknown plot ID -> record skipped.

.badm_spec <- function(prefix, reference, ne = TRUE, polar = TRUE, ne_prefix = prefix,
                       out_prefix = prefix, plot = reference == "plot", northref = TRUE) {
    list(
        group = paste0("GRP_", prefix),
        date_cols = paste0(prefix, c("_DATE", "_DATE_START")),
        plot_col = if (plot) paste0(prefix, "_PLOT") else NULL,
        ne = if (ne) list(paste0(ne_prefix, c("_EASTWARD_DIST", "_NORTHWARD_DIST"))) else list(),
        polar = if (polar) list(paste0(prefix, c("_DISTANCE_POLAR", "_ANGLE_POLAR"))) else list(),
        latlong = NULL,
        northref_cols = if (northref) paste0(prefix, "_NORTHREF") else character(0),
        reference = reference,
        out_lat = paste0(out_prefix, "_LOCATION_CALC_LAT"),
        out_lon = paste0(out_prefix, "_LOCATION_CALC_LONG"),
        extra_out = character(0),
        db_cols = paste0(out_prefix, c("_LOCATION_CALC_LAT", "_LOCATION_CALC_LONG"))
    )
}

# Replace fields of a spec (modifyList() would merge nested unnamed lists).
.set_fields <- function(spec, fields) {
    for (nm in names(fields)) spec[nm] <- list(fields[[nm]])
    spec
}

GROUP_SPECS <- list(
    # 1. Sampling plots (kplot table). Rectangular CPs: PLOT_RECT_X = short side,
    #    PLOT_RECT_Y = long side, PLOT_RECT_DIR = azimuth of the long side (degrees,
    #    corrected with the declination when PLOT_NORTHREF is "Magnetic"); the CP
    #    point is the rectangle centre. See rect_corners_ne() in convert.R.
    GRP_PLOT = list(
        group = "GRP_PLOT", date_cols = "PLOT_DATE", plot_col = "PLOT_ID",
        ne = list(c("PLOT_EASTWARD_DIST", "PLOT_NORTHWARD_DIST")),
        polar = list(c("PLOT_DISTANCE_POLAR", "PLOT_ANGLE_POLAR")),
        latlong = c("PLOT_LOCATION_LAT", "PLOT_LOCATION_LONG"),
        northref_cols = "PLOT_NORTHREF", reference = "hierarchy",
        rect_cols = c(x = "PLOT_RECT_X", y = "PLOT_RECT_Y", dir = "PLOT_RECT_DIR"),
        out_lat = "PLOT_LOCATION_LAT", out_lon = "PLOT_LOCATION_LONG",
        extra_out = character(0), db_cols = c("PLOT_LOCATION_LAT", "PLOT_LOCATION_LONG")
    ),
    # 2. Tower positions: reference provider only, never converted.
    GRP_LOCATION = list(
        group = "GRP_LOCATION", date_cols = "LOCATION_DATE", reference = "provider",
        latlong = c("LOCATION_LAT", "LOCATION_LONG")
    ),
    # 3-6. Tower-referenced groups (N/E offsets only)
    GRP_EC      = .set_fields(.badm_spec("EC", "tower", polar = FALSE),
                             list(extra_out = "EC_SA_NORTH_MAGDEC",
                                  db_cols = c("EC_LOCATION_CALC_LAT", "EC_LOCATION_CALC_LONG", "EC_SA_NORTH_MAGDEC"))),
    GRP_BM      = .badm_spec("BM", "tower", polar = FALSE),
    GRP_STO     = .badm_spec("STO", "tower", polar = FALSE, ne_prefix = "STO_PROF", out_prefix = "STO_PROF"),
    GRP_INSTMAN = .badm_spec("INSTMAN", "tower", polar = FALSE),
    # 7. SPPS: polar (own column names) or submitted lat/long; no N/E offsets
    GRP_SPPS    = .set_fields(.badm_spec("SPPS", "plot", ne = FALSE),
                             list(polar = list(c("SPPS_LOCATION_DIST", "SPPS_LOCATION_ANG")),
                                  latlong = c("SPPS_LOCATION_LAT", "SPPS_LOCATION_LON"))),
    # 8. GAI: hemispherical photo (DHP) or ceptometer (CEPT) positions. The pair is chosen
    #    from GAI_METHOD (see GAI_METHOD_PAIRS); within the pair polar > N/E. Unknown method:
    #    first complete pair, DHP then CEPT (the order of ne/polar below).
    GRP_GAI     = .set_fields(.badm_spec("GAI", "plot"),
                             list(ne = list(c("GAI_DHP_EASTWARD_DIST", "GAI_DHP_NORTHWARD_DIST"),
                                            c("GAI_CEPT_EASTWARD_DIST", "GAI_CEPT_NORTHWARD_DIST")),
                                  polar = list(c("GAI_DHP_DISTANCE_POLAR", "GAI_DHP_ANGLE_POLAR"),
                                               c("GAI_CEPT_DISTANCE_POLAR", "GAI_CEPT_ANGLE_POLAR")),
                                  method_col = "GAI_METHOD",
                                  instrument_pairs = list(
                                      DHP  = list(ne = c("GAI_DHP_EASTWARD_DIST", "GAI_DHP_NORTHWARD_DIST"),
                                                  polar = c("GAI_DHP_DISTANCE_POLAR", "GAI_DHP_ANGLE_POLAR")),
                                      CEPT = list(ne = c("GAI_CEPT_EASTWARD_DIST", "GAI_CEPT_NORTHWARD_DIST"),
                                                  polar = c("GAI_CEPT_DISTANCE_POLAR", "GAI_CEPT_ANGLE_POLAR"))),
                                  northref_cols = c("GAI_NORTHREF", "GAI_DHP_NORTHREF", "GAI_CEPT_NORTHREF"))),
    # 9-12. Plot-referenced groups (N/E or polar)
    GRP_TREE      = .badm_spec("TREE", "plot"),
    GRP_LITTERPNT = .badm_spec("LITTERPNT", "plot"),
    GRP_WTDPNT    = .badm_spec("WTDPNT", "plot"),
    GRP_D_SNOW    = .badm_spec("D_SNOW", "plot")
)

#' GAI_METHOD -> instrument pair whose position describes the measurement.
#'
#' Matching is case-insensitive on the trimmed value. Values found in the API (2026-10):
#'
#' | GAI_METHOD           | pair    | why                                                  |
#' |----------------------|---------|------------------------------------------------------|
#' | DHP                  | DHP     | digital hemispherical photography                    |
#' | SUNSCAN              | CEPT    | Delta-T SunScan ceptometer                           |
#' | ACCUPAR              | CEPT    | METER AccuPAR ceptometer                             |
#' | Modified VGA         | unknown | no instrument position (none in the data either)     |
#' | Visual estimation    | unknown | idem                                                 |
#' | Destructive          | unknown | idem                                                 |
#' | Spectral reflectance | unknown | idem                                                 |
#'
#' Any other value is "unknown": the first complete pair is used and GAI_METHOD_UNKNOWN is flagged.
GAI_METHOD_PAIRS <- c("DHP" = "DHP", "SUNSCAN" = "CEPT", "ACCUPAR" = "CEPT",
                      "MODIFIED VGA" = "unknown", "VISUAL ESTIMATION" = "unknown",
                      "DESTRUCTIVE" = "unknown", "SPECTRAL REFLECTANCE" = "unknown")

#' BADM groups that are downloaded (all specs except the plot table).
BADM_GROUPS <- setdiff(names(GROUP_SPECS), "GRP_PLOT")

#' BADM groups that are converted (all except the reference provider).
CONVERTED_GROUPS <- setdiff(BADM_GROUPS, "GRP_LOCATION")

#' Plot-area rules (radius / side in metres).
PLOT_AREA_RULES <- list(
    spi_radius   = c(forest = 15, other = 10),   # SP-I circular area
    cp_forest_r  = sqrt(2000 / pi),             # circle of 2000 m^2 (~25.23 m)
    cp_crop_r    = sqrt(5^2 + 10^2) / 2,        # half-diagonal of a 5 x 10 m plot (~5.59 m)
    cp_mire_side = 0.6                          # north-aligned 0.6 x 0.6 m square
)
