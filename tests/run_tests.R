# run_tests.R -- offline test suite (base R only; geosphere used as an optional cross-check).
#
# Usage: Rscript tests/run_tests.R
# Exit status 1 if any test fails. Each test prints PASS / FAIL / SKIP.

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
ROOT  <- normalizePath(file.path(dirname(if (length(.file)) .file else "tests/run_tests.R"), ".."))
for (f in list.files(file.path(ROOT, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
check_packages()
logger_init(NULL, echo = FALSE)

.results <- data.frame(test = character(0), result = character(0), note = character(0))
test <- function(name, expr) {
    res <- tryCatch({ force(expr); c("PASS", "") },
                    skip = function(s) c("SKIP", conditionMessage(s)),
                    error = function(e) c("FAIL", conditionMessage(e)))
    .results[nrow(.results) + 1, ] <<- c(name, res)
    cat(sprintf("%-4s  %s%s\n", res[1], name, if (nzchar(res[2])) paste0("  -- ", res[2]) else ""))
}
skip <- function(msg) stop(structure(class = c("skip", "condition"), list(message = msg, call = NULL)))
expect_lt <- function(x, tol, what = "value") {
    if (any(is.na(x)) || any(x >= tol)) stop(sprintf("%s: max %.3g >= %.3g", what, max(x), tol))
}
has_geosphere <- requireNamespace("geosphere", quietly = TRUE)

# Geodesic distance in metres (geosphere if available, else AEQD == geodesic)
geod_m <- function(lat1, lon1, lat2, lon2) {
    if (has_geosphere) geosphere::distGeo(cbind(lon1, lat1), cbind(lon2, lat2))
    else distance_m(lat1, lon1, lat2, lon2)
}

# Test grid: 3 latitudes x distances 1 m .. 2 km x 8 directions
set.seed(1)
GRID <- expand.grid(lat0 = c(5, 45, 78), d = c(1, 10, 100, 500, 2000), az = seq(10, 350, by = 45))
GRID$lon0 <- c(-60, 10, 25)[match(GRID$lat0, c(5, 45, 78))]

# ---- 1. Conversions ---------------------------------------------------------------

test("round trip lat/long -> N/E -> lat/long (< 1 mm)", {
    tgt <- ne_to_latlon(GRID$lat0, GRID$lon0, GRID$d * sin(GRID$az * pi / 180), GRID$d * cos(GRID$az * pi / 180))
    ne  <- latlon_to_ne(GRID$lat0, GRID$lon0, tgt$lat, tgt$lon)
    back <- ne_to_latlon(GRID$lat0, GRID$lon0, ne$east, ne$north)
    expect_lt(geod_m(tgt$lat, tgt$lon, back$lat, back$lon), 1e-3, "error (m)")
})

test("round trip polar -> lat/long -> polar (< 1 mm)", {
    ll  <- polar_to_latlon(GRID$lat0, GRID$lon0, GRID$d, GRID$az)
    pol <- latlon_to_polar(GRID$lat0, GRID$lon0, ll$lat, ll$lon)
    # error as position difference implied by the recovered polar pair
    back <- polar_to_latlon(GRID$lat0, GRID$lon0, pol$distance, pol$azimuth)
    expect_lt(geod_m(ll$lat, ll$lon, back$lat, back$lon), 1e-3, "position error (m)")
    expect_lt(abs(pol$distance - GRID$d), 1e-3, "distance error (m)")
    # azimuth error expressed as the arc length it implies at that distance
    daz <- ((pol$azimuth - GRID$az + 180) %% 360) - 180
    expect_lt(abs(daz * pi / 180 * GRID$d), 1e-3, "azimuth error as arc (m)")
})

test("AEQD polar vs geosphere geodesic destPoint/distGeo (< 1 mm)", {
    if (!has_geosphere) skip("geosphere not installed")
    ours <- polar_to_latlon(GRID$lat0, GRID$lon0, GRID$d, GRID$az)
    ref  <- geosphere::destPoint(cbind(GRID$lon0, GRID$lat0), GRID$az, GRID$d)
    expect_lt(geosphere::distGeo(cbind(ours$lon, ours$lat), ref), 1e-3, "position difference (m)")
    dg <- geosphere::distGeo(cbind(GRID$lon0, GRID$lat0), cbind(ours$lon, ours$lat))
    expect_lt(abs(dg - GRID$d), 1e-3, "distance difference (m)")
})

test("polar <-> N/E quadrants and azimuth wrap", {
    ne <- polar_to_ne(c(10, 10, 10, 10), c(45, 135, 225, 315))
    stopifnot(all.equal(sign(ne$east), c(1, 1, -1, -1)), all.equal(sign(ne$north), c(1, -1, -1, 1)))
    p <- ne_to_polar(c(0, 1, 0, -1, -1e-12), c(1, 0, -1, 0, 1))
    dang <- ((p$azimuth - c(0, 90, 180, 270, 0) + 180) %% 360) - 180
    stopifnot(all(abs(dang) < 1e-9), all(p$azimuth >= 0 & p$azimuth < 360))
})

# ---- 2. Magnetic corrections ---------------------------------------------------------

test("rotate_magnetic_ne known values (D = 0, +5, -5, quadrants)", {
    r0 <- rotate_magnetic_ne(c(3, -4), c(7, 2), 0)
    stopifnot(isTRUE(all.equal(r0$east, c(3, -4))), isTRUE(all.equal(r0$north, c(7, 2))))
    r <- rotate_magnetic_ne(c(0, 0, 100, 0), c(100, 100, 0, -100), c(5, -5, 5, 5))
    s <- 100 * sin(5 * pi / 180); co <- 100 * cos(5 * pi / 180)
    stopifnot(isTRUE(all.equal(r$east,  c(s, -s, co, -s), tolerance = 1e-12)),
              isTRUE(all.equal(r$north, c(co, co, -s, -co), tolerance = 1e-12)))
    stopifnot(isTRUE(all.equal(round(r$east[1], 3), 8.716)), isTRUE(all.equal(round(r$north[1], 3), 99.619)))
    rr <- rotate_magnetic_ne(c(12, -7), c(-3, 9), 17)
    stopifnot(isTRUE(all.equal(sqrt(rr$east^2 + rr$north^2), sqrt(c(12, -7)^2 + c(-3, 9)^2))))
})

test("magnetic point: polar and N/E inputs give the same position (< 1 mm)", {
    D <- c(4, -6.5); d <- 50; az_mag <- 30
    ne_mag <- polar_to_ne(d, az_mag)
    base <- data.frame(ref_lat = 45, ref_lon = 10, in_lat = NA_real_, in_lon = NA_real_, magdec = D)
    a <- compute_positions(transform(base, in_distance = d, in_azimuth = az_mag, in_east = NA_real_, in_north = NA_real_))
    b <- compute_positions(transform(base, in_distance = NA_real_, in_azimuth = NA_real_,
                                     in_east = ne_mag$east, in_north = ne_mag$north))
    expect_lt(geod_m(a$lat, a$lon, b$lat, b$lon), 1e-3, "difference (m)")
    stopifnot(isTRUE(all.equal(a$azimuth, (az_mag + D) %% 360)))
})

test("compute_positions format precedence polar > N/E > lat/long", {
    rec <- data.frame(ref_lat = 45, ref_lon = 10, magdec = 0,
                      in_distance = c(10, NA, NA, NA), in_azimuth = c(90, NA, NA, NA),
                      in_east = c(5, 5, NA, NA), in_north = c(0, 0, NA, NA),
                      in_lat = c(45.1, 45.1, 45.1, NA), in_lon = c(10, 10, 10, NA))
    r <- compute_positions(rec)
    stopifnot(identical(r$input_format, c("polar", "ne", "latlong", NA)),
              identical(r$n_formats, c(3L, 2L, 1L, 0L)),
              !is.na(r$format_disagreement_m[1]), !is.na(r$format_disagreement_m[2]),
              is.na(r$format_disagreement_m[3]))
    # row 2: N/E (5, 0) vs the lat/long 0.1 deg north -> ~11.1 km apart
    stopifnot(abs(r$format_disagreement_m[2] - 11120) < 50)
})

# ---- 3. Declination ----------------------------------------------------------------

test("declination vs NOAA NCEI IGRF-14 reference (45N 10E 2024-07-01, tol 0.05 deg)", {
    if (!oce_available()) skip("oce not available")
    # NOAA NCEI geomag calculator, model IGRF (DGRF2020/IGRF-14), 2024.4973: 3.57317 deg
    D <- declination(45, 10, as.Date("2024-07-01"))
    stopifnot(abs(D - 3.57317) < 0.05)
})

test("declination regression vs oce 1.8-4 value (tol 1e-3 deg)", {
    if (!oce_available()) skip("oce not available")
    stopifnot(abs(declination(45, 10, as.Date("2024-07-01")) - 3.573298466) < 1e-3)
})

test("declination: one site, several dates -> one value per date (oce pitfall)", {
    if (!oce_available()) skip("oce not available")
    # oce::magneticField silently uses only the first time when lat/lon are scalars
    d <- declination(45, 10, as.Date(c("2002-01-01", "2024-07-01")))
    stopifnot(length(d) == 2, abs(d[2] - 3.573298466) < 1e-3, abs(d[1] - d[2]) > 0.5)
})

test("declination: missing oce gives an explicit, actionable error", {
    old <- options(coordconv.force_no_oce = TRUE)
    msg <- tryCatch({ declination(45, 10, as.Date("2024-07-01")); "" }, error = conditionMessage)
    options(old)
    stopifnot(grepl("oce", msg), grepl("magnetic", msg, ignore.case = TRUE))
})

test("declination: date checks (< 1900 error, > 2030 warning)", {
    if (!oce_available()) skip("oce not available")
    e <- tryCatch(declination(45, 10, as.Date("1850-01-01")), error = function(e) "err")
    stopifnot(identical(e, "err"))
    w <- tryCatch(declination(45, 10, as.Date("2031-06-01")), warning = function(w) "warn")
    stopifnot(identical(w, "warn"))
})

# ---- 4. Utilities --------------------------------------------------------------------

test("parse_badm_date (YYYY, YYYYMM, YYYYMMDD, YYYYMMDDHHMM, invalid)", {
    d <- parse_badm_date(c("2015", "201503", "20150317", "201503171230", "2015-03-17", "", NA, "20151399"))
    stopifnot(identical(d, as.Date(c("2015-07-01", "2015-03-15", "2015-03-17", "2015-03-17",
                                     "2015-03-17", NA, NA, NA))))
    stopifnot(identical(parse_badm_date(20150317), as.Date("2015-03-17")))
})

test("eco recode: known codes map, unknown code errors", {
    stopifnot(identical(recode_eco(c("ENF", "cro", "WET", "GRA", NA, "URB", "WAT")),
                        c("forest", "crop", "mire", "grassland", NA, "urban", "water")))
    e <- tryCatch(recode_eco("XYZ"), error = function(e) "err")
    stopifnot(identical(e, "err"))
})

test("rect_corners_ne: sides, right angles, 180-degree symmetry, north alignment", {
    m <- rect_corners_ne(5, 10, 30)
    side <- function(a, b) sqrt(sum((m[a, ] - m[b, ])^2))
    stopifnot(isTRUE(all.equal(c(side("TR", "BR"), side("BR", "BL"), side("BL", "TL"), side("TL", "TR")),
                               c(10, 5, 10, 5))))
    v1 <- m["BR", ] - m["TR", ]; v2 <- m["TL", ] - m["TR", ]
    stopifnot(abs(sum(v1 * v2)) < 1e-9)
    # long side points to the azimuth: TR - BR has azimuth 30
    stopifnot(isTRUE(all.equal(ne_to_polar(v1[1] * -1, v1[2] * -1)$azimuth, 30)))
    m180 <- rect_corners_ne(5, 10, 210)
    stopifnot(isTRUE(all.equal(sort(round(m180, 9)), sort(round(m, 9)))))
    m0 <- rect_corners_ne(5, 10, 0)
    stopifnot(isTRUE(all.equal(unname(m0["TR", ]), c(2.5, 5))), isTRUE(all.equal(unname(m0["BL", ]), c(-2.5, -5))))
})

# ---- 5. Reference points and group specs ------------------------------------------------

TOWERS <- data.frame(location_id = c("L1", "L2"), date = as.Date(c("2015-01-01", "2020-01-01")),
                     lat = c(46, 46.0001), lon = c(11, 11.0001))

test("tower_at: before first, between, after last, no date (flag)", {
    t <- tower_at(TOWERS, as.Date(c("2010-05-01", "2017-01-01", "2021-01-01", NA)))
    stopifnot(identical(t$ref_id, c("TOWER:L1", "TOWER:L1", "TOWER:L2", "TOWER:L2")),
              identical(t$no_date, c(FALSE, FALSE, FALSE, TRUE)))
    t2 <- tower_at(TOWERS, as.Date("2020-01-01"))           # on the date itself -> that position
    stopifnot(identical(t2$ref_id, "TOWER:L2"))
})

test("sp2_parent: SP-II -> SP-I (hyphen/underscore, zero padding, missing parent)", {
    p <- sp2_parent(c("SP-II_01-02", "SP-II_03_01", "SP-II_09-01"), c("SP-I_01", "SP-I_3"))
    stopifnot(identical(p, c("SP-I_01", "SP-I_3", NA)))
})

test("record reference: plot, empty/outside/EC -> tower fallback, unknown plot", {
    pts <- data.frame(plot_id = c("CP_01", "SP-I_02"), lat = c(46.001, 46.002), lon = c(11.001, 11.002))
    r <- resolve_record_reference(c("CP_01", NA, "Outside CP", "EC", "cp_01 ", "CP_99", "SP-I_02"),
                                  as.Date(rep("2021-06-01", 7)), pts, TOWERS, "plot")
    stopifnot(identical(r$ref_rule, c("plot", "tower_fallback", "tower_fallback", "tower_fallback",
                                      "plot", "unknown_plot", "plot")),
              identical(r$ref_id[c(2, 6)], c("TOWER:L2", "CP_99")), is.na(r$ref_lat[6]),
              isTRUE(all.equal(r$ref_lat[7], 46.002)))
    t <- resolve_record_reference(c("CP_01", NA), as.Date(c("2016-01-01", NA)), pts, TOWERS, "tower")
    stopifnot(identical(t$ref_id, c("TOWER:L1", "TOWER:L2")))
})

test("group specs: required fields and legacy output column names", {
    req <- c("group", "date_cols", "ne", "polar", "northref_cols", "reference", "out_lat", "out_lon")
    for (g in c("GRP_PLOT", CONVERTED_GROUPS)) {
        miss <- setdiff(req, names(GROUP_SPECS[[g]]))
        if (length(miss)) stop(g, " lacks ", paste(miss, collapse = ","))
        if (GROUP_SPECS[[g]]$group != g) stop("group name mismatch in ", g)
    }
    legacy <- c(GRP_EC = "EC", GRP_BM = "BM", GRP_STO = "STO_PROF", GRP_INSTMAN = "INSTMAN", GRP_SPPS = "SPPS",
                GRP_GAI = "GAI", GRP_TREE = "TREE", GRP_LITTERPNT = "LITTERPNT", GRP_WTDPNT = "WTDPNT",
                GRP_D_SNOW = "D_SNOW")
    stopifnot(setequal(names(legacy), CONVERTED_GROUPS))
    for (g in names(legacy))
        stopifnot(identical(c(GROUP_SPECS[[g]]$out_lat, GROUP_SPECS[[g]]$out_lon),
                            paste0(legacy[[g]], c("_LOCATION_CALC_LAT", "_LOCATION_CALC_LONG"))))
    stopifnot(identical(GROUP_SPECS$GRP_EC$extra_out, "EC_SA_NORTH_MAGDEC"),
              length(GROUP_SPECS$GRP_GAI$ne) == 2, length(GROUP_SPECS$GRP_GAI$polar) == 2)
})

# ---- 6. Processing -----------------------------------------------------------------------

test("process_group: magnetic N/E record rotated at its tower/date, update columns as legacy", {
    if (!oce_available()) skip("oce not available")
    df <- data.frame(ID = c("1", "2", "3", ""), EC_DATE = c("20190301", "20210301", "", "20210301"),
                     EC_EASTWARD_DIST = c("10", "10", "5", "1"), EC_NORTHWARD_DIST = c("0", "0", "5", "1"),
                     EC_NORTHREF = c("Geographic", "Magnetic", "Magnetic", NA), stringsAsFactors = FALSE)
    r <- process_group(df, GROUP_SPECS$GRP_EC, TOWERS, NULL, "XX-Tst")
    stopifnot(identical(names(r$update), c("ID", "EC_SA_NORTH_MAGDEC", "EC_LOCATION_CALC_LONG", "EC_LOCATION_CALC_LAT")),
              identical(r$update$ID, c("1", "2")),                     # 3: magnetic without date skipped; "" dropped
              r$full$status[3] == "skipped: magnetic record without date")
    D <- declination(46.0001, 11.0001, as.Date("2021-03-01"))
    e <- rotate_magnetic_ne(10, 0, D)
    ll <- ne_to_latlon(46.0001, 11.0001, e$east, e$north)
    stopifnot(abs(r$update$EC_SA_NORTH_MAGDEC[2] - D) < 1e-9, r$update$EC_SA_NORTH_MAGDEC[1] == 0)
    expect_lt(geod_m(r$update$EC_LOCATION_CALC_LAT[2], r$update$EC_LOCATION_CALC_LONG[2], ll$lat, ll$lon), 1e-3, "pos (m)")
    ll1 <- ne_to_latlon(46, 11, 10, 0)                                  # 2019 -> first tower position
    expect_lt(geod_m(r$update$EC_LOCATION_CALC_LAT[1], r$update$EC_LOCATION_CALC_LONG[1], ll1$lat, ll1$lon), 1e-3, "pos (m)")
})

test("process_group: GAI DHP/CEPT pairs coalesced, plot reference and unknown plot", {
    pts <- data.frame(plot_id = "CP_01", lat = 46.001, lon = 11.001)
    df <- data.frame(ID = c("1", "2", "3"), GAI_DATE = "2021", GAI_PLOT = c("CP_01", "CP_01", "CP_07"),
                     GAI_DHP_EASTWARD_DIST = c("3", NA, "1"), GAI_DHP_NORTHWARD_DIST = c("4", NA, "1"),
                     GAI_CEPT_DISTANCE_POLAR = c(NA, "10", NA), GAI_CEPT_ANGLE_POLAR = c(NA, "90", NA),
                     stringsAsFactors = FALSE)
    r <- process_group(df, GROUP_SPECS$GRP_GAI, TOWERS, pts, "XX-Tst")
    stopifnot(identical(r$full$input_format, c("ne", "polar", "ne")), identical(r$update$ID, c("1", "2")),
              r$full$status[3] == "skipped: unknown reference plot")
    p2 <- polar_to_latlon(46.001, 11.001, 10, 90)
    expect_lt(geod_m(r$update$GAI_LOCATION_CALC_LAT[2], r$update$GAI_LOCATION_CALC_LONG[2], p2$lat, p2$lon), 1e-3, "pos (m)")
})

test("process_plots: hierarchy (SP-II from resolved SP-I), update only N/E/polar plots, areas", {
    pl <- data.frame(ID = as.character(1:6), PLOT_ID = c("TOWER", "CP_01", "SP-I_01", "SP-II_01-01", "SP-II_01-02", "CP_02"),
                     PLOT_LOCATION_LAT = c("46", "46.0003", NA, "46.0012", NA, NA),
                     PLOT_LOCATION_LONG = c("11", "11.0002", NA, "11.0001", NA, NA),
                     PLOT_EASTWARD_DIST = c(NA, NA, "100", NA, "5", NA), PLOT_NORTHWARD_DIST = c(NA, NA, "50", NA, "-3", NA),
                     PLOT_DISTANCE_POLAR = c(NA, NA, NA, NA, NA, "40"), PLOT_ANGLE_POLAR = c(NA, NA, NA, NA, NA, "210"),
                     PLOT_RECT_X = c(NA, NA, NA, NA, NA, "10"), PLOT_RECT_Y = c(NA, NA, NA, NA, NA, "20"),
                     PLOT_RECT_DIR = c(NA, NA, NA, NA, NA, "30"), PLOT_DATE = "20210601", stringsAsFactors = FALSE)
    r <- process_plots(pl, TOWERS, "forest", "XX-Tst")
    stopifnot(identical(sort(r$update$ID), c("3", "5", "6")))     # TOWER, CP_01, SP-II_01-01 are lat/long
    spi <- ne_to_latlon(46.0001, 11.0001, 100, 50)
    sp2 <- ne_to_latlon(spi$lat, spi$lon, 5, -3)
    u <- r$update[r$update$ID == "5", ]
    expect_lt(geod_m(u$PLOT_LOCATION_LAT, u$PLOT_LOCATION_LONG, sp2$lat, sp2$lon), 1e-3, "SP-II pos (m)")
    a <- r$areas
    area_of <- function(id) { g <- a[a$ID == id, ]; c <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(g)))
                              as.numeric(sf::st_area(sf::st_transform(g, sf::st_crs(aeqd_crs(c[2], c[1]))))) }
    stopifnot(abs(area_of("CP_01") - 2000) / 2000 < 0.005,         # forest circle 2000 m2 (64-gon)
              abs(area_of("CP_02") - 200) < 0.01,                  # rectangle 10 x 20
              abs(area_of("SP-I_01") - pi * 15^2) / (pi * 15^2) < 0.005)
})

# ---- 7. QC flags -------------------------------------------------------------------------

test("QC: N/E vs polar disagreement metric ignores lat/long", {
    rec <- data.frame(ref_lat = 46, ref_lon = 11, magdec = 0, in_distance = c(5, 5, NA), in_azimuth = c(90, 90, NA),
                      in_east = c(5, 2, 3), in_north = c(0, 0, 4), in_lat = c(NA, 46.01, 46.01), in_lon = c(NA, 11, 11))
    r <- compute_positions(rec)
    stopifnot(abs(r$ne_polar_disagreement_m[1]) < 1e-9, abs(r$ne_polar_disagreement_m[2] - 3) < 1e-9,
              is.na(r$ne_polar_disagreement_m[3]))
})

test("QC flags: FORMAT_DISAGREEMENT, LARGE_SHIFT (API and --previous), NOT_COMPUTED, OUTSIDE_TA", {
    tw <- list(lat = 46, lon = 11)
    pos <- ne_to_latlon(46, 11, c(10, 10, 100, 10), c(0, 0, 0, 0))
    prev <- ne_to_latlon(46, 11, c(10, 40, 100, NA), c(0, 0, 0, NA))
    full <- data.frame(ID = c("1", "2", "3", "4", "5"), status = c("ok", "ok", "ok", "ok", "skipped: unknown reference plot"),
                       input_format = c("polar", "ne", "ne", "ne", "ne"), magdec_applied = c(0, 2.5, 0, 0, NA),
                       ne_polar_disagreement_m = c(3.2, NA, NA, NA, NA),
                       calc_lat = c(pos$lat, NA), calc_lon = c(pos$lon, NA),
                       db_TREE_LOCATION_CALC_LAT = c(prev$lat, NA), db_TREE_LOCATION_CALC_LONG = c(prev$lon, NA))
    ta <- sf::st_sf(ID = "TA", geometry = sf::st_sfc(polygon_from_ne(46, 11, square_ne(100)), crs = 4326))  # +-50 m
    prv <- data.frame(ID = "4", lat = ne_to_latlon(46, 11, 10, 50)$lat, lon = ne_to_latlon(46, 11, 10, 50)$lon)
    f <- qc_group_flags("XX-Tst", "GRP_TREE", full, QC_DEFAULTS, tw, ta = ta, previous = prv)
    key <- paste(f$ID, f$flag)
    stopifnot(setequal(key, c("1 FORMAT_DISAGREEMENT", "2 LARGE_SHIFT", "4 LARGE_SHIFT", "5 NOT_COMPUTED", "3 OUTSIDE_TA")))
    stopifnot(grepl("shift_m=30.00; D_applied=2.500; previous=api", f$detail[key == "2 LARGE_SHIFT"], fixed = TRUE),
              grepl("previous=previous_dir", f$detail[key == "4 LARGE_SHIFT"], fixed = TRUE),
              abs(as.numeric(sub(".*=", "", f$detail[key == "3 OUTSIDE_TA"])) - 50) < 0.5)
    # thresholds come from settings
    f2 <- qc_group_flags("XX-Tst", "GRP_TREE", full, qc_settings(list(qc = list(format_disagreement_max_m = 5))), tw)
    stopifnot(!any(f2$flag == "FORMAT_DISAGREEMENT"))
})

test("GAI instrument rule: method picks the pair; fallback, mismatch and within-pair disagreement", {
    pts <- data.frame(plot_id = "CP_01", lat = 46.001, lon = 11.001)
    df <- data.frame(
        ID     = c("1", "2", "3", "4", "5"),
        GAI_DATE = "20220926", GAI_PLOT = "CP_01",
        GAI_METHOD = c("SUNSCAN", "DHP", "Leaf counting", "DHP", "sunscan "),
        GAI_DHP_EASTWARD_DIST  = c(NA, "3", NA, NA, NA), GAI_DHP_NORTHWARD_DIST = c(NA, "4", NA, NA, NA),
        GAI_DHP_DISTANCE_POLAR = c("11.53", NA, "11.53", NA, NA), GAI_DHP_ANGLE_POLAR = c("270.97", NA, "270.97", NA, NA),
        GAI_CEPT_EASTWARD_DIST = c("0.773", "7", "2", NA, "5"), GAI_CEPT_NORTHWARD_DIST = c("21.086", "7", "2", NA, "0"),
        GAI_CEPT_DISTANCE_POLAR = c(NA, NA, NA, "6", "10"), GAI_CEPT_ANGLE_POLAR = c(NA, NA, NA, "90", "90"),
        stringsAsFactors = FALSE)
    r <- process_group(df, GROUP_SPECS$GRP_GAI, TOWERS, pts, "XX-Tst")
    f <- r$full
    # 1: SUNSCAN with a copied DHP polar and its own CEPT N/E -> CEPT N/E position
    p1 <- ne_to_latlon(46.001, 11.001, 0.773, 21.086)
    stopifnot(f$instrument_pair[1] == "CEPT", f$input_format[1] == "ne", is.na(f$pair_flag[1]))
    expect_lt(geod_m(f$calc_lat[1], f$calc_lon[1], p1$lat, p1$lon), 1e-3, "SUNSCAN pos (m)")
    # 2: DHP record -> DHP pair even though CEPT N/E is also filled
    stopifnot(f$instrument_pair[2] == "DHP", f$in_east_dist[2] == 3, f$in_north_dist[2] == 4)
    # 3: unknown method -> fallback (first complete pair, polar > N/E) + flag; no cross-pair disagreement
    stopifnot(f$pair_flag[3] == "GAI_METHOD_UNKNOWN", f$instrument_pair[3] == "DHP", f$input_format[3] == "polar",
              is.na(f$ne_polar_disagreement_m[3]))
    # 4: DHP method but only the CEPT pair is filled -> CEPT used + flag
    stopifnot(f$pair_flag[4] == "GAI_PAIR_MISMATCH", f$instrument_pair[4] == "CEPT", grepl("pair_used=CEPT", f$pair_detail[4]))
    # 5: within-pair N/E (5,0) vs polar (10 @ 90) -> 5 m disagreement; case-insensitive method
    stopifnot(f$instrument_pair[5] == "CEPT", abs(f$ne_polar_disagreement_m[5] - 5) < 1e-6)
    q <- qc_group_flags("XX-Tst", "GRP_GAI", f, QC_DEFAULTS, list(lat = 46, lon = 11))
    stopifnot(setequal(paste(q$ID, q$flag), c("3 GAI_METHOD_UNKNOWN", "4 GAI_PAIR_MISMATCH", "5 FORMAT_DISAGREEMENT")))
    stopifnot(all(c("1", "2", "3", "4", "5") %in% r$update$ID))   # flagged records are not withheld
})

# ---- 8. End to end -----------------------------------------------------------------------

COORDCONV_ROOT <- ROOT
source(file.path(ROOT, "examples", "run_example.R"))

test("end to end: examples/ (offline) reproduces examples/expected", {
    if (!oce_available()) skip("oce not available (the example has magnetic records)")
    r <- run_example()
    if (r$status != "PASS") stop(paste(r$details, collapse = "; "))
})

test("end to end without oce: magnetic records stop with the explicit declination error", {
    old <- options(coordconv.force_no_oce = TRUE)
    r <- tryCatch(run_example(), error = function(e) list(status = "FAIL", details = conditionMessage(e)))
    options(old)
    if (r$status != "SKIP") stop(paste(r$details, collapse = "; "))
})

# ---- Summary ---------------------------------------------------------------------------

cat(sprintf("\n%d tests: %d passed, %d failed, %d skipped\n", nrow(.results),
            sum(.results$result == "PASS"), sum(.results$result == "FAIL"), sum(.results$result == "SKIP")))
if (any(.results$result == "FAIL")) quit(status = 1)
