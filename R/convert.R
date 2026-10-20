# convert.R -- pure coordinate conversions (no I/O, no globals).
#
# Conventions
#   - lat/long: WGS84 geographic coordinates, decimal degrees.
#   - north/east offsets (east, north): metres, TRUE north, measured in the local
#     azimuthal equidistant (AEQD) plane centred on the reference point.
#   - polar (distance, azimuth): metres; degrees clockwise from TRUE north, [0, 360).
#
# The AEQD projection centred on the reference point maps every point to
# (E, N) = d * (sin a, cos a), where d and a are the geodesic distance and the
# initial geodesic azimuth from the reference (PROJ uses geodesic computations
# for the ellipsoidal AEQD). Polar <-> N/E is therefore exact, and N/E offsets
# are true-north offsets in the tangent plane at the reference point.
# Only PROJ.4-style strings are used, so this works from PROJ 6.3 onwards.

LONGLAT_CRS <- "+proj=longlat +datum=WGS84 +no_defs"

#' PROJ string of the AEQD CRS centred on a point.
#'
#' @param lat0,lon0 centre (decimal degrees).
#' @return character(1).
aeqd_crs <- function(lat0, lon0) {
    sprintf("+proj=aeqd +lat_0=%.10f +lon_0=%.10f +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
            lat0, lon0)
}

# Apply fn(rows, lat0, lon0) once per distinct reference point; rows with NA
# reference or NA inputs are skipped. Returns a 2-column matrix.
.by_reference <- function(ref_lat, ref_lon, ok, fn) {
    n <- length(ref_lat)
    out <- matrix(NA_real_, n, 2)
    ok <- ok & !is.na(ref_lat) & !is.na(ref_lon)
    if (!any(ok)) return(out)
    key <- paste(sprintf("%.10f", ref_lat), sprintf("%.10f", ref_lon))
    for (k in unique(key[ok])) {
        rows <- which(ok & key == k)
        out[rows, ] <- fn(rows, ref_lat[rows[1]], ref_lon[rows[1]])
    }
    out
}

#' North/east offsets -> lat/long.
#'
#' @param ref_lat,ref_lon reference point(s), recycled to the length of east.
#' @param east,north offsets in metres (true north).
#' @return data.frame(lat, lon).
ne_to_latlon <- function(ref_lat, ref_lon, east, north) {
    n <- length(east)
    ref_lat <- rep_len(ref_lat, n); ref_lon <- rep_len(ref_lon, n)
    m <- .by_reference(ref_lat, ref_lon, !is.na(east) & !is.na(north), function(rows, la, lo) {
        sf::sf_project(aeqd_crs(la, lo), LONGLAT_CRS, cbind(east[rows], north[rows]))
    })
    data.frame(lat = m[, 2], lon = m[, 1])
}

#' Lat/long -> north/east offsets from the reference point.
#'
#' @return data.frame(east, north) in metres.
latlon_to_ne <- function(ref_lat, ref_lon, lat, lon) {
    n <- length(lat)
    ref_lat <- rep_len(ref_lat, n); ref_lon <- rep_len(ref_lon, n)
    m <- .by_reference(ref_lat, ref_lon, !is.na(lat) & !is.na(lon), function(rows, la, lo) {
        sf::sf_project(LONGLAT_CRS, aeqd_crs(la, lo), cbind(lon[rows], lat[rows]))
    })
    data.frame(east = m[, 1], north = m[, 2])
}

#' Polar -> north/east: E = d sin(a), N = d cos(a).
#'
#' @param distance metres; @param azimuth degrees clockwise from north.
#' @return data.frame(east, north).
polar_to_ne <- function(distance, azimuth) {
    a <- azimuth * pi / 180
    data.frame(east = distance * sin(a), north = distance * cos(a))
}

#' North/east -> polar: d = sqrt(E^2 + N^2), a = atan2(E, N) in [0, 360).
#'
#' @return data.frame(distance, azimuth).
ne_to_polar <- function(east, north) {
    az <- (atan2(east, north) * 180 / pi) %% 360
    az[!is.na(az) & az >= 360] <- 0
    data.frame(distance = sqrt(east^2 + north^2), azimuth = az)
}

#' Polar -> lat/long (composition of polar_to_ne and ne_to_latlon).
polar_to_latlon <- function(ref_lat, ref_lon, distance, azimuth) {
    ne <- polar_to_ne(distance, azimuth)
    ne_to_latlon(ref_lat, ref_lon, ne$east, ne$north)
}

#' Lat/long -> polar (composition of latlon_to_ne and ne_to_polar).
latlon_to_polar <- function(ref_lat, ref_lon, lat, lon) {
    ne <- latlon_to_ne(ref_lat, ref_lon, lat, lon)
    ne_to_polar(ne$east, ne$north)
}

#' Geodesic-equivalent distance (m) between two sets of points.
#'
#' Computed as the AEQD distance from the first point, which equals the
#' geodesic distance (used for comparisons and tests).
distance_m <- function(lat1, lon1, lat2, lon2) {
    ne <- latlon_to_ne(lat1, lon1, lat2, lon2)
    sqrt(ne$east^2 + ne$north^2)
}

# ---- Positions of records -----------------------------------------------------

#' Compute positions of records given their reference point and inputs.
#'
#' The input format of each record is chosen with precedence
#' polar > north/east > lat/long (pairs must be complete).
#' `magdec` (degrees) is the declination applied to the polar azimuth and to the
#' N/E offsets (0 = true north); it must already be NA-free for records whose
#' format is polar or ne and whose north reference is magnetic.
#'
#' @param rec data.frame with columns ref_lat, ref_lon, in_east, in_north,
#'   in_distance, in_azimuth, in_lat, in_lon, magdec.
#' @return rec with added columns: input_format ("polar", "ne", "latlong" or
#'   NA), n_formats, east, north, distance, azimuth (true north, relative to the
#'   reference), lat, lon, format_disagreement_m (max distance between the
#'   position from the chosen format and positions from the other formats
#'   present; NA if only one format).
compute_positions <- function(rec) {
    has_pol <- !is.na(rec$in_distance) & !is.na(rec$in_azimuth)
    has_ne  <- !is.na(rec$in_east) & !is.na(rec$in_north)
    has_ll  <- !is.na(rec$in_lat) & !is.na(rec$in_lon)
    fmt <- ifelse(has_pol, "polar", ifelse(has_ne, "ne", ifelse(has_ll, "latlong", NA_character_)))
    D <- rec$magdec
    D[is.na(D)] <- 0

    # True-north offsets of every available format
    pol_ne <- polar_to_ne(rec$in_distance, rec$in_azimuth + D)
    ne_ne  <- rotate_magnetic_ne(rec$in_east, rec$in_north, D)   # identity where D == 0
    ll_ne  <- latlon_to_ne(rec$ref_lat, rec$ref_lon, rec$in_lat, rec$in_lon)

    pick <- function(col) ifelse(fmt == "polar", pol_ne[[col]],
                          ifelse(fmt == "ne", ne_ne[[col]], ll_ne[[col]]))
    rec$input_format <- fmt
    rec$n_formats <- has_pol + has_ne + has_ll
    rec$east  <- pick("east")
    rec$north <- pick("north")
    pol <- ne_to_polar(rec$east, rec$north)
    rec$distance <- pol$distance
    rec$azimuth  <- pol$azimuth

    ll <- ne_to_latlon(rec$ref_lat, rec$ref_lon, rec$east, rec$north)
    rec$lat <- ifelse(fmt %in% "latlong", rec$in_lat, ll$lat)
    rec$lon <- ifelse(fmt %in% "latlong", rec$in_lon, ll$lon)

    dis <- function(alt, present) ifelse(present, sqrt((alt$east - rec$east)^2 + (alt$north - rec$north)^2), NA)
    d_all <- cbind(dis(pol_ne, has_pol & fmt != "polar"), dis(ne_ne, has_ne & fmt != "ne"),
                   dis(ll_ne, has_ll & fmt != "latlong"))
    rec$format_disagreement_m <- suppressWarnings(apply(d_all, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)))
    # N/E vs polar only (both user-submitted; plot lat/long are earlier uploads, not inputs)
    rec$ne_polar_disagreement_m <- ifelse(has_pol & has_ne,
                                          sqrt((pol_ne$east - ne_ne$east)^2 + (pol_ne$north - ne_ne$north)^2), NA_real_)
    rec
}

# ---- Plot-area shapes (offsets in metres, true north) --------------------------

#' Corners of a rectangle centred on the origin.
#'
#' @param x_short length of the short side (PLOT_RECT_X), m.
#' @param y_long length of the long side (PLOT_RECT_Y), m.
#' @param dir_deg azimuth of the long side (PLOT_RECT_DIR), degrees from true north.
#' @return 4x2 matrix (east, north) of the corners, rows TR, BR, BL, TL, where
#'   "top" is the end of the long side the azimuth points to and "right" is 90
#'   degrees clockwise from it. With u = (sin a, cos a) along the long side and
#'   v = (cos a, -sin a) along the short side: corners = +-(Y/2) u +- (X/2) v.
#'   Azimuths a and a + 180 give the same rectangle.
rect_corners_ne <- function(x_short, y_long, dir_deg) {
    a <- dir_deg * pi / 180
    u <- c(sin(a), cos(a))
    v <- c(cos(a), -sin(a))
    hy <- y_long / 2; hx <- x_short / 2
    m <- rbind(TR =  hy * u + hx * v, BR = -hy * u + hx * v,
               BL = -hy * u - hx * v, TL =  hy * u - hx * v)
    colnames(m) <- c("east", "north")
    m
}

#' Vertices of a circle of radius r (m) centred on the origin.
circle_ne <- function(r, n = 64) {
    t <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
    cbind(east = r * sin(t), north = r * cos(t))
}

#' Vertices of a north-aligned square of side `side` (m) centred on the origin.
square_ne <- function(side) rect_corners_ne(side, side, 0)

#' Build a lon/lat polygon (sfg) from offsets around a reference point.
#'
#' @param ring_ne 2-column matrix of vertex offsets (east, north), not closed.
polygon_from_ne <- function(ref_lat, ref_lon, ring_ne) {
    ll <- ne_to_latlon(ref_lat, ref_lon, ring_ne[, 1], ring_ne[, 2])
    xy <- cbind(ll$lon, ll$lat)
    sf::st_polygon(list(rbind(xy, xy[1, ])))
}
