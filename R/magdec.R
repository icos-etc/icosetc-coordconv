# magdec.R -- magnetic declination and magnetic-north corrections.
#
# Convention: declination D is east-positive (IGRF via oce::magneticField).
#   true azimuth = magnetic azimuth + D
#   true offsets: E_t = E cos D + N sin D ;  N_t = N cos D - E sin D
# 'oce' is loaded lazily: everything else runs without it; only records with a
# magnetic north reference need it.

IGRF_VERSION    <- 14
IGRF_VALID_FROM <- as.Date("1900-01-01")
IGRF_VALID_TO   <- as.Date("2030-12-31")

#' Is a north-reference value "Magnetic" (case-insensitive; NA -> FALSE)?
is_magnetic <- function(northref) {
    !is.na(northref) & tolower(trimws(northref)) == "magnetic"
}

#' Magnetic declination (degrees, east-positive) from IGRF-14.
#'
#' Evaluated at 12:00 UTC of each date. Dates before 1900 are an error; dates
#' after the IGRF-14 validity (2030-12-31) raise a warning (extrapolated).
#'
#' @param lat,lon location(s), decimal degrees (recycled).
#' @param date Date vector.
#' @return numeric vector (NA where date or location is NA).
declination <- function(lat, lon, date) {
    if (!oce_available()) {
        stop("Magnetic declination needs the R package 'oce' (>= ", OCE_MIN_VERSION, "), ",
             "which is not available. Install it (install.packages(\"oce\")) to process ",
             "records with a magnetic north reference.", call. = FALSE)
    }
    n <- max(length(lat), length(lon), length(date))
    lat <- rep_len(lat, n); lon <- rep_len(lon, n); date <- rep_len(as.Date(date), n)
    if (any(date < IGRF_VALID_FROM, na.rm = TRUE))
        stop("Declination requested for a date before ", IGRF_VALID_FROM, call. = FALSE)
    if (any(date > IGRF_VALID_TO, na.rm = TRUE))
        warning("Declination extrapolated beyond IGRF-14 validity (", IGRF_VALID_TO, ")", call. = FALSE)
    out <- rep(NA_real_, n)
    ok <- !is.na(lat) & !is.na(lon) & !is.na(date)
    if (any(ok)) {
        tm <- as.POSIXct(paste(format(date[ok]), "12:00:00"), tz = "UTC")
        out[ok] <- oce::magneticField(longitude = lon[ok], latitude = lat[ok], time = tm,
                                      version = IGRF_VERSION)$declination
    }
    out
}

#' Rotate north/east offsets measured from magnetic north to true north.
#'
#' @param east,north offsets (m) measured with magnetic north.
#' @param D declination (degrees, east-positive).
#' @return data.frame(east, north) relative to true north.
rotate_magnetic_ne <- function(east, north, D) {
    d <- D * pi / 180
    data.frame(east  = east * cos(d) + north * sin(d),
               north = north * cos(d) - east * sin(d))
}

#' Correct a magnetic azimuth to true north: (az + D) mod 360.
correct_magnetic_azimuth <- function(azimuth, D) {
    (azimuth + D) %% 360
}
