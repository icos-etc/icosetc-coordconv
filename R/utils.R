# utils.R -- logging, date parsing, ecosystem recoding and small helpers.
#
# Inputs : messages, raw BADM date strings, IGBP codes.
# Outputs: log lines (console + optional log file), Date vectors, eco labels.

# ---- Logging ----------------------------------------------------------------

.log_state <- new.env()
.log_state$file  <- NULL
.log_state$site  <- NULL
.log_state$echo  <- TRUE
.log_state$count <- list()

#' Initialise the logger.
#'
#' @param file path of the log file (appended), or NULL for console only.
#' @param echo also print to the console.
logger_init <- function(file = NULL, echo = TRUE) {
    .log_state$file  <- file
    .log_state$echo  <- echo
    .log_state$site  <- NULL
    .log_state$count <- list()
    if (!is.null(file)) dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    invisible(TRUE)
}

#' Set (or clear with NULL) the station code prefixed to log lines.
log_set_site <- function(site = NULL) {
    .log_state$site <- site
    invisible(site)
}

#' Write one log line.
#'
#' @param level "INFO", "WARN" or "ERROR".
#' @param ... pieces pasted together (no separator).
log_msg <- function(level, ...) {
    site <- .log_state$site
    line <- sprintf("%s [%-5s] %s%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level,
                    if (is.null(site)) "" else paste0(site, ": "), paste0(..., collapse = ""))
    key <- if (is.null(site)) "_run" else site
    cnt <- .log_state$count[[key]]
    if (is.null(cnt)) cnt <- c(INFO = 0L, WARN = 0L, ERROR = 0L)
    cnt[[level]] <- cnt[[level]] + 1L
    .log_state$count[[key]] <- cnt
    if (isTRUE(.log_state$echo)) message(line)
    if (!is.null(.log_state$file)) cat(line, "\n", sep = "", file = .log_state$file, append = TRUE)
    invisible(line)
}
log_info  <- function(...) log_msg("INFO", ...)
log_warn  <- function(...) log_msg("WARN", ...)
log_error <- function(...) log_msg("ERROR", ...)

#' Number of log lines of a level written for a station (or the run, "_run").
log_count <- function(key, level) {
    cnt <- .log_state$count[[key]]
    if (is.null(cnt)) 0L else cnt[[level]]
}

# ---- Dates ------------------------------------------------------------------

#' Parse BADM dates.
#'
#' YYYY -> YYYY-07-01, YYYYMM -> YYYY-MM-15, YYYYMMDD[HHMM[SS]] -> that day.
#' Non-digit characters are removed first. Invalid or empty values give NA.
#'
#' @param x character or numeric vector.
#' @return Date vector.
parse_badm_date <- function(x) {
    if (is.numeric(x)) x <- format(x, scientific = FALSE, trim = TRUE)
    x <- gsub("[^0-9]", "", as.character(x))
    out <- rep(as.Date(NA), length(x))
    n <- nchar(x)
    ok4 <- !is.na(x) & n == 4
    ok6 <- !is.na(x) & n == 6
    ok8 <- !is.na(x) & n >= 8
    out[ok4] <- as.Date(paste0(x[ok4], "0701"), "%Y%m%d")
    out[ok6] <- as.Date(paste0(x[ok6], "15"), "%Y%m%d")
    out[ok8] <- as.Date(substr(x[ok8], 1, 8), "%Y%m%d")
    out
}

# ---- Ecosystem types ----------------------------------------------------------

#' The single place where IGBP codes are mapped to the eco labels used for
#' plot areas: forest, crop, grassland, mire, urban, water.
ECO_TYPES <- c(
    DBF = "forest", EBF = "forest", ENF = "forest", MF = "forest", DNF = "forest",
    CRO = "crop",
    GRA = "grassland", CSH = "grassland", OSH = "grassland", SAV = "grassland",
    WET = "mire",
    URB = "urban",
    WAT = "water"
)

#' Recode IGBP codes to eco labels; unknown (non-NA) codes are an error.
#'
#' @param igbp character vector of IGBP codes.
#' @return character vector of eco labels (NA where igbp is NA/empty).
recode_eco <- function(igbp) {
    code <- toupper(trimws(as.character(igbp)))
    code[!nzchar(code)] <- NA
    unknown <- unique(code[!is.na(code) & !code %in% names(ECO_TYPES)])
    if (length(unknown))
        stop("Unknown IGBP code(s): ", paste(unknown, collapse = ", "),
             ". Add them to ECO_TYPES in R/utils.R.", call. = FALSE)
    unname(ECO_TYPES[code])
}

# ---- Small helpers ------------------------------------------------------------

#' Convert a character column to numeric ("" and "NA" -> NA, ',' decimal accepted).
as_num <- function(x) {
    if (is.null(x)) return(NULL)
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", "-9999")] <- NA
    suppressWarnings(as.numeric(sub(",", ".", x, fixed = TRUE)))
}

#' Get a column as character, or NA if the column does not exist.
col_or_na <- function(df, col) {
    if (!is.null(col) && col %in% names(df)) as.character(df[[col]]) else rep(NA_character_, nrow(df))
}

#' Row-wise first non-NA value over several columns (missing columns ignored).
coalesce_cols <- function(df, cols) {
    out <- rep(NA_character_, nrow(df))
    for (cl in cols) {
        v <- col_or_na(df, cl)
        v[!is.na(v) & !nzchar(trimws(v))] <- NA
        out[is.na(out)] <- v[is.na(out)]
    }
    out
}
