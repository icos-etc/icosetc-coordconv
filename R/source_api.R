# source_api.R -- adapter for the internal API (all network and raw-file I/O).
#
# Inputs : settings (base URL, endpoint templates), station codes.
# Outputs: data.frames (all columns as character) and paths of downloaded zips.
#
# A "source" is either live (HTTP GET, TLS verification on by default) or
# offline (a raw directory with the cache layout). With cache = TRUE every live
# download is kept under <raw_dir> so a run can be reproduced with --offline.
# Raw layout:
#   <raw_dir>/_network/station_class.tsv, station_igbp.tsv
#   <raw_dir>/<SITE>/plot.tsv, <SITE>/<GRP_X>.tsv, <SITE>/features_<kml|shape>.zip
# Error messages never contain the API host (it is replaced by "<api>").

#' Create a data source.
#'
#' @param settings settings list (see config/settings_example.R).
#' @param mode "live" or "offline".
#' @param raw_dir cache directory (live, cache = TRUE) or input directory (offline).
#' @param cache keep live downloads under raw_dir.
#' @param pause seconds to wait after each live request.
make_source <- function(settings, mode = c("live", "offline"), raw_dir = NULL, cache = FALSE, pause = 0) {
    mode <- match.arg(mode)
    if (mode == "offline" && (is.null(raw_dir) || !dir.exists(raw_dir)))
        stop("Offline raw directory not found: ", raw_dir %||% "<none>", call. = FALSE)
    if (mode == "live" && (is.na(settings$api_base_url) || !nzchar(settings$api_base_url)))
        stop("settings$api_base_url is not set (needed for live downloads)", call. = FALSE)
    if (mode == "live" && cache && is.null(raw_dir)) stop("cache = TRUE needs raw_dir", call. = FALSE)
    list(settings = settings, mode = mode, raw_dir = raw_dir, cache = cache, pause = pause)
}

# Relative raw path of a resource.
.raw_rel <- function(site, kind, var = NULL) {
    switch(kind,
           network = file.path("_network", paste0("station_", var, ".tsv")),
           plot    = file.path(site, "plot.tsv"),
           badm    = file.path(site, paste0(var, ".tsv")),
           feature = file.path(site, paste0("features_", var, ".zip")),
           stop("unknown resource kind: ", kind))
}

# Build the URL of a resource from the endpoint templates.
.api_url <- function(settings, site, kind, var = NULL) {
    fill <- function(tpl, site, var = "", type = "")
        gsub("{type}", type, gsub("{var}", var, gsub("{site}", site, tpl, fixed = TRUE), fixed = TRUE), fixed = TRUE)
    base <- sub("/+$", "", settings$api_base_url)
    path <- switch(kind,
                   network = fill(settings$endpoints$badm, "select_all", var = settings$station_list_vars[[var]]),
                   plot    = fill(settings$endpoints$plot, tolower(site), type = "plot"),
                   badm    = fill(settings$endpoints$badm, tolower(site), var = var),
                   feature = fill(settings$endpoints$plot, tolower(site), type = var))
    paste0(base, "/", path)
}

# Download url to dest; on failure stop with a sanitised message.
.download <- function(url, dest, settings) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    sanitise <- function(m) gsub(settings$api_base_url, "<api>", m, fixed = TRUE)
    res <- tryCatch({
        if (isTRUE(settings$tls_verify)) {
            suppressWarnings(utils::download.file(url, dest, method = "libcurl", quiet = TRUE, mode = "wb"))
        } else {
            if (!requireNamespace("curl", quietly = TRUE))
                stop("tls_verify = FALSE needs the 'curl' package")
            h <- curl::new_handle(ssl_verifypeer = 0L, ssl_verifyhost = 0L)
            curl::curl_download(url, dest, handle = h, quiet = TRUE)
            0L
        }
    }, error = function(e) structure(sanitise(conditionMessage(e)), class = "dl_error"))
    if (inherits(res, "dl_error") || !identical(as.integer(res), 0L)) {
        if (file.exists(dest)) unlink(dest)
        msg <- if (inherits(res, "dl_error")) unclass(res) else paste("status", res)
        if (grepl("certificate|SSL", msg, ignore.case = TRUE))
            msg <- paste0(msg, " (TLS verification failed; see tls_verify in config/settings_example.R)")
        stop("download failed: ", msg, call. = FALSE)
    }
    invisible(dest)
}

#' Get the local path of a raw resource, downloading it in live mode.
#'
#' @return path, or NULL when an offline file does not exist.
fetch_raw <- function(src, site, kind, var = NULL) {
    rel <- .raw_rel(site, kind, var)
    if (src$mode == "offline") {
        p <- file.path(src$raw_dir, rel)
        return(if (file.exists(p)) p else NULL)
    }
    ext <- if (kind == "feature") ".zip" else ".tsv"
    dest <- if (src$cache) file.path(src$raw_dir, rel) else tempfile(fileext = ext)
    on.exit(if (src$pause > 0) Sys.sleep(src$pause))
    .download(.api_url(src$settings, site, kind, var), dest, src$settings)
}

#' Read an API tab-separated table; every column as character.
#'
#' Parsed explicitly (no quoting): rows with fewer fields than the header are
#' padded; extra trailing empty fields are dropped; rows with extra non-empty
#' fields are misaligned (e.g. a tab inside a comment) and are dropped with a
#' warning, since their values cannot be assigned to columns reliably.
#' Empty files, or files with a header only, give a data.frame with 0 rows.
#'
#' @param path file path (NULL allowed).
#' @param label name used in the warning.
read_api_tsv <- function(path, label = basename(path %||% "")) {
    if (is.null(path) || !file.exists(path) || file.info(path)$size == 0) return(data.frame())
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8", skipNul = TRUE)
    lines <- sub("\r$", "", lines)
    lines <- lines[nzchar(trimws(lines))]
    if (!length(lines)) return(data.frame())
    hdr <- trimws(strsplit(lines[1], "\t", fixed = TRUE)[[1]])
    nc <- length(hdr)
    fields <- strsplit(lines[-1], "\t", fixed = TRUE)
    nf <- lengths(fields)
    extra <- which(nf > nc)
    bad <- extra[vapply(fields[extra], function(f) any(nzchar(trimws(f[(nc + 1):length(f)]))), logical(1))]
    if (length(bad)) {
        log_warn(label, ": ", length(bad), " row(s) with more fields than the header (misaligned) dropped")
        fields <- fields[-bad]
    }
    m <- vapply(fields, function(f) { length(f) <- nc; trimws(f) }, character(nc))
    m <- matrix(m, ncol = nc, byrow = TRUE)
    m[!is.na(m) & (m == "" | m == "NA")] <- NA
    df <- as.data.frame(m, stringsAsFactors = FALSE)
    names(df) <- hdr
    df
}

#' Network station list: one row per station with code, class and IGBP code.
#'
#' The latest DATE_LABEL / record per station is used.
#' @return data.frame(site_code, class, igbp).
fetch_station_list <- function(src) {
    cls <- read_api_tsv(fetch_raw(src, NULL, "network", "class"))
    eco <- read_api_tsv(fetch_raw(src, NULL, "network", "igbp"))
    need <- function(df, cols, what) if (!all(cols %in% names(df)))
        stop("station list (", what, ") lacks column(s): ", paste(setdiff(cols, names(df)), collapse = ", "), call. = FALSE)
    need(cls, c("site_code", "CLASS_ICOS"), "class"); need(eco, c("site_code", "IGBP"), "igbp")
    if ("DATE_LABEL" %in% names(cls)) {
        cls <- cls[!is.na(cls$DATE_LABEL), ]
        cls <- cls[order(cls$site_code, cls$DATE_LABEL), ]
    }
    cls <- cls[!duplicated(cls$site_code, fromLast = TRUE), ]
    eco <- eco[!is.na(eco$IGBP) & nzchar(eco$IGBP), ]
    eco <- eco[!duplicated(eco$site_code, fromLast = TRUE), ]
    out <- data.frame(site_code = cls$site_code, class = cls$CLASS_ICOS, stringsAsFactors = FALSE)
    out$igbp <- eco$IGBP[match(out$site_code, eco$site_code)]
    out[order(out$site_code), ]
}

#' Sampling-plot table of a station (downloaded once, shared by all processing).
#'
#' Normalises the identifier/coordinate headers to PLOT_ID, PLOT_LOCATION_LAT,
#' PLOT_LOCATION_LONG and drops the "File download success" marker row.
fetch_plot_table <- function(src, site) {
    df <- read_api_tsv(fetch_raw(src, site, "plot"), "plot table")
    if (!nrow(df)) return(df)
    nm <- names(df)
    ren <- function(pattern, to) {
        i <- which(grepl(pattern, nm, ignore.case = TRUE) & !nm %in% c("PLOT_ID", "PLOT_LOCATION_LAT", "PLOT_LOCATION_LONG"))
        if (!to %in% nm && length(i) == 1) nm[i] <<- to
    }
    ren("FEAT", "PLOT_ID"); ren("LAT$", "PLOT_LOCATION_LAT"); ren("LON(G)?$", "PLOT_LOCATION_LONG")
    names(df) <- nm
    marker <- apply(df, 1, function(r) any(grepl("File download success", r, fixed = TRUE)))
    df[!marker, , drop = FALSE]
}

#' One BADM variable group of a station (0-row data.frame if empty / absent).
fetch_badm <- function(src, site, group) {
    read_api_tsv(fetch_raw(src, site, "badm", group), group)
}

#' Download and unzip the polygon features of a station.
#'
#' @param type "kml" or "shape".
#' @param exdir directory to unzip into.
#' @return unzipped file paths (character(0) if nothing available).
fetch_feature_files <- function(src, site, type, exdir) {
    zip <- fetch_raw(src, site, "feature", type)
    if (is.null(zip) || !file.exists(zip) || file.info(zip)$size <= 100) return(character(0))
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    files <- tryCatch(utils::unzip(zip, exdir = exdir, overwrite = TRUE), error = function(e) character(0))
    if (!src$cache && src$mode == "live") unlink(zip)
    files
}

#' Upload calculated coordinates to the database -- NOT IMPLEMENTED (out of scope v1).
#'
#' Intended design: DBI + odbc connection built from settings$db_* (password
#' from an environment variable), one transaction per BADM group, updating
#' <X>_LOCATION_CALC_LAT/LONG by record ID from the merged update files.
upload_coordinates <- function(...) {
    stop("upload_coordinates(): not implemented (DB upload is out of scope for v1)", call. = FALSE)
}
