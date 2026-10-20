# run_coordconv.R -- command-line entry point of icosetc-coordconv.
#
# Converts ICOS ecosystem station coordinates (sampling plots and BADM groups)
# between lat/long, north/east offsets and polar coordinates, and writes
# coordinate tables and KML files.
#
# Usage:
#   Rscript run_coordconv.R --sites <codes|all> [--what plots,badm,features] [--kml]
#                           [--out <dir>] [--cache] [--offline <raw_dir>] [--settings <file>]
# Inputs : internal API (live) or a raw directory (--offline); settings file.
# Outputs: <out>/<run_id>/ with coordconv.log, summary.tsv, sites/, merged/ (and raw/ with --cache).

# Repository root: COORDCONV_ROOT if set by a calling script, else this script's folder.
.coordconv_root <- function() {
    if (exists("COORDCONV_ROOT", envir = globalenv())) return(normalizePath(get("COORDCONV_ROOT", envir = globalenv())))
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(f)) return(normalizePath(dirname(f[1])))
    normalizePath(".")
}
ROOT <- .coordconv_root()
for (.f in list.files(file.path(ROOT, "R"), pattern = "\\.R$", full.names = TRUE)) source(.f)
rm(.f)

USAGE <- "Usage: Rscript run_coordconv.R --sites <codes|all> [options]

  --sites <list|all>       comma-separated station codes, or 'all'
  --what <list>            any of plots,badm,features (default: plots,badm)
  --kml                    also write KML layers of plots and BADM points
  --out <dir>              output root (default: settings out_dir, under the repository)
  --cache                  keep raw downloads under <out>/<run_id>/raw/
  --offline <raw_dir>      read raw files instead of the API (same layout as raw/)
  --settings <file>        settings file (default: COORDCONV_SETTINGS or config/settings_local.R)
  --previous <dir>         merged files of a previous run (GRP_<X>.tsv), used for the LARGE_SHIFT
                           check where the API returns no previous calculated values
  --help                   show this help
"

#' Parse command-line arguments (--key value, --key=value, or flags).
#'
#' @param args character vector (commandArgs(TRUE)).
#' @return list of options.
parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
    flags <- c("kml", "cache", "help")
    opts <- list(sites = NULL, what = "plots,badm", kml = FALSE, out = NULL, cache = FALSE,
                 offline = NULL, settings = NULL, previous = NULL, help = FALSE)
    i <- 1
    while (i <= length(args)) {
        a <- args[i]
        if (!grepl("^--", a)) stop("Unexpected argument: ", a, call. = FALSE)
        kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
        key <- kv[1]
        if (!key %in% names(opts)) stop("Unknown option --", key, call. = FALSE)
        if (key %in% flags) {
            opts[[key]] <- TRUE
        } else if (length(kv) > 1) {
            opts[[key]] <- paste(kv[-1], collapse = "=")
        } else {
            if (i == length(args)) stop("Option --", key, " needs a value", call. = FALSE)
            i <- i + 1
            opts[[key]] <- args[i]
        }
        i <- i + 1
    }
    opts$what <- trimws(strsplit(opts$what, ",", fixed = TRUE)[[1]])
    bad <- setdiff(opts$what, c("plots", "badm", "features"))
    if (length(bad)) stop("Unknown --what value(s): ", paste(bad, collapse = ", "), call. = FALSE)
    if (!is.null(opts$sites)) opts$sites <- trimws(strsplit(opts$sites, ",", fixed = TRUE)[[1]])
    opts
}

#' Process one station: download once, convert, write tables and KML.
#'
#' @return list(updates = named list of update tables, counts = data.frame).
process_station <- function(site, info, src, opts, site_dir, raw_dir = NULL) {
    eco <- recode_eco(info$igbp)
    plot_df <- tryCatch(fetch_plot_table(src, site), error = function(e) {
        log_warn("plot table: ", conditionMessage(e)); data.frame() })
    location <- fetch_badm(src, site, "GRP_LOCATION")
    towers <- tower_positions(location)
    if (!nrow(towers)) log_warn("GRP_LOCATION: no tower position")
    else if (nrow(towers) > 1) log_info("GRP_LOCATION: ", nrow(towers), " tower positions over time")

    results <- list(); counts <- list()
    plots <- tryCatch(process_plots(plot_df, towers, eco, site), error = function(e) {
        log_error("GRP_PLOT: ", conditionMessage(e))
        counts$GRP_PLOT <<- transform(.counts(site, "GRP_PLOT", .rec_template(0), 0),
                                      status = "error", message = conditionMessage(e))
        NULL })
    if (!is.null(plots)) towers <- plots$towers
    plot_points <- if (is.null(plots)) data.frame(plot_id = character(0), lat = numeric(0), lon = numeric(0)) else plots$points
    if ("plots" %in% opts$what && !is.null(plots)) {
        results$GRP_PLOT <- plots; counts$GRP_PLOT <- plots$counts
    }

    badm_res <- NULL
    if ("badm" %in% opts$what) {
        badm <- list(GRP_LOCATION = location)
        for (g in CONVERTED_GROUPS) {
            badm[[g]] <- tryCatch(fetch_badm(src, site, g), error = function(e) {
                log_warn(g, ": ", conditionMessage(e))
                counts[[g]] <<- transform(.counts(site, g, .rec_template(0), 0), status = "error",
                                          message = conditionMessage(e))
                NULL })
        }
        badm_res <- process_badm(badm, towers, plot_points, site)
        for (g in names(badm_res)) { results[[g]] <- badm_res[[g]]; counts[[g]] <- badm_res[[g]]$counts }
    }

    features <- NULL
    if ("features" %in% opts$what) {
        if (identical(info$class, "Associated")) {
            log_info("features: Associated station, polygon features not processed")
        } else {
            # Polygon features are optional: a failure is logged and never discards the
            # station's plot and BADM outputs.
            exdir <- if (!is.null(raw_dir)) file.path(raw_dir, site, "features") else tempfile("features_")
            get_files <- function(type) tryCatch(fetch_feature_files(src, site, type, file.path(exdir, type)),
                                                 error = function(e) {
                                                     log_warn("features (", type, "): ", conditionMessage(e))
                                                     character(0) })
            files <- c(get_files("kml"), get_files("shape"))
            if (!length(files)) {
                log_warn("features: no polygon files available")
            } else {
                tw <- tower_at(towers, as.Date(NA))
                features <- tryCatch(prepare_features(read_feature_layers(files), tw$ref_lat, tw$ref_lon),
                                     error = function(e) { log_warn("features: ", conditionMessage(e)); NULL })
                if (!is.null(features))
                    log_info("features: TA ", if (is.null(features$TA)) "missing" else "ok", ", EA ",
                             if (is.null(features$EA)) "none" else "ok", ", LCT ", length(features$LCT))
            }
        }
    }

    # QC flags (records are never withheld; flags are for human review)
    tw <- tower_at(towers, as.Date(NA))
    flags <- if (is.na(tw$ref_lat)) NULL else
        qc_station_flags(site, results, opts$qc %||% QC_DEFAULTS, list(lat = tw$ref_lat, lon = tw$ref_lon),
                         ta = features$TA, previous = opts$previous_values)
    if (!is.null(flags) && nrow(flags)) {
        for (fl in intersect(QC_FLAGS, unique(flags$flag))) {
            ff <- flags[flags$flag == fl, ]
            log_warn("QC ", fl, ": ", nrow(ff), " record(s) flagged (",
                     paste(sprintf("%s %d", names(table(ff$group)), as.integer(table(ff$group))), collapse = ", "), ")")
        }
        write_tsv(flags, file.path(site_dir, sprintf("%s-QC_flags.tsv", site)))
    }
    if (length(counts)) {
        for (g in names(counts)) for (fl in QC_FLAGS)
            counts[[g]][[paste0("n_flag_", fl)]] <- if (is.null(flags)) 0L else sum(flags$group == g & flags$flag == fl)
    }

    write_site_tables(site, results, site_dir)
    write_site_kml(site, site_dir, plots = plots, badm = badm_res, features = features,
                   plot_kml = opts$kml && "plots" %in% opts$what, badm_kml = opts$kml && "badm" %in% opts$what)
    list(updates = lapply(results, `[[`, "update"),
         counts = if (length(counts)) do.call(rbind, counts) else NULL,
         flags = flags)
}

#' Run the conversion.
#'
#' @param opts options from parse_cli().
#' @return invisible list(run_dir, summary, n_failed).
main <- function(opts, echo = TRUE) {
    if (isTRUE(opts$help)) { cat(USAGE); return(invisible(NULL)) }
    if (is.null(opts$sites)) stop("--sites is required\n", USAGE, call. = FALSE)
    check_packages()
    settings <- load_settings(opts$settings, root = ROOT, required = is.null(opts$offline))
    opts$qc <- qc_settings(settings)
    opts$previous_values <- read_previous_values(opts$previous)

    out_root <- opts$out %||% file.path(ROOT, settings$out_dir)
    run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
    run_dir <- file.path(out_root, run_id)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    logger_init(file.path(run_dir, "coordconv.log"), echo = echo)
    log_info("icosetc-coordconv run ", run_id, " | what: ", paste(opts$what, collapse = ","),
             " | kml: ", opts$kml, " | source: ", if (is.null(opts$offline)) "live" else "offline",
             " | cache: ", opts$cache)
    writeLines(c(paste("run_id:", run_id), paste("sites:", paste(opts$sites, collapse = ",")),
                 paste("what:", paste(opts$what, collapse = ",")), paste("kml:", opts$kml),
                 paste("source:", if (is.null(opts$offline)) "live" else "offline"), paste("cache:", opts$cache),
                 paste("oce available:", oce_available())),
               file.path(run_dir, "run_settings.txt"))

    raw_dir <- NULL
    src <- if (!is.null(opts$offline)) {
        raw_dir <- normalizePath(opts$offline)
        make_source(settings, "offline", raw_dir = raw_dir)
    } else {
        if (opts$cache) raw_dir <- file.path(run_dir, "raw")
        make_source(settings, "live", raw_dir = raw_dir, cache = opts$cache)
    }

    stations <- fetch_station_list(src)
    if (identical(tolower(opts$sites), "all")) {
        sites <- setdiff(stations$site_code, settings$excluded_sites)
    } else {
        sites <- opts$sites
        unknown <- setdiff(sites, stations$site_code)
        if (length(unknown)) log_warn("station(s) not in the station list: ", paste(unknown, collapse = ", "))
    }
    log_info(length(sites), " station(s) to process")

    summary <- list(); updates <- list(); station_rows <- list(); all_flags <- list()
    for (i in seq_along(sites)) {
        s <- sites[i]
        log_set_site(s)
        log_info("processing (", i, " of ", length(sites), ")")
        info <- stations[match(s, stations$site_code), ]
        if (is.na(info$site_code)) info <- data.frame(site_code = s, class = NA, igbp = NA)
        t0 <- Sys.time()
        res <- tryCatch(withCallingHandlers(
                            process_station(s, info, src, opts, file.path(run_dir, "sites", s), raw_dir),
                            warning = function(w) { log_warn("R warning: ", conditionMessage(w))
                                                    invokeRestart("muffleWarning") }),
                        error = function(e) { log_error("station failed: ", conditionMessage(e)); e })
        failed <- inherits(res, "error")
        if (!failed) { summary[[s]] <- res$counts; updates[[s]] <- res$updates; all_flags[[s]] <- res$flags }
        grp_err <- !failed && !is.null(res$counts) && any(res$counts$status == "error")
        station_rows[[s]] <- data.frame(
            site = s, status = if (failed) "error" else if (grp_err) "partial" else "ok",
            n_warn = log_count(s, "WARN"), n_error = log_count(s, "ERROR"),
            seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
            message = if (failed) conditionMessage(res) else "", stringsAsFactors = FALSE)
        log_set_site(NULL)
    }

    merge_group_files(updates, file.path(run_dir, "merged"))
    flags <- do.call(rbind, all_flags)
    if (is.null(flags)) flags <- data.frame(site = character(0), group = character(0), ID = character(0),
                                            flag = character(0), detail = character(0))
    write_tsv(flags, file.path(run_dir, "merged", "QC_flags.tsv"))
    log_info("QC flags: ", if (nrow(flags)) paste(sprintf("%s %d", names(table(flags$flag)), as.integer(table(flags$flag))), collapse = ", ") else "none",
             " (merged/QC_flags.tsv)")
    counts <- if (length(summary)) do.call(rbind, summary) else NULL
    if (!is.null(counts)) write_summary(counts, file.path(run_dir, "summary.tsv"))
    st <- do.call(rbind, station_rows)
    write_summary(st, file.path(run_dir, "stations.tsv"))

    log_info("summary: ", sum(st$status == "ok"), " ok, ", sum(st$status == "partial"), " partial, ",
             sum(st$status == "error"), " failed")
    w <- max(nchar(st$site), 4)
    log_info(sprintf("%-*s  %-7s  %6s  %6s  %s", w, "site", "status", "warn", "error", "message"))
    for (k in seq_len(nrow(st)))
        log_info(sprintf("%-*s  %-7s  %6d  %6d  %s", w, st$site[k], st$status[k], st$n_warn[k],
                         st$n_error[k], st$message[k]))
    log_info("outputs in ", run_dir)
    invisible(list(run_dir = run_dir, summary = counts, stations = st, n_failed = sum(st$status == "error")))
}

if (sys.nframe() == 0L) {
    res <- tryCatch(main(parse_cli()), error = function(e) {
        message("ERROR: ", conditionMessage(e)); quit(status = 2) })
    if (!is.null(res) && res$n_failed > 0) quit(status = 1)
}
