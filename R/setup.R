# setup.R -- package checks and settings loading.
#
# Inputs : installed packages; a settings file (R script defining `settings`).
# Outputs: errors with actionable messages; a validated settings list.
# Packages are never installed at runtime.

REQUIRED_PACKAGES <- c(sf = "1.0-0", data.table = "1.14.0")
OCE_MIN_VERSION   <- "1.8-4"   # first version with IGRF-14 as default model

#' Check that required packages are installed with a minimum version.
#'
#' @param req named character vector: package = minimum version.
#' @return invisible TRUE; stops listing every missing/outdated package.
check_packages <- function(req = REQUIRED_PACKAGES) {
    problems <- character(0)
    for (p in names(req)) {
        if (!requireNamespace(p, quietly = TRUE)) {
            problems <- c(problems, sprintf("  - %s (>= %s): not installed", p, req[[p]]))
        } else if (utils::packageVersion(p) < req[[p]]) {
            problems <- c(problems, sprintf("  - %s (>= %s): found %s", p, req[[p]],
                                            as.character(utils::packageVersion(p))))
        }
    }
    if (length(problems)) {
        stop("Required R packages are missing or too old:\n", paste(problems, collapse = "\n"),
             "\nInstall them (e.g. install.packages(c(", paste0('"', names(req), '"', collapse = ", "),
             "))) and re-run.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Is a suitable 'oce' available (used only for magnetic declination)?
#'
#' The option coordconv.force_no_oce = TRUE simulates its absence (tests).
#' @return logical(1).
oce_available <- function(min_version = OCE_MIN_VERSION) {
    if (isTRUE(getOption("coordconv.force_no_oce"))) return(FALSE)
    requireNamespace("oce", quietly = TRUE) && utils::packageVersion("oce") >= min_version
}

# Defaults used when no settings file is needed (offline runs).
DEFAULT_SETTINGS <- list(
    api_base_url = NA_character_, endpoints = list(), station_list_vars = c(class = NA, igbp = NA),
    excluded_sites = character(0), out_dir = "output", tls_verify = TRUE,
    qc = list(format_disagreement_max_m = 1.0, large_shift_m = 20)
)

SETTINGS_REQUIRED_KEYS <- c("api_base_url", "endpoints", "station_list_vars",
                            "excluded_sites", "out_dir", "tls_verify")

#' Locate and load the settings file.
#'
#' Lookup order: `path` argument (--settings), environment variable
#' COORDCONV_SETTINGS, <root>/config/settings_local.R.
#'
#' @param path optional explicit path.
#' @param root repository root directory.
#' @param required if FALSE and no file is found, DEFAULT_SETTINGS are returned.
#' @return settings list.
load_settings <- function(path = NULL, root = ".", required = TRUE) {
    candidates <- c(path, Sys.getenv("COORDCONV_SETTINGS", unset = ""),
                    file.path(root, "config", "settings_local.R"))
    candidates <- candidates[nzchar(candidates)]
    found <- candidates[file.exists(candidates)][1]
    if (is.na(found)) {
        if (!required) return(DEFAULT_SETTINGS)
        stop("No settings file found. Copy ", file.path(root, "config", "settings_example.R"),
             " to config/settings_local.R and fill in the values, or pass --settings <file> ",
             "or set COORDCONV_SETTINGS.", call. = FALSE)
    }
    env <- new.env(parent = baseenv())
    sys.source(found, envir = env)
    if (!exists("settings", envir = env, inherits = FALSE))
        stop("Settings file does not define a list called `settings`: ", found, call. = FALSE)
    s <- get("settings", envir = env)
    missing_keys <- setdiff(SETTINGS_REQUIRED_KEYS, names(s))
    if (length(missing_keys))
        stop("Settings file is missing key(s): ", paste(missing_keys, collapse = ", "),
             " (see config/settings_example.R)", call. = FALSE)
    if (!all(c("badm", "plot") %in% names(s$endpoints)))
        stop("settings$endpoints must define 'badm' and 'plot'", call. = FALSE)
    if (!all(c("class", "igbp") %in% names(s$station_list_vars)))
        stop("settings$station_list_vars must be named c(class = ..., igbp = ...)", call. = FALSE)
    s
}
