# run_example.R -- run the synthetic example offline and compare with the expected outputs.
#
# Usage: Rscript examples/run_example.R
# Inputs : examples/raw/ (invented station XX-Tst, invented record IDs)
# Outputs: a temporary output directory; comparison with examples/expected/*.tsv
#          (numeric tolerance 1e-9 after rounding). Needs no network and no settings file.
# The example contains magnetic records, so it needs 'oce'; without it the run
# must report the explicit declination error and the comparison is skipped.

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grepl("^--file=", .args)])
COORDCONV_ROOT <- normalizePath(file.path(dirname(if (length(.file)) .file else "examples/run_example.R"), ".."))
source(file.path(COORDCONV_ROOT, "run_coordconv.R"))

#' Run the example and compare it with examples/expected.
#'
#' @param out output root (default: a temporary directory).
#' @param echo print log lines on the console.
#' @return list(status = "PASS"|"FAIL"|"SKIP", details = character()).
run_example <- function(out = tempfile("coordconv_example_"), echo = FALSE) {
    ex <- file.path(COORDCONV_ROOT, "examples")
    args <- c("--sites", "XX-Tst", "--what", "plots,badm", "--kml", "--offline", file.path(ex, "raw"), "--out", out)
    opts <- parse_cli(args)
    res <- main(opts, echo = echo)
    logger_init(NULL, echo = echo)
    log_lines <- readLines(file.path(res$run_dir, "coordconv.log"))
    if (!oce_available()) {
        ok <- any(grepl("needs the R package 'oce'", log_lines, fixed = TRUE))
        return(list(status = if (ok) "SKIP" else "FAIL",
                    details = if (ok) "oce not available: explicit declination error reported, comparison skipped"
                              else "oce not available and no explicit declination error in the log"))
    }
    site_dir <- file.path(res$run_dir, "sites", "XX-Tst")
    expected <- list.files(file.path(ex, "expected"), pattern = "\\.tsv$")
    produced <- list.files(site_dir, pattern = "\\.tsv$")
    details <- character(0)
    miss <- setdiff(expected, produced); extra <- setdiff(produced, expected)
    if (length(miss)) details <- c(details, paste("missing:", paste(miss, collapse = ", ")))
    if (length(extra)) details <- c(details, paste("unexpected:", paste(extra, collapse = ", ")))
    for (f in intersect(expected, produced)) {
        a <- utils::read.delim(file.path(ex, "expected", f), colClasses = "character", na.strings = "")
        b <- utils::read.delim(file.path(site_dir, f), colClasses = "character", na.strings = "")
        if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) {
            details <- c(details, paste0(f, ": different shape or columns")); next
        }
        for (cl in names(a)) {
            na <- suppressWarnings(as.numeric(a[[cl]])); nb <- suppressWarnings(as.numeric(b[[cl]]))
            numeric_col <- all(is.na(na) == is.na(a[[cl]])) && all(is.na(nb) == is.na(b[[cl]]))
            same <- if (numeric_col) isTRUE(all.equal(na, nb, tolerance = 1e-9)) else identical(a[[cl]], b[[cl]])
            if (!same) details <- c(details, paste0(f, ": column ", cl, " differs"))
        }
    }
    kml <- list.files(file.path(site_dir, "kml"), pattern = "\\.kml$")
    if (length(kml) < 5) details <- c(details, "fewer KML layers than expected")
    list(status = if (length(details)) "FAIL" else "PASS",
         details = c(details, sprintf("%d table(s) compared, %d KML layer(s)", length(intersect(expected, produced)), length(kml))))
}

if (sys.nframe() == 0L) {
    r <- run_example()
    cat(r$status, "\n", paste0("  ", r$details, collapse = "\n"), "\n", sep = "")
    if (r$status == "FAIL") quit(status = 1)
}
