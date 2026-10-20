# write_outputs.R -- writers: rounded TSV tables, merged group files, KML layers, run summary.
#
# Inputs : processed tables, sf layers.
# Outputs: files under output/<run>/ (see README):
#   sites/<SITE>/<SITE>-GRP_<X>_update.tsv, _full.tsv, kml/<SITE>_<layer>.kml
#   merged/GRP_<X>.tsv, summary.tsv
# Numbers are rounded at write time: lat/long 8 dp, angles (incl. declination)
# 3 dp, distances 2 dp.

#' Number of decimals for a column name (NA = leave unchanged).
.digits_for <- function(name) {
    if (grepl("(angle|azimuth|magdec)", name, ignore.case = TRUE)) return(3)
    if (grepl("(lat|long?)$", name, ignore.case = TRUE)) return(8)
    if (grepl("(dist|distance|east|north|_m$)", name, ignore.case = TRUE)) return(2)
    NA
}

#' Round the numeric (double) columns of a table according to their names.
round_for_output <- function(df) {
    for (nm in names(df)) {
        if (is.double(df[[nm]])) {
            d <- .digits_for(nm)
            if (!is.na(d)) df[[nm]] <- round(df[[nm]], d)
        }
    }
    df
}

#' Write a tab-separated table (rounded, NA as empty, no scientific notation).
write_tsv <- function(df, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(round_for_output(df), path, sep = "\t", na = "", quote = FALSE, scipen = 100,
                       eol = "\n")   # same bytes on Windows and Linux
    invisible(path)
}

#' Write the update and full tables of one station.
#'
#' @param site station code.
#' @param results named list (group -> list(update, full)).
#' @param dir station output directory.
write_site_tables <- function(site, results, dir) {
    for (g in names(results)) {
        r <- results[[g]]
        if (!is.null(r$full) && nrow(r$full)) write_tsv(r$full, file.path(dir, sprintf("%s-%s_full.tsv", site, g)))
        if (!is.null(r$update) && nrow(r$update)) write_tsv(r$update, file.path(dir, sprintf("%s-%s_update.tsv", site, g)))
    }
}

#' Merge the update tables of all stations into one file per group.
#'
#' @param updates list (one per station) of named lists group -> update table.
#' @param dir output directory (merged/).
merge_group_files <- function(updates, dir) {
    groups <- unique(unlist(lapply(updates, names)))
    for (g in groups) {
        tabs <- Filter(function(x) !is.null(x) && nrow(x), lapply(updates, `[[`, g))
        if (length(tabs)) write_tsv(data.table::rbindlist(tabs, use.names = TRUE, fill = TRUE),
                                    file.path(dir, paste0(g, ".tsv")))
    }
}

#' KML colour (aabbggrr) from "#RRGGBB" and alpha (0-255).
.kml_colour <- function(hex, alpha = 255) {
    h <- sub("^#", "", hex)
    sprintf("%02x%s%s%s", alpha, substr(h, 5, 6), substr(h, 3, 4), substr(h, 1, 2))
}

#' Write one KML layer with the GDAL KML driver and a single per-layer style.
#'
#' The plain KML driver ignores feature styles, so one <Style> is inserted in
#' the document and referenced by every placemark (<styleUrl>).
#'
#' @param x sf object (any geometry), in or transformable to WGS84.
#' @param path output .kml file.
#' @param name_col column used as placemark name.
#' @param colour "#RRGGBB".
#' @param layer layer (folder) name.
write_kml_layer <- function(x, path, name_col = "ID", colour = "#FF0000", layer = NULL) {
    if (is.null(x) || !nrow(x)) return(invisible(NULL))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    x <- sf::st_transform(sf::st_zm(x, drop = TRUE), 4326)
    out <- sf::st_sf(Name = as.character(x[[name_col]]),
                     Description = if ("type" %in% names(x)) as.character(x$type) else "",
                     geometry = sf::st_geometry(x))
    if (file.exists(path)) unlink(path)
    layer <- gsub("[^A-Za-z0-9_]", "_", layer %||% tools::file_path_sans_ext(basename(path)))  # valid XML id
    sf::st_write(out, path, driver = "KML", layer = layer, quiet = TRUE)
    txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    style <- sprintf(paste0('<Style id="layer_style"><IconStyle><color>%s</color><scale>0.8</scale></IconStyle>',
                            '<LabelStyle><scale>0.7</scale></LabelStyle>',
                            '<LineStyle><color>%s</color><width>2</width></LineStyle>',
                            '<PolyStyle><color>%s</color></PolyStyle></Style>'),
                     .kml_colour(colour), .kml_colour(colour), .kml_colour(colour, 80))
    txt <- sub("(<Document[^>]*>)", paste0("\\1\n", style), txt, perl = TRUE)
    txt <- gsub("(?s)(<Placemark[^>]*>\\s*(<name>.*?</name>\\s*)?(<description>.*?</description>\\s*)?)",
                "\\1<styleUrl>#layer_style</styleUrl>\n", txt, perl = TRUE)
    writeLines(txt, path, useBytes = TRUE)
    invisible(path)
}

#' Points data.frame(ID, lat, lon[, type]) -> sf.
points_sf <- function(df, id_col = "ID", lat_col = "lat", lon_col = "lon") {
    df <- df[!is.na(df[[lat_col]]) & !is.na(df[[lon_col]]), , drop = FALSE]
    if (!nrow(df)) return(NULL)
    sf::st_as_sf(data.frame(ID = df[[id_col]], type = if ("type" %in% names(df)) df$type else "",
                            lon = df[[lon_col]], lat = df[[lat_col]]),
                 coords = c("lon", "lat"), crs = 4326)
}

# Per-layer colours (simple, fixed)
KML_COLOURS <- c(tower = "#000000", TA = "#2CA02C", EA = "#D62728", LCT = "#BCBD22",
                 CP = "#1F77B4", SPI = "#FF7F0E", SPII = "#9467BD", BADM = "#17BECF")

#' Write the KML layers of one station.
#'
#' @param site station code.
#' @param dir station output directory (kml/ is created inside).
#' @param plots result of process_plots() (or NULL).
#' @param badm result of process_badm() (or NULL; written when badm_kml = TRUE).
#' @param features result of prepare_features() (or NULL).
write_site_kml <- function(site, dir, plots = NULL, badm = NULL, features = NULL, plot_kml = TRUE, badm_kml = TRUE) {
    kdir <- file.path(dir, "kml")
    f <- function(layer) file.path(kdir, sprintf("%s_%s.kml", site, layer))
    if (plot_kml && !is.null(plots)) {
        pts <- plots$points
        write_kml_layer(points_sf(transform(pts[pts$type == "TOWER", ], ID = "EC tower")), f("EC-tower"), colour = KML_COLOURS[["tower"]])
        for (tp in c("CP", "SP-I", "SP-II")) {
            sel <- pts[pts$type == tp, ]
            if (nrow(sel)) write_kml_layer(points_sf(transform(sel, ID = plot_id)), f(gsub("-", "", tp)),
                                           colour = KML_COLOURS[[gsub("-", "", tp)]])
        }
        a <- plots$areas
        if (!is.null(a)) {
            if (any(a$type == "CP area"))   write_kml_layer(a[a$type == "CP area", ], f("CP_areas"), colour = KML_COLOURS[["CP"]])
            if (any(a$type == "SP-I area")) write_kml_layer(a[a$type == "SP-I area", ], f("SPI_areas"), colour = KML_COLOURS[["SPI"]])
        }
    }
    if (badm_kml && !is.null(badm)) {
        for (g in names(badm)) {
            fu <- badm[[g]]$full
            if (is.null(fu) || !nrow(fu)) next
            fu <- fu[fu$status == "ok" & fu$input_format %in% c("ne", "polar"), , drop = FALSE]
            if (nrow(fu)) write_kml_layer(points_sf(data.frame(ID = fu$ID, lat = fu$calc_lat, lon = fu$calc_lon, type = g)),
                                          f(g), colour = KML_COLOURS[["BADM"]])
        }
    }
    if (!is.null(features)) {
        if (!is.null(features$TA)) write_kml_layer(features$TA, f("TA"), colour = KML_COLOURS[["TA"]])
        if (!is.null(features$EA)) write_kml_layer(features$EA, f("EA"), colour = KML_COLOURS[["EA"]])
        for (i in seq_along(features$LCT))
            write_kml_layer(features$LCT[[i]], f(sprintf("LCT%02d", i)), colour = KML_COLOURS[["LCT"]])
    }
    invisible(kdir)
}

#' Write the run summary table.
write_summary <- function(rows, path) {
    write_tsv(rows, path)
}
