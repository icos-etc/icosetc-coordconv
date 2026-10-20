# features.R -- spatial polygon features of a station: Target Area (TA),
# Exclusion Areas (EA) and Land Cover Typologies (LCT).
#
# Inputs : files unzipped from the station's KML and shapefile downloads.
# Outputs: list(TA = sf, EA = sf (union of all EA), LCT = list of sf), WGS84.
# KML files are preferred over shapefiles for the same layer kind. Invalid
# geometries are repaired with st_make_valid(); EA polygons are unioned in a
# local AEQD plane centred on the tower.

#' Layer kind of a feature file from its name: "TA", "EA", "LCT" or NA.
feature_kind <- function(path) {
    b <- toupper(tools::file_path_sans_ext(basename(path)))
    ifelse(grepl("(^|[^A-Z])LCT", b), "LCT",
    ifelse(grepl("(^|[^A-Z])EA([^A-Z]|$)", b), "EA",
    ifelse(grepl("(^|[^A-Z])TA([^A-Z]|$)", b), "TA", NA_character_)))
}

# Read one KML/SHP file as polygons in WGS84 (geometry + Name).
.read_feature_file <- function(path, fallback_crs = NULL) {
    lyr <- sf::st_layers(path)$name[1]
    x <- sf::st_read(path, layer = lyr, quiet = TRUE)
    x <- sf::st_zm(x, drop = TRUE)
    if (is.na(sf::st_crs(x))) {
        if (is.null(fallback_crs)) stop("no CRS in ", basename(path), " and no fallback CRS")
        log_warn("features: ", basename(path), " has no CRS; using the Target Area CRS")
        sf::st_crs(x) <- fallback_crs
    }
    x <- sf::st_transform(x, 4326)
    nm <- if ("Name" %in% names(x)) as.character(x$Name) else tools::file_path_sans_ext(basename(path))
    sf::st_sf(ID = nm, geometry = sf::st_geometry(x))
}

#' Read the TA, EA and LCT layers from a set of files.
#'
#' @param files paths of unzipped files (.kml, .shp and side-car files).
#' @return list(TA, EA (list), LCT (list)) of sf objects; elements may be empty.
read_feature_layers <- function(files) {
    files <- files[grepl("\\.(kml|shp)$", files, ignore.case = TRUE)]
    kind <- feature_kind(files)
    pick <- function(k) {
        f <- files[kind %in% k]
        if (any(grepl("\\.kml$", f, ignore.case = TRUE))) f[grepl("\\.kml$", f, ignore.case = TRUE)] else f
    }
    for (f in pick("TA")[grepl("\\.shp$", pick("TA"), ignore.case = TRUE)]) {
        sidecars <- paste0(tools::file_path_sans_ext(f), c(".shx", ".dbf"))
        if (!all(file.exists(sidecars))) stop("shapefile ", basename(f), " lacks .shx/.dbf side-car files")
    }
    ta_files <- pick("TA")
    TA <- if (length(ta_files)) do.call(rbind, lapply(ta_files, .read_feature_file)) else NULL
    ta_crs <- if (!is.null(TA)) sf::st_crs(TA) else NULL
    EA  <- lapply(pick("EA"), .read_feature_file, fallback_crs = ta_crs)
    LCT <- lapply(sort(pick("LCT")), .read_feature_file, fallback_crs = ta_crs)
    ign <- sum(is.na(kind))
    if (ign) log_info("features: ", ign, " file(s) not recognised as TA/EA/LCT ignored")
    list(TA = TA, EA = EA, LCT = LCT)
}

# Repair invalid geometries (logged).
.make_valid <- function(x, what) {
    if (is.null(x) || !nrow(x)) return(x)
    bad <- !sf::st_is_valid(x)
    if (any(bad, na.rm = TRUE)) {
        log_warn("features: ", sum(bad, na.rm = TRUE), " invalid ", what, " geometr(y/ies) repaired")
        x <- sf::st_make_valid(x)
    }
    x
}

#' Clean the layers and union the exclusion areas.
#'
#' @param layers output of read_feature_layers().
#' @param tower_lat,tower_lon AEQD centre for planar operations.
#' @return list(TA = sf|NULL, EA = sf|NULL, LCT = list of sf).
prepare_features <- function(layers, tower_lat, tower_lon) {
    old <- sf::sf_use_s2(FALSE); on.exit(suppressMessages(sf::sf_use_s2(old)))
    crs_local <- sf::st_crs(aeqd_crs(tower_lat, tower_lon))
    TA <- .make_valid(layers$TA, "TA")
    EA <- NULL
    if (length(layers$EA)) {
        ea <- do.call(rbind, layers$EA)
        ea <- .make_valid(sf::st_transform(ea, crs_local), "EA")
        u <- tryCatch(sf::st_union(sf::st_buffer(ea, 0)), error = function(e) {
            log_warn("features: EA union failed (", conditionMessage(e), "); EA polygons kept separate")
            NULL
        })
        EA <- if (is.null(u)) sf::st_transform(ea, 4326) else sf::st_sf(ID = "EA", geometry = sf::st_transform(u, 4326))
    }
    LCT <- lapply(layers$LCT, .make_valid, what = "LCT")
    list(TA = TA, EA = EA, LCT = LCT)
}
