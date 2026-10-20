# settings_example.R -- template for the local settings file.
#
# Copy this file to config/settings_local.R (git-ignored) and fill in the real
# values, or point to another file with --settings <file> or the environment
# variable COORDCONV_SETTINGS. Never commit real values.
#
# Endpoint templates are appended to api_base_url with a "/" and may use the
# placeholders {site} (lower-case station code), {var} (BADM variable group)
# and {type} (kplot download type: plot, kml or shape).

settings <- list(
    api_base_url      = "https://<host>:<port>",
    endpoints         = list(
        badm = "<BADM records endpoint template, uses {site} and {var}>",
        plot = "<plot/spatial-data endpoint template, uses {site} and {type}>"
    ),
    # BADM variable names of the two network-wide station lists
    station_list_vars = c(class = "<station class variable>", igbp = "<IGBP variable>"),
    excluded_sites    = c("XX-Aaa", "XX-Bbb"),   # test or non-network stations
    out_dir           = "output",
    tls_verify        = TRUE,                    # FALSE requires the 'curl' package

    # QC flag thresholds (records are flagged for review, never withheld)
    qc = list(
        format_disagreement_max_m = 1.0,   # N/E vs polar positions of the same record
        large_shift_m             = 20     # new vs previous calculated position
    ),

    # Future DB upload (not implemented in v1): placeholders only
    db_driver   = "<odbc driver>",
    db_server   = "<host,port>",
    db_name     = "<database>",
    db_user     = "<user>",
    db_password = "<read from an environment variable, never stored here>"
)
