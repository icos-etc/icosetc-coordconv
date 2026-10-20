# icosetc-coordconv

Coordinate conversion for ICOS ecosystem stations. The positions of sampling plots and of
sensors and ancillary measurements are submitted as **lat/long**, **north/east offsets** or
**polar coordinates** (distance and azimuth) from a reference point. These scripts compute
lat/long for every record, flag doubtful records for review, and write coordinate tables
(TSV) and KML layers.

## Conventions

- **Units:** metres and decimal degrees, WGS84.
- **North/east offsets are true-north offsets.** Each record is computed in a local
  azimuthal equidistant projection centred on its reference point, so distances and
  azimuths are geodesic. UTM grid north is not used.
- **Azimuths** are in degrees, clockwise from true north, `[0, 360)`.
  Polar to offsets: `E = d sin(a)`, `N = d cos(a)`.
- **Magnetic north.** When a record's north reference is "Magnetic", the declination `D`
  (east-positive, IGRF-14 via the R package `oce`) is evaluated at the record's reference
  point and date:
  - azimuths: `true = magnetic + D`;
  - offsets: `E_t = E cos D + N sin D`, `N_t = N cos D - E sin D`.
- **References:**
  - the tower is the `GRP_LOCATION` position valid at the record date;
  - CP and SP-I plots are referenced to the tower;
  - SP-II plots are referenced to their SP-I;
  - BADM records use the plot named in `<X>_PLOT`, or the tower if it is empty or "outside".
- **Format precedence:** polar > N/E > lat/long.
- **GAI:** the instrument pair (hemispherical photo or ceptometer) is chosen from `GAI_METHOD`.

## Requirements

- R >= 4.4 with `sf` (GDAL, GEOS and PROJ >= 6.3), `data.table` and `oce` (>= 1.8-4).
- `oce` is only needed for magnetic records. Without it everything else runs, and
  magnetic records stop with a clear error.
- `geosphere` is optional; the tests use it for a cross-check.
- Packages are never installed automatically.

## Setup

```sh
git clone https://github.com/icos-etc/icosetc-coordconv.git
cd icosetc-coordconv
cp config/settings_example.R config/settings_local.R   # edit; this file is git-ignored
```

`config/settings_local.R` holds the API address and endpoints, the excluded stations and the
QC thresholds. It is never committed. Offline runs do not need it.

## Usage

```sh
Rscript run_coordconv.R --sites <codes|all> [options]
```

| Option | Meaning |
|---|---|
| `--sites` | comma-separated station codes, or `all` (excluded stations are skipped) |
| `--what` | any of `plots`, `badm`, `features` (default `plots,badm`) |
| `--kml` | also write KML layers of plots and BADM points |
| `--out` | output root (default `output/`) |
| `--cache` | keep the raw downloads under `<out>/<run_id>/raw/` (a reproducible snapshot) |
| `--offline` | read a raw directory instead of the API |
| `--previous` | merged files of a previous run, used for `LARGE_SHIFT` where the API has no previous value |
| `--settings` | settings file (default: `COORDCONV_SETTINGS`, then `config/settings_local.R`) |

**Offline example.** This works on a fresh clone, with no settings file and no network:

```sh
Rscript run_coordconv.R --sites XX-Tst --what plots,badm --kml --offline examples/raw --out example_output
Rscript examples/run_example.R   # same run, compared with examples/expected
Rscript tests/run_tests.R        # unit tests and the example end to end
```

The example station `XX-Tst` and all its data are synthetic.

**Live mode** downloads from the ICOS Ecosystem Thematic Centre internal API, so it needs
access to that API. Anyone else can use offline mode with their own raw files in this layout
(tab-separated, header row, one file per table):

```
<raw_dir>/_network/station_class.tsv   site_code, CLASS_ICOS [, DATE_LABEL]
<raw_dir>/_network/station_igbp.tsv    site_code, IGBP
<raw_dir>/<SITE>/plot.tsv              ID, FEAT (plot ID: TOWER, CP_nn, SP-I_nn, SP-II_nn-mm), LAT, LON,
                                       PLOT_EASTWARD_DIST, PLOT_NORTHWARD_DIST, PLOT_DISTANCE_POLAR,
                                       PLOT_ANGLE_POLAR, PLOT_RECT_X/Y/DIR, PLOT_NORTHREF, PLOT_DATE
<raw_dir>/<SITE>/GRP_LOCATION.tsv      ID, LOCATION_LAT, LOCATION_LONG, LOCATION_DATE
<raw_dir>/<SITE>/GRP_<X>.tsv           ID plus the BADM columns of the group (see R/group_specs.R)
<raw_dir>/<SITE>/features_kml.zip      optional: Target Area, Exclusion Area and Land Cover
                                       Typology polygons (file names containing _TA, _EA, _LCT)
```

## Outputs

Each run writes to `<out>/<run_id>/`:

- `sites/<SITE>/<SITE>-GRP_<X>_update.tsv`: record `ID` plus the calculated columns, only for
  records submitted as N/E or polar.
- `sites/<SITE>/<SITE>-GRP_<X>_full.tsv`: every record with its input format, reference,
  declination, N/E, polar, lat/long and status.
- `merged/GRP_<X>.tsv`: the update files of all stations, one per group.
- `merged/QC_flags.tsv` (and `<SITE>-QC_flags.tsv`): records flagged for review. Flagged
  records are still included in the update files.
  - `FORMAT_DISAGREEMENT`: N/E and polar of one instrument pair differ.
  - `OUTSIDE_TA`: point outside the target area.
  - `LARGE_SHIFT`: point moved far from its previous value.
  - `NOT_COMPUTED`: position could not be computed.
  - `GAI_METHOD_UNKNOWN`, `GAI_PAIR_MISMATCH`: GAI instrument pair problems.
- `sites/<SITE>/kml/*.kml`: KML layers of the tower, plots, plot areas, BADM points (with
  `--kml`) and polygon features (with `--what features`).
- `coordconv.log`, `summary.tsv`, `stations.tsv`: log, counts per group, status per station.

Rounding: lat/long 8 decimals, distances 2, angles 3. Files are tab-separated with LF line
endings.

## License

BSD 3-Clause, see [LICENSE](LICENSE).

Provided as is. No support is provided by the ICOS Ecosystem Thematic Centre.
