# Data access and placing the datasets

**Raw datasets not included.** `data/raw/` is git-ignored. To reproduce this analysis, download the files from the websites listed below and organize them according to the folder structure shown.

## 1. 2023 NNS Public Use Files (required) – © DOST-FNRI
Download from the DOST-FNRI National Nutrition Survey Public Use Files portal (eNutrition site,
https://enutrition.fnri.dost.gov.ph/ ; registration/terms may apply). **Read and follow DOST-FNRI's terms of use
and request permission before redistributing any file.** Required folders:
```
data/raw/nns_puf/2023_Biochem/2023_Biochem_Dataset.csv
data/raw/nns_puf/2023_Biochem/2023_Biochem_Data Dictionary.csv
data/raw/nns_puf/2023_Clinical and Health/2023_Clinical_Dataset.csv
```

## 2. Province/city boundaries (required only with `--recompute-elevation`)
BetterGov GeoJSON: unzip `geojson_province.zip` **and** `geojson_city.zip` (and municipality files for Davao
Occidental) anywhere under `data/raw/bettergov/`. Files are named `ph.<region>.<province>.<unit>.any.geo.json`.
Check the source's licence and cite it.

## 3. Elevation grid (required only with `--recompute-elevation`)
`data/raw/figshare/PHILIPPINES_altitude_data.csv` (columns x, y, alt; ~0.0083° cells, ~84 MB, sea cells = NA).
It is identical to the HDX raster `PHL_msk_alt`, i.e. the same source. **Record its original source, DOI and
licence yourself** (not documented in the file). Too large to commit comfortably; do not add it to Git.

## 4. JAXA 30 m DSM (optional validation)
Unzip tiles `ALPSMLC30_N0xxE1xx_DSM.tif` under `data/raw/jaxa/`; run `src/validate_elevation.py`. Follow JAXA's
registration/terms and cite JAXA.

## What you DO commit
Code (`src/`, `tests/`), `requirements.txt`, docs, `data/derived/province_elevation.csv` and aggregated
`outputs/`. Never commit individual-level records (including a merged analytic file).
