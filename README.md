# Anemia Prevalence by Age Group, Sex, and Severity in the Philippines: A Cross-Sectional Analysis of the 2023 National Nutrition Survey Public Use Files

Reproducible code for Tables 1-3 and Supplementary Table S1: weighted anemia prevalence by population group
(under-5 children, school-age children, adolescents, adults, older persons, non-pregnant non-lactating women,
pregnant women by trimester, lactating mothers), classified with the WHO 2024 haemoglobin guideline, with
elevation and smoking adjustment. 

## Acknowledgement
I gratefully acknowledge the **Department of Science and Technology – Food and Nutrition Research Institute
(DOST-FNRI)** for conducting the 2023 National Nutrition Survey and making the Public Use Files available.
The analysis, interpretation and any errors are solely those of the authors; DOST-FNRI is not responsible for them.

> Suggested data citation (complete the bracketed parts): DOST-FNRI. *2023 National Nutrition Survey, Public Use
> Files*. Taguig City: Department of Science and Technology – Food and Nutrition Research Institute;
> October 2026.

## Repository layout
```
src/nns_anemia_who2024.py   full pipeline (Steps 0-7, commented)
src/validate_elevation.py   optional check of the elevation against a 30 m DSM
tests/test_who2024.py       unit tests of the WHO rules
data/raw/                   YOUR copies of the source data (git-ignored; see DATA_ACCESS.md)
data/derived/               province_elevation.csv (derived, safe to commit)
outputs/                    tables, participant flow, EDA summaries
```

## Quick start
```bash
python -m venv .venv && source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
# place data as described in DATA_ACCESS.md, then:
python src/nns_anemia_who2024.py                          # Tables 1-3, S1 (cached elevation)
python src/nns_anemia_who2024.py --recompute-elevation --sensitivity --subnational
# tested with Python 3.12
python -m pytest -q                                      # optional unit tests
```

## What each step does (all in `src/nns_anemia_who2024.py`)
| Step | Purpose |
|---|---|
| 0 | Libraries, paths, constants (region labels, weight variable) |
| 1 | Load Biochem (Hb, age, sex, `psc`, weights, design) and Clinical (`smoke_status`), merge on `hhnum`+`member_code` |
| 2 | Eligibility: valid Hb, positive weight, plausible Hb (1-30 g/dL), age ≥6 months; writes `participant_flow.csv` (STROBE item 13) |
| 2b | EDA and data-quality flags (psc anomalies, smoking missingness) in `outputs/eda/` |
| 3 | Elevation: match 117 NNS provinces/HUCs to boundary polygons, summarise a ~1 km grid (mean, median, percentiles, % ≥500 m) |
| 4 | WHO 2024 rules: cutoffs by age/sex/trimester, elevation equation (≥500 m), smoking 3 g/L, severity (Table 3), WHO 2011 comparator |
| 5 | Population groups (adolescents/adults exclude pregnant, lactating and females with `psc` ≠ 0) |
| 6 | Survey estimation: weights, stratified PSU-clustered Taylor linearisation, logit 95% CIs |
| 7 | Build Tables 1-3, S1; optional sensitivity (elevation scenarios) and region/province tables |

## Elevation: read this before citing results
The PUF contains **no coordinates**, so elevation is assigned at **province/HUC level** (median of the 1 km grid
cells inside each polygon). This is an ecological approximation that mainly affects highland provinces
(e.g. Baguio, Benguet, Mountain Province, Ifugao, Kalinga, Abra, Nueva Vizcaya, Bukidnon). Elevation is weighted by
land area, not population. Use `--sensitivity` to see the effect. Best practice: obtain PSU/barangay (PSGC)
identifiers from DOST-FNRI and extract elevation at sampling-unit level.

## Key assumptions (edit in the script if documentation says otherwise)
- Strata = province/HUC; PSU = province × `ms_psucode`; weight = `fwgti_hemoglobin` (confirm with FNRI documentation).
- Current smoker = `smoke_status == 1`, age ≥10; cigarettes/day unknown, so 3 g/L; missing smoking = non-smoker.
- WHO cutoffs end at age 65: adults >65 use the adult cutoffs; lactating mothers use non-pregnant women cutoffs.
- WHO 2011 comparator uses no altitude/smoking adjustment.

## Data availability statement (template)
The 2023 NNS Public Use Files are available from DOST-FNRI [https://enutrition.fnri.dost.gov.ph/]. Boundary and elevation data are
from https://figshare.com/articles/dataset/Philippines_elevation_data/11777733?file=21494226, https://data.bettergov.ph/datasets/7,https://data.humdata.org/dataset/philippines-elevation-model, and https://portal.opentopography.org/raster?opentopoID=OTALOS.112016.4326.2. Code, derived province elevation and aggregate results are in this repository
[November 2026].
