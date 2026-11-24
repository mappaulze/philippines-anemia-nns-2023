#!/usr/bin/env python3
"""

Anemia Prevalence by Age Group, Sex, and Severity in the Philippines: 
A Cross-Sectional Analysis of the 2023 National Nutrition Survey 
Public Use Files

Reproduces: Table 1 (prevalence), Table 2 (sex), Table 3 (severity),
Supplementary Table S1 (WHO 2011 vs WHO 2024), plus optional sensitivity and
regional/provincial tables.

Run from the repository root:
    python src/nns_anemia_who2024.py                       # tables (uses cached elevation)
    python src/nns_anemia_who2024.py --recompute-elevation
    python src/nns_anemia_who2024.py --sensitivity --subnational

Data acknowledgement: the 2023 National Nutrition Survey Public Use Files were produced
by the Department of Science and Technology - Food and Nutrition Research Institute
(DOST-FNRI). See README.md.
"""
# ----------------------------------------------------------------------------------
# STEP 0. LIBRARIES AND SETTINGS
# numpy/pandas: data handling | scipy: logit confidence intervals |
# matplotlib.path: point-in-polygon (used instead of geopandas so only common
# packages are required).
# ----------------------------------------------------------------------------------
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from scipy.special import expit, logit
from scipy.stats import t as t_dist

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
OUT = ROOT / "outputs"

NNS_DIR = RAW / "nns_puf"                       # 2023 NNS PUF folders (see DATA_ACCESS.md)
BIOCHEM = NNS_DIR / "2023_Biochem" / "2023_Biochem_Dataset.csv"
BIOCHEM_DICT = NNS_DIR / "2023_Biochem" / "2023_Biochem_Data Dictionary.csv"
CLINICAL = NNS_DIR / "2023_Clinical and Health" / "2023_Clinical_Dataset.csv"
BOUNDARIES = RAW / "bettergov"                  # BetterGov GeoJSON (province + city files)
ALT_GRID = RAW / "figshare" / "PHILIPPINES_altitude_data.csv"   # ~1 km elevation grid (x, y, alt)
ELEV_CSV = DERIVED / "province_elevation.csv"   # cached province/HUC elevation statistics

# Survey design (ASSUMPTION - confirm with DOST-FNRI documentation):
#   strata = province/HUC (provhuc); primary sampling unit = provhuc x ms_psucode.
WEIGHT = "fwgti_hemoglobin"

REGIONS = {13: "National Capital Region (NCR)", 14: "Cordillera Administrative Region (CAR)",
           1: "Region I (Ilocos Region)", 2: "Region II (Cagayan Valley)", 3: "Region III (Central Luzon)",
           41: "Region IV-A (CALABARZON)", 42: "Region IV-B (MIMAROPA)", 5: "Region V (Bicol Region)",
           6: "Region VI (Western Visayas)", 7: "Region VII (Central Visayas)", 8: "Region VIII (Eastern Visayas)",
           9: "Region IX (Zamboanga Peninsula)", 10: "Region X (Northern Mindanao)", 11: "Region XI (Davao Region)",
           12: "Region XII (SOCCSKSARGEN)", 16: "Region XIII (CARAGA)", 15: "BARMM"}


# ----------------------------------------------------------------------------------
# STEP 1. LOAD AND MERGE
# Biochem file = haemoglobin, age, sex, physiological status code (psc), weights, design.
# Clinical file = smoking status. Both keyed by household number + member code.
# psc: 0 = not pregnant/lactating (females), 1-9 = months pregnant, 10-12 = lactating,
# 99 = not applicable (males).
# ----------------------------------------------------------------------------------
def load_data():
    bio = pd.read_csv(BIOCHEM, low_memory=False)
    cli = pd.read_csv(CLINICAL, low_memory=False)[["hhnum", "member_code", "smoke_status"]]
    assert not bio.duplicated(["hhnum", "member_code"]).any(), "duplicate person keys in Biochem"
    assert not cli.duplicated(["hhnum", "member_code"]).any(), "duplicate person keys in Clinical"
    return bio.merge(cli, on=["hhnum", "member_code"], how="left")


# ----------------------------------------------------------------------------------
# STEP 2. ELIGIBILITY CHECK (participant flow, STROBE item 13)
# Keep persons with a measured Hb, a positive weight and a plausible Hb
# (1-30 g/dL, as in Hernandez-Vasquez et al. 2024) aged >= 6 months.
# ----------------------------------------------------------------------------------
def eligibility(df):
    flow = [("Records in Biochem file", len(df))]
    d = df[df["hemoglobin"].notna()]; flow.append(("With a haemoglobin value", len(d)))
    d = d[d[WEIGHT] > 0]; flow.append(("With a positive Hb weight", len(d)))
    d = d[d["hemoglobin"].between(1, 30)]; flow.append(("With plausible Hb (1-30 g/dL)", len(d)))
    d = d[d["age"] >= 0.5].reset_index(drop=True); flow.append(("Aged >=6 months (analytic sample)", len(d)))
    flow = pd.DataFrame(flow, columns=["step", "n"])
    flow["excluded"] = -flow["n"].diff().fillna(0).astype(int)
    # current smoking: smoke_status == 1; WHO smoking adjustment only applied from age 10
    d["smoker"] = ((d["smoke_status"] == 1) & (d["age"] >= 10)).astype(int)
    d["female"] = d["sex"] == 2
    return d, flow


# ----------------------------------------------------------------------------------
# STEP 2b. EXPLORATORY DATA ANALYSIS
# Ranges, psc/age/sex cross-checks and data-quality flags. Inspect these CSVs before
# trusting any table (e.g. psc anomalies such as pregnancy at <15 or >=50 years).
# ----------------------------------------------------------------------------------
def eda(d):
    e = OUT / "eda"; e.mkdir(parents=True, exist_ok=True)
    d[["age", "hemoglobin", WEIGHT]].describe().T.to_csv(e / "numeric_summary.csv")
    pd.crosstab(pd.cut(d.age, [0.5, 5, 10, 20, 50, 60, 120], right=False), [d.sex, d.psc]).to_csv(e / "age_by_sex_psc.csv")
    d.groupby("provhuc").agg(n=("age", "size"), n_psu=("ms_psucode", "nunique")).to_csv(e / "sample_by_province.csv")
    flags = {"pregnant <15 y": ((d.psc.between(1, 9)) & (d.age < 15)).sum(),
             "pregnant >=50 y": ((d.psc.between(1, 9)) & (d.age >= 50)).sum(),
             "lactating >=50 y": ((d.psc.between(10, 12)) & (d.age >= 50)).sum(),
             "smoking status missing, age>=10": ((d.age >= 10) & d.smoke_status.isna()).sum(),
             "current smokers (age>=10)": d.smoker.sum()}
    pd.DataFrame({"check": list(flags), "n": [int(v) for v in flags.values()]}).to_csv(e / "data_quality_flags.csv", index=False)
    print("EDA written to", e)


# ----------------------------------------------------------------------------------
# STEP 3. ELEVATION (province/HUC level) - see README "Elevation" for limitations.
# The PUF has NO coordinates, so each NNS province/HUC is matched to a BetterGov polygon
# and summarised from a ~1 km elevation grid (cells whose centre falls inside it).
# ----------------------------------------------------------------------------------
HUC = {  # NNS label -> (province slug, unit slug) of the BetterGov city polygon
    "Angeles City": ("pampanga", "angeles-city"), "Olongapo City": ("zambales", "olongapo-city"),
    "Lucena City": ("quezon", "lucena-city"), "Iloilo City": ("iloilo", "iloilo-city"),
    "Bacolod City": ("negros-occidental", "bacolod-city"), "Cebu City": ("cebu", "cebu-city"),
    "Lapu-Lapu City (Opon)": ("cebu", "lapu-lapu-city"), "Mandaue City": ("cebu", "mandaue-city"),
    "Tacloban City": ("leyte", "tacloban-city"), "Zamboanga City": ("zamboanga-del-sur", "zamboanga-city"),
    "City of Isabela": ("basilan", "isabela-city"), "Iligan City": ("lanao-del-norte", "iligan-city"),
    "Cagayan De Oro City": ("misamis-oriental", "cagayan-de-oro-city"), "Davao City": ("davao-del-sur", "davao-city"),
    "General Santos City": ("south-cotabato", "general-santos-city"), "Cotabato City": ("maguindanao", "cotabato-city"),
    "Baguio City": ("benguet", "baguio-city"), "Butuan City": ("agusan-del-norte", "butuan-city"),
    "Puerto Princesa City": ("palawan", "puerto-princesa-city"), "City of Manila": ("metropolitan-manila", "manila"),
    "Mandaluyong City": ("metropolitan-manila", "mandaluyong"), "City of Marikina": ("metropolitan-manila", "marikina"),
    "City of Pasig": ("metropolitan-manila", "pasig-city"), "Quezon City": ("metropolitan-manila", "quezon-city"),
    "San Juan City": ("metropolitan-manila", "san-juan"), "Caloocan City": ("metropolitan-manila", "kalookan-city"),
    "City of Malabon": ("metropolitan-manila", "malabon"), "City of Navotas": ("metropolitan-manila", "navotas"),
    "City of Valenzuela": ("metropolitan-manila", "valenzuela"), "Las Piñas City": ("metropolitan-manila", "las-pias"),
    "Makati City": ("metropolitan-manila", "makati-city"), "City of Muntinlupa": ("metropolitan-manila", "muntinlupa"),
    "City of Parañaque": ("metropolitan-manila", "paraaque"), "Pasay City": ("metropolitan-manila", "pasay-city"),
    "Pateros": ("metropolitan-manila", "pateros"), "Taguig City": ("metropolitan-manila", "taguig")}
PROV_SLUG = {"Cagayan Province": "cagayan", "Isabela Province": "isabela", "Quezon Province": "quezon",
             "Western Samar": "samar", "North Cotabato": "north-cotabato"}
DAVAO_OCC = [("davao-del-sur", u) for u in ["don-marcelino", "jose-abad-santos", "malita", "santa-maria", "sarangani"]]


def index_geojson():
    """Index BetterGov files by (province slug, unit slug); names look like
    ph.<region>.<province>.<unit>.any.geo.json ('any' unit = whole province)."""
    idx = {}
    for f in BOUNDARIES.rglob("*.geo.json"):
        p = f.name.split(".")
        if len(p) == 7 and p[4] == "any":
            idx[(p[2], p[3])] = f
    return idx


def load_polygons(path):
    g = json.load(open(path)); geom = g["geometry"]
    polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
    return [[np.array(r) for r in poly] for poly in polys]          # polygon -> [outer ring, holes...]


def points_inside(pts, polys):
    """Ray-casting point-in-polygon (outer rings minus holes), bounding-box prefiltered."""
    mask = np.zeros(len(pts), bool)
    for rings in polys:
        o = rings[0]; x0, y0 = o.min(0); x1, y1 = o.max(0)
        sel = np.flatnonzero((pts[:, 0] >= x0) & (pts[:, 0] <= x1) & (pts[:, 1] >= y0) & (pts[:, 1] <= y1))
        if sel.size == 0:
            continue
        m = MplPath(o).contains_points(pts[sel])
        for h in rings[1:]:
            m &= ~MplPath(h).contains_points(pts[sel])
        mask[sel[m]] = True
    return mask


def province_labels():
    d = pd.read_csv(BIOCHEM_DICT, encoding="latin1")
    p = d[d.variable_name == "provhuc"][["var_value", "value_label"]].copy()
    p.columns = ["code", "label"]; p["code"] = p["code"].astype(int); p["label"] = p["label"].str.strip()
    return p.reset_index(drop=True)


def build_province_elevation():
    idx = index_geojson(); lab = province_labels()
    grid = pd.read_csv(ALT_GRID, na_values="NA").dropna(subset=["alt"]).reset_index(drop=True)
    pts, alt = grid[["x", "y"]].values, grid["alt"].values
    rows = []
    for r in lab.itertuples():
        if r.label in HUC:
            polys, excl = load_polygons(idx[HUC[r.label]]), []
        elif r.label == "Davao Occidental":
            polys = sum([load_polygons(idx[k]) for k in DAVAO_OCC], []); excl = []
        else:
            prov = PROV_SLUG.get(r.label, r.label.lower().replace(" ", "-"))
            keys = [(prov, "any")] + ([("shariff-kabunsuan", "any")] if r.label == "Maguindanao" else [])
            polys = sum([load_polygons(idx[k]) for k in keys], [])
            # remove highly urbanized cities / Davao Occidental that sit inside this province polygon
            excl = [load_polygons(idx[k]) for k in list(HUC.values()) + DAVAO_OCC if k[0] == prov]
        m = points_inside(pts, polys)
        for E in excl:
            ii = np.flatnonzero(m)
            if ii.size:
                m[ii] &= ~points_inside(pts[ii], E)
        a = alt[m]; ok = a.size > 0
        rows.append(dict(code=r.code, label=r.label, n_cells=int(m.sum()),
                         mean=a.mean() if ok else np.nan, median=np.median(a) if ok else np.nan,
                         p10=np.percentile(a, 10) if ok else np.nan, p90=np.percentile(a, 90) if ok else np.nan,
                         min=a.min() if ok else np.nan, max=a.max() if ok else np.nan,
                         pct_ge500=(a >= 500).mean() * 100 if ok else np.nan,
                         pct_ge1000=(a >= 1000).mean() * 100 if ok else np.nan))
    e = pd.DataFrame(rows)
    # tiny NCR cities have <20 grid cells: use the median of the other NCR units (all <500 m)
    stat = ["mean", "median", "p10", "p90", "min", "max"]
    ncr = e[e.code.between(77, 93) & (e.n_cells >= 20)]
    for c in e[e.code.between(77, 93) & (e.n_cells < 20)].code:
        e.loc[e.code == c, stat] = ncr[stat].median().values
    DERIVED.mkdir(parents=True, exist_ok=True); e.to_csv(ELEV_CSV, index=False)
    return e


def get_elevation(recompute):
    if recompute or not ELEV_CSV.exists():
        print("Computing province elevation from grid + boundaries (a few minutes)...")
        return build_province_elevation()
    return pd.read_csv(ELEV_CSV)


# ----------------------------------------------------------------------------------
# STEP 4. WHO CLASSIFICATION RULES
# WHO 2024 guideline (Tables 2-5): cutoffs by age/sex/pregnancy trimester; Hb is adjusted
# DOWN for elevation (>=500 m) and current smoking before comparing with the cutoff.
# All values in g/L (observed Hb in g/dL x 10).
# ----------------------------------------------------------------------------------
def elevation_adjustment(alt_m):
    """WHO 2024 equation: 0.0056384*alt + 0.0000003*alt^2 (g/L), zero below 500 m."""
    alt_m = np.asarray(alt_m, float)
    return np.where(alt_m < 500, 0.0, 0.0056384 * alt_m + 0.0000003 * alt_m ** 2)


SMOKING_ADJ = 3.0     # g/L, smoker with unknown cigarettes/day (WHO Table 5); PUF has no quantity


def _trimesters(psc):
    psc = np.asarray(psc)
    return (psc >= 1) & (psc <= 3), (psc >= 4) & (psc <= 6), (psc >= 7) & (psc <= 9)


def who2024_cutoff(age, female, psc):
    t1, t2, t3 = _trimesters(psc); age = np.asarray(age)
    return np.select([t1, t2, t3, age < 2, age < 5, age < 12, age < 15, np.asarray(female)],
                     [110, 105, 110, 105, 110, 115, 120, 120], 130).astype(float)


def who2024_severity_bounds(age, psc):
    """Lower bounds (g/L) of 'mild' and 'moderate'; below the moderate bound = severe (WHO Table 3)."""
    t1, t2, t3 = _trimesters(psc); age = np.asarray(age)
    mild = np.select([t1, t2, t3, age < 2, age < 5], [100, 95, 100, 95, 100], 110).astype(float)
    mod = np.where(t1 | t2 | t3 | (age < 5), 70, 80).astype(float)
    return mild, mod


def who2011_cutoff(age, female, psc):
    """WHO 2011 comparator (no altitude/smoking adjustment): all <5 y <110; pregnant <110."""
    psc, age = np.asarray(psc), np.asarray(age)
    return np.select([(psc >= 1) & (psc <= 9), age < 5, age < 12, age < 15, np.asarray(female)],
                     [110, 110, 115, 120, 120], 130).astype(float)


def classify(d, alt):
    """Add WHO 2024 anemia (0/1), severity (0 none,1 mild,2 moderate,3 severe) and WHO 2011 anemia."""
    hb = d["hemoglobin"].values * 10
    adj_hb = hb - elevation_adjustment(alt) - SMOKING_ADJ * d["smoker"].values
    cut = who2024_cutoff(d.age, d.female, d.psc)
    mild, mod = who2024_severity_bounds(d.age, d.psc)
    d = d.copy()
    d["anemia24"] = (adj_hb < cut).astype(int)
    d["severity"] = np.select([adj_hb >= cut, adj_hb >= mild, adj_hb >= mod], [0, 1, 2], 3)
    d["anemia11"] = (hb < who2011_cutoff(d.age, d.female, d.psc)).astype(int)
    return d


# ----------------------------------------------------------------------------------
# STEP 5. POPULATION GROUPS (as pre-specified)
# Adolescents/adults exclude pregnant, lactating and females with psc != 0.
# ----------------------------------------------------------------------------------
def groups(d):
    a, f, p = d.age.values, d.female.values, d.psc.values
    return {"Overall": np.ones(len(d), bool),
            "Under 5 children": (a >= .5) & (a < 5), "School-age children": (a >= 5) & (a < 10),
            "Adolescents": (a >= 10) & (a < 20) & (~f | (p == 0)), "Adults": (a >= 20) & (a < 60) & (~f | (p == 0)),
            "Older persons": a >= 60, "NPNL women": f & (p == 0) & (a >= 15) & (a < 50),
            "Pregnant women": (p >= 1) & (p <= 9), "First trimester": (p >= 1) & (p <= 3),
            "Second trimester": (p >= 4) & (p <= 6), "Third trimester": (p >= 7) & (p <= 9),
            "Lactating mothers": (p >= 10) & (p <= 12)}


# ----------------------------------------------------------------------------------
# STEP 6. SURVEY ESTIMATION
# Weighted proportion (ratio estimator) with Taylor-linearised variance: PSUs clustered
# within strata; subgroup (domain) estimation keeps the full design; logit 95% CI.
# ----------------------------------------------------------------------------------
class Design:
    def __init__(self, d):
        self.w = d[WEIGHT].values
        self.psu = pd.factorize(d.provhuc * 10000 + d.ms_psucode)[0]
        self.n_psu = self.psu.max() + 1
        stratum = pd.factorize(d.provhuc)[0]
        s_of_psu = np.zeros(self.n_psu, int); s_of_psu[self.psu] = stratum
        self.order = np.argsort(s_of_psu, kind="stable")
        self.bounds = np.r_[0, np.cumsum(np.bincount(s_of_psu))]
        self.df = self.n_psu - (stratum.max() + 1)
        self.tq = t_dist.ppf(.975, self.df)

    def prop(self, y, dom):
        """Return n, cases, % and 95% CI limits for indicator y within boolean domain."""
        y = np.asarray(y, float); dom = np.asarray(dom, bool); n = int(dom.sum())
        if n == 0:
            return dict(n=0, cases=0, pct=np.nan, lo=np.nan, hi=np.nan)
        wd = self.w * dom; W = wd.sum(); p = (wd * y).sum() / W
        z = np.bincount(self.psu, weights=wd * (y - p) / W, minlength=self.n_psu)[self.order]
        var = 0.0
        for s in range(len(self.bounds) - 1):
            g = z[self.bounds[s]:self.bounds[s + 1]]
            if len(g) > 1:
                var += len(g) / (len(g) - 1) * ((g - g.mean()) ** 2).sum()
        se = np.sqrt(var)
        if 0 < p < 1:
            sl = se / (p * (1 - p)); lo, hi = expit(logit(p) - self.tq * sl), expit(logit(p) + self.tq * sl)
        else:                       # zero/all cases: exact bound that ignores the design effect
            lo, hi = (0, 1 - 0.025 ** (1 / n)) if p == 0 else (0.025 ** (1 / n), 1)
        return dict(n=n, cases=int(y[dom].sum()), pct=p * 100, lo=lo * 100, hi=hi * 100)


def fmt(r):
    return "—" if np.isnan(r["pct"]) else f"{r['pct']:.1f} ({r['lo']:.1f}–{r['hi']:.1f})"


# ----------------------------------------------------------------------------------
# STEP 7. TABLES
# ----------------------------------------------------------------------------------
def table1_and_s1(d, G, des):
    t1, s1 = [], []
    for g, m in G.items():
        r24, r11 = des.prop(d.anemia24, m), des.prop(d.anemia11, m)
        t1.append({"Population": g, "n": r24["n"], "Anemia cases (n)": r24["cases"],
                   "Prevalence % (95% CI), WHO 2024": fmt(r24)})
        s1.append({"Population": g, "WHO 2011 % (95% CI), unadjusted": fmt(r11), "WHO 2011 %": round(r11["pct"], 1),
                   "WHO 2024 %": round(r24["pct"], 1), "Difference (pp)": round(r24["pct"] - r11["pct"], 1)})
    return pd.DataFrame(t1), pd.DataFrame(s1)


def table2(d, G, des):
    rows = []
    for g in ["Overall", "Under 5 children", "School-age children", "Adolescents", "Adults", "Older persons"]:
        for sex, flag in (("Male", ~d.female.values), ("Female", d.female.values)):
            r = des.prop(d.anemia24, G[g] & flag)
            rows.append({"Age group": g, "Sex": sex, "n": r["n"], "Anemia cases (n)": r["cases"],
                         "Prevalence % (95% CI)": fmt(r)})
    return pd.DataFrame(rows)


def table3(d, G, des):
    rows = []
    for g, m in G.items():
        row = {"Population": g}
        for name, y in (("Any anemia", d.severity > 0), ("Mild", d.severity == 1),
                        ("Moderate", d.severity == 2), ("Severe", d.severity == 3)):
            row[f"{name} %"] = round(des.prop(y, m)["pct"], 1)
        rows.append(row)
    return pd.DataFrame(rows)


def sensitivity(d0, e, G, des):
    """Prevalence under alternative elevation assumptions (one province statistic for everyone)."""
    m = d0.merge(e[["code", "mean", "median", "p10", "p90", "min", "max"]].rename(columns={"code": "provhuc"}),
                 on="provhuc", how="left")
    scen = {"Sea level (no elevation adjustment)": np.zeros(len(d0)), "Province minimum": m["min"],
            "Province 10th percentile": m.p10, "Province median (PRIMARY)": m["median"], "Province mean": m["mean"],
            "Province 90th percentile": m.p90, "Province maximum": m["max"]}
    out = {}
    for name, alt in scen.items():
        c = classify(d0, np.asarray(alt, float))
        out[name] = {g: round(des.prop(c.anemia24, gm)["pct"], 1) for g, gm in G.items()}
    return pd.DataFrame(out).T


def subnational(d, G, des):
    """Prevalence by region and province/HUC for overall, sex and population groups."""
    base = {k: G[k] for k in ["Overall", "Under 5 children", "School-age children", "Adolescents", "Adults",
                              "Older persons", "NPNL women", "Pregnant women", "Lactating mothers"]}
    base["Male"], base["Female"] = ~d.female.values, d.female.values
    lab = province_labels().set_index("code")["label"]
    units = [("Region", rc, name, (d.regcode == rc).values) for rc, name in REGIONS.items()] + \
            [("Province/HUC", pc, lab[pc], (d.provhuc == pc).values) for pc in sorted(d.provhuc.unique())]
    rows = []
    for level, code, name, um in units:
        for g, gm in base.items():
            r = des.prop(d.anemia24, um & gm)
            if r["n"] < 30:                      # do not publish unstable cell estimates
                r.update(pct=np.nan, lo=np.nan, hi=np.nan)
            flag = "suppressed (n<30)" if r["n"] < 30 else ("caution (n 30-49)" if r["n"] < 50 else "")
            rows.append(dict(level=level, code=code, unit=name, group=g, **r, reliability=flag))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recompute-elevation", action="store_true",
                    help="rebuild province elevation from grid + boundaries")
    ap.add_argument("--sensitivity", action="store_true")
    ap.add_argument("--subnational", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)

    raw = load_data()
    d, flow = eligibility(raw)
    flow.to_csv(OUT / "participant_flow.csv", index=False); print(flow.to_string(index=False))
    eda(d)

    e = get_elevation(args.recompute_elevation)
    d = d.merge(e[["code", "median"]].rename(columns={"code": "provhuc", "median": "elev_m"}), on="provhuc", how="left")
    assert d.elev_m.notna().all(), "some provinces lack elevation"
    d0 = d.drop(columns="elev_m")
    d = classify(d, d.elev_m.values)

    G, des = groups(d), Design(d)
    print(f"Design: {des.n_psu} PSUs, {len(des.bounds) - 1} strata, df = {des.df}")
    T1, S1 = table1_and_s1(d, G, des)
    T2, T3 = table2(d, G, des), table3(d, G, des)
    with pd.ExcelWriter(OUT / "tables_1-3_S1.xlsx") as xw:
        T1.to_excel(xw, sheet_name="Table 1", index=False); T2.to_excel(xw, sheet_name="Table 2", index=False)
        T3.to_excel(xw, sheet_name="Table 3", index=False); S1.to_excel(xw, sheet_name="Table S1", index=False)
        e.to_excel(xw, sheet_name="Province elevation", index=False)
        if args.sensitivity:
            sensitivity(d0, e, G, des).to_excel(xw, sheet_name="Sensitivity")
    for n, t in (("table1", T1), ("table2", T2), ("table3", T3), ("tableS1", S1)):
        t.to_csv(OUT / f"{n}.csv", index=False)
    if args.subnational:
        subnational(d, G, des).to_csv(OUT / "subnational_long.csv", index=False)
    print(T1.to_string(index=False)); print("Outputs in", OUT)


if __name__ == "__main__":
    main()
