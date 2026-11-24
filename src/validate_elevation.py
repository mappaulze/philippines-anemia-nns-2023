#!/usr/bin/env python3
"""
OPTIONAL: validate the 1 km province elevation against an independent 30 m DSM (JAXA ALOS tiles).

For selected provinces this draws random points inside the province polygon, reads the 30 m DSM
at those points and compares mean/median/%>=500 m with data/derived/province_elevation.csv.
Expected input: GeoTIFF tiles named like ALPSMLC30_N016E120_DSM.tif placed anywhere under
data/raw/jaxa/ (unzip the JAXA downloads there). Requires: tifffile.

    python src/validate_elevation.py --provinces benguet mountain-province ifugao kalinga abra nueva-vizcaya quirino
"""
import argparse
import numpy as np
import pandas as pd
import tifffile
import nns_anemia_who2024 as m

ap = argparse.ArgumentParser(); ap.add_argument("--provinces", nargs="+", required=True)
ap.add_argument("--npoints", type=int, default=15000); args = ap.parse_args()
rng = np.random.default_rng(1)                      # fixed seed -> reproducible
idx = m.index_geojson(); tiles_dir = m.RAW / "jaxa"
files = {f.name: f for f in tiles_dir.rglob("ALPSMLC30_*_DSM.tif")}
cache = {}


def tile(lat, lon):
    k = f"ALPSMLC30_N{lat:03d}E{lon:03d}_DSM.tif"
    if k not in cache:
        cache[k] = tifffile.imread(files[k]) if k in files else None
    return cache[k]


e = pd.read_csv(m.ELEV_CSV).set_index("label"); rows = []
lab = {p.lower().replace(" ", "-"): p for p in e.index}
for prov in args.provinces:
    polys = m.load_polygons(idx[(prov, "any")])
    allp = np.vstack([r[0] for r in polys]); x0, y0 = allp.min(0); x1, y1 = allp.max(0)
    pts = np.c_[rng.uniform(x0, x1, 8 * args.npoints), rng.uniform(y0, y1, 8 * args.npoints)]
    pts = pts[m.points_inside(pts, polys)][: args.npoints]
    la, lo = np.floor(pts[:, 1]).astype(int), np.floor(pts[:, 0]).astype(int); v = np.full(len(pts), np.nan)
    for a, b in set(zip(la, lo)):
        t = tile(a, b)
        if t is None:
            continue
        sel = (la == a) & (lo == b); n = t.shape[0]
        r = np.clip(((a + 1 - pts[sel, 1]) * n).astype(int), 0, n - 1); c = np.clip(((pts[sel, 0] - b) * n).astype(int), 0, n - 1)
        v[sel] = t[r, c]
    ok = np.isfinite(v) & (v > -1000)
    g = e.loc[lab[prov]] if prov in lab else None
    rows.append(dict(province=prov, points_with_tile=int(ok.sum()), points_requested=len(pts), jaxa_mean=v[ok].mean(),
                     grid_mean=None if g is None else g["mean"], jaxa_median=np.median(v[ok]),
                     grid_median=None if g is None else g["median"], jaxa_pct_ge500=(v[ok] >= 500).mean() * 100,
                     grid_pct_ge500=None if g is None else g["pct_ge500"]))
out = pd.DataFrame(rows); out["mean_diff_m"] = out.jaxa_mean - out.grid_mean
out.round(1).to_csv(m.OUT / "elevation_validation.csv", index=False); print(out.round(1).to_string(index=False))
print("NOTE: where points_with_tile << points_requested, tiles are missing and the comparison is partial.")
