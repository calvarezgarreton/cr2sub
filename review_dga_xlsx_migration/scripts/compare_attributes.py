#!/usr/bin/env python3
"""Extended v1.1 vs v1.2 comparison of attributes.csv, with focus on CLSoilMaps.

Reads ../versions/{v1.1,v1.2}/*_attributes.csv and (optionally)
../comparison_v1.1_vs_v1.2/clsoilmaps_reextraction.csv (made by
validate_clsoilmaps_extraction.R). Writes CSVs/PNGs to ../comparison_v1.1_vs_v1.2/
and prints a summary.

Usage (needs pandas, numpy, matplotlib):
    python3 compare_attributes.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "comparison_v1.1_vs_v1.2"
old = pd.read_csv(ROOT / "versions/v1.1/cr2sub_v1.1_attributes.csv").set_index("cr2sub_id")
new = pd.read_csv(ROOT / "versions/v1.2/cr2sub_v1.2_attributes.csv").set_index("cr2sub_id")
common = old.index.intersection(new.index)
o, n = old.loc[common], new.loc[common]
SOIL = [c for c in new.columns if "clsoilmap" in c]

# 1. Column by column ---------------------------------------------------------
rows = []
for c in new.columns:
    a, b = o[c], n[c]
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        both = a.notna() & b.notna()
        d = (a - b).abs()[both]
        differ = int((d > 1e-6).sum())
        rows.append(dict(column=c, compared=int(both.sum()), differ=differ,
                         na_to_value=int((a.isna() & b.notna()).sum()),
                         value_to_na=int((a.notna() & b.isna()).sum()),
                         median_abs_diff=float(d[d > 1e-6].median()) if differ else 0.0,
                         max_abs_diff=float(d.max()) if len(d) else 0.0))
    else:
        both = a.notna() & b.notna()
        rows.append(dict(column=c, compared=int(both.sum()), differ=int((a[both] != b[both]).sum()),
                         na_to_value=int((a.isna() & b.notna()).sum()),
                         value_to_na=int((a.notna() & b.isna()).sum()),
                         median_abs_diff=np.nan, max_abs_diff=np.nan))
col_diffs = pd.DataFrame(rows)
col_diffs.to_csv(OUT / "attribute_column_diffs.csv", index=False)
print(f"wells: v1.1={len(old)} v1.2={len(new)} common={len(common)}")
print("\n== Columns that differ (common wells) ==")
print(col_diffs[col_diffs.differ + col_diffs.na_to_value + col_diffs.value_to_na > 0]
      .to_string(index=False, float_format=lambda x: f"{x:.4g}"))
print("columns identical:", int((col_diffs.differ + col_diffs.na_to_value + col_diffs.value_to_na == 0).sum()),
      "of", len(col_diffs))

# 2. Basin assignment ---------------------------------------------------------
print("\n== Basin assignment ==")
for c, label in (("cr2sub_in_basin_camels", "CAMELS (v1 -> v2)"), ("cr2sub_in_basin_bna", "BNA")):
    a, b = o[c], n[c]
    print(f"{label}: same={int((a == b).sum())}  changed={int((a.notna() & b.notna() & (a != b)).sum())}  "
          f"NA->id={int((a.isna() & b.notna()).sum())}  id->NA={int((a.notna() & b.isna()).sum())}")

# 3. CLSoilMaps: v1.1 vs v1.2 -----------------------------------------------------
print("\n== CLSoilMaps v1.1 vs v1.2 (common wells) ==")
so, sn = o[SOIL], n[SOIL]
same = (so.eq(sn) | (so.isna() & sn.isna())).all(axis=1)
print(f"wells with identical soil values: {int(same.sum())} of {len(common)}")
print("wells with soil NA: v1.1=%d v1.2=%d (all wells: %d)" %
      (so.isna().any(axis=1).sum(), sn.isna().any(axis=1).sum(), new[SOIL].isna().any(axis=1).sum()))
stats = pd.DataFrame({"min": new[SOIL].min(), "median": new[SOIL].median(),
                      "max": new[SOIL].max(), "n_na": new[SOIL].isna().sum()}).round(3)
stats.to_csv(OUT / "clsoilmaps_summary_v1.2.csv")
print(stats.to_string())

# physical consistency
chk = {
    "clay+sand > 100 (0-100)": int(((new.cr2sub_clsoilmap_clay_0_100cm + new.cr2sub_clsoilmap_sand_0_100cm) > 100).sum()),
    "clay+sand > 100 (100-200)": int(((new.cr2sub_clsoilmap_clay_100_200cm + new.cr2sub_clsoilmap_sand_100_200cm) > 100).sum()),
    "negative values": int((new[SOIL] < 0).any(axis=1).sum()),
    "bulk density outside 0.5-2.0": int(((new.filter(like="bulkd") < 0.5) | (new.filter(like="bulkd") > 2.0)).any(axis=1).sum()),
}
print("\nphysical checks (wells failing):", chk)

# wells whose coordinates changed but soil was not recomputed
moved = col_diffs.set_index("column").loc["dga_well_lat", "differ"]
chg = common[((o.dga_well_lat - n.dga_well_lat).abs() > 0.01) | ((o.dga_well_lon - n.dga_well_lon).abs() > 0.01)]
print(f"\nwells with changed coordinates (>0.01 deg): {len(chg)} ->", list(chg))
print("  their soil values are identical v1.1 vs v1.2 (not recomputed from new location):",
      bool(same.loc[chg].all()) if len(chg) else "n/a")

# 4. Spatial coherence between neighbouring wells ---------------------------------
lat = np.radians(new.cr2sub_lat.values)
lon = np.radians(new.cr2sub_lon.values)
dlat = lat[:, None] - lat[None, :]
dlon = lon[:, None] - lon[None, :]
h = np.sin(dlat / 2) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2
dist = 2 * 6371 * np.arcsin(np.sqrt(h))
np.fill_diagonal(dist, np.inf)
nn = dist.argmin(axis=1)
nn_d = dist[np.arange(len(new)), nn]
vals = new[SOIL].values
coh = []
for k, c in enumerate(SOIL):
    ok = ~np.isnan(vals[:, k]) & ~np.isnan(vals[nn, k]) & (nn_d < 5)
    rnd = np.random.default_rng(0).permutation(len(new))
    ok_r = ~np.isnan(vals[:, k]) & ~np.isnan(vals[rnd, k])
    coh.append(dict(variable=c, pairs_within_5km=int(ok.sum()),
                    corr_nearest=float(np.corrcoef(vals[ok, k], vals[nn, k][ok])[0, 1]),
                    corr_random_pairs=float(np.corrcoef(vals[ok_r, k], vals[rnd, k][ok_r])[0, 1])))
coh = pd.DataFrame(coh)
coh.to_csv(OUT / "clsoilmaps_neighbour_coherence.csv", index=False)
print("\n== Spatial coherence: soil value vs nearest well (<5 km) ==")
print(coh.round(3).to_string(index=False))

print("\n== New wells 05120007 / 05120011 vs their 3 nearest wells ==")
for wid in (5120007, 5120011):
    i = new.index.get_loc(wid)
    order = np.argsort(dist[i])[:3]
    cols = ["cr2sub_clsoilmap_awc_0_100cm", "cr2sub_clsoilmap_clay_0_100cm", "cr2sub_clsoilmap_ksat_0_100cm"]
    t = new.iloc[[i, *order]][cols].copy()
    t.insert(0, "dist_km", [0.0, *np.round(dist[i][order], 1)])
    print(f"well {wid}:"); print(t.round(2).to_string())

# 5. Re-extraction vs stored values ---------------------------------------------
rx = OUT / "clsoilmaps_reextraction.csv"
if rx.exists():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    re = pd.read_csv(rx).set_index("cr2sub_id")
    ids = new.index.intersection(re.index)
    res = []
    fig, axs = plt.subplots(2, 5, figsize=(18, 7))
    for ax, c in zip(axs.ravel(), SOIL):
        a, b = new.loc[ids, c], re.loc[ids, c]
        ok = a.notna() & b.notna()
        d = (a - b).abs()[ok]
        rel = (d / a[ok].abs().clip(lower=1e-9))
        res.append(dict(variable=c, n=int(ok.sum()), n_na_stored=int(a.isna().sum()),
                        n_na_reextracted=int(b.isna().sum()),
                        corr=float(np.corrcoef(a[ok], b[ok])[0, 1]),
                        median_abs_diff=float(d.median()), p95_abs_diff=float(d.quantile(.95)),
                        max_abs_diff=float(d.max()), median_rel_diff_pct=float(rel.median() * 100),
                        n_identical=int((d < 1e-3).sum())))
        ax.scatter(a[ok], b[ok], s=4, alpha=.5)
        lim = [min(a[ok].min(), b[ok].min()), max(a[ok].max(), b[ok].max())]
        ax.plot(lim, lim, "k--", lw=.8)
        ax.set_title(c.replace("cr2sub_clsoilmap_", ""), fontsize=9)
        ax.set_xlabel("stored (v1.1/v1.2 CSV)", fontsize=7); ax.set_ylabel("re-extracted", fontsize=7)
    fig.suptitle("CLSoilMaps: stored values vs re-extraction with scripts/functions/functions_clsoilmaps.R")
    fig.tight_layout(); fig.savefig(OUT / "clsoilmaps_stored_vs_reextracted.png", dpi=110)
    res = pd.DataFrame(res)
    res.to_csv(OUT / "clsoilmaps_reextraction_diffs.csv", index=False)
    print("\n== Re-extraction vs stored values ==")
    print(res.round(3).to_string(index=False))
    worst = (re.loc[ids, SOIL] - new.loc[ids, SOIL]).abs().div(new[SOIL].std()).max(axis=1).sort_values(ascending=False).head(10)
    print("\n10 wells with the largest standardized difference:")
    print(worst.round(2).to_string())
    na_diff = ids[(new.loc[ids, SOIL].isna().any(axis=1)) != (re.loc[ids, SOIL].isna().any(axis=1))]
    print("wells where NA status differs (stored vs re-extracted):", len(na_diff), list(na_diff)[:10])
