#!/usr/bin/env python3
"""Compare CR2SUB v1.2 against the DGA source (AM006T0034616_1.xlsx) it was built from.

Independent re-computation in Python (the pipeline is R): monthly series, coordinates,
UTM/lat-lon consistency, elevation and basin code. Only v1.2 is used.
Writes fig_g_vs_dga.png and vs_dga_*.csv to ../comparison_v1.1_vs_v1.2/.
Needs pandas, numpy, matplotlib, openpyxl, pyproj.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Transformer

ROOT = Path(__file__).resolve().parent.parent
REPO = ROOT.parent
OUT = ROOT / "comparison_v1.1_vs_v1.2"
BLUE, ORANGE, GREY = "#1f4e79", "#e07b00", "#999999"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": .25})

raw = pd.read_excel(REPO / "input/DGA_GWL_obs_descarga/AM006T0034616_1.xlsx", sheet_name=0, dtype=str)
raw.columns = ["name", "code", "prov", "comm", "basin", "subbasin", "subsub", "utm_n", "utm_e",
               "lat_raw", "lon_raw", "elev", "fecha", "valor"]
for c in raw.columns:
    raw[c] = raw[c].str.strip()
att = pd.read_csv(ROOT / "versions/v1.2/cr2sub_v1.2_attributes.csv").set_index("cr2sub_id")
gwl = pd.read_csv(ROOT / "versions/v1.2/cr2sub_v1.2_gwl_mon.csv", index_col="date", parse_dates=True)
gwl.columns = gwl.columns.astype(int)

# 1. Wells ---------------------------------------------------------------------------
n_all = len(raw)
raw = raw[raw.code.notna()].copy()
print(f"Filas del export: {n_all}; sin código de estación: {n_all - len(raw)}")
raw["id"] = raw.code.str.replace(r"-[0-9Kk]$", "", regex=True).astype(int)
src_wells = set(raw.id.unique())
print(f"Pozos en el export DGA: {len(src_wells)}; en attributes v1.2: {len(att)}; "
      f"en el export y no en v1.2: {sorted(src_wells - set(att.index))}; en v1.2 y no en el export: {sorted(set(att.index) - src_wells)}")

# 2. Monthly series re-computed from the raw records -------------------------------------
raw["fecha"] = pd.to_datetime(raw.fecha, format="%d-%m-%Y", errors="coerce")
raw["valor"] = pd.to_numeric(raw.valor, errors="coerce")
n_all = len(raw)
valid = raw.dropna(subset=["code", "fecha", "valor"])
print(f"Registros: {n_all}; válidos (código, fecha y valor): {len(valid)}; con valor <= 0: {int((valid.valor <= 0).sum())}")
daily = valid.groupby(["id", "fecha"]).valor.mean().reset_index()
daily = daily[daily.valor > 0]
daily["mes"] = daily.fecha.values.astype("datetime64[M]")
monthly = (-daily.groupby(["id", "mes"]).valor.mean()).unstack(0)
monthly = monthly.reindex(columns=[c for c in monthly.columns if c in gwl.columns])
common_cols = monthly.columns
a = monthly.reindex(gwl.index)[common_cols]
b = gwl[common_cols]
both = a.notna() & b.notna()
diff = (a - b).abs().where(both)
only_src = int((a.notna() & b.isna()).sum().sum())
only_csv = int((a.isna() & b.notna()).sum().sum())
mism = int((diff > 1e-6).sum().sum())
print(f"Series: celdas pozo-mes con dato en ambos={int(both.sum().sum())}; distintas (>1e-6 m)={mism}; "
      f"solo en el cálculo independiente={only_src}; solo en gwl_mon.csv={only_csv}; dif. máx={float(diff.max().max()):.2g} m")
st = diff.stack(); st[st > 1e-6].rename("abs_diff").to_csv(OUT / "vs_dga_series_mismatches.csv")

# 3. Coordinates: DGA raw string -> degrees, UTM -> lat/lon ----------------------------------
def dms(s):
    s = str(s).strip()
    return int(s[:3]) + int(s[3:5]) / 60 + int(s[5:7]) / 3600
w = raw.drop_duplicates("id").set_index("id")
lat_i = -w.lat_raw.map(dms); lon_i = -w.lon_raw.map(dms)
w = w.loc[w.index.intersection(att.index)]
dl = pd.concat([(lat_i - att.dga_well_lat).abs(), (lon_i - att.dga_well_lon).abs()], axis=1).dropna()
dl.columns = ["dlat", "dlon"]
print(f"Coordenadas: DGA (DDDMMSS) convertidas de forma independiente vs dga_well_lat/lon: máx |dlat|={dl.dlat.max():.2g}°, máx |dlon|={dl.dlon.max():.2g}°")
print(f"cr2sub_lat/lon distintas de dga_well_lat/lon: {int(((att.cr2sub_lat != att.dga_well_lat) | (att.cr2sub_lon != att.dga_well_lon)).sum())} pozos")
res = {}
for z in (18, 19):
    t = Transformer.from_crs(f"EPSG:327{z}", "EPSG:4326", always_xy=True)
    e = pd.to_numeric(w.utm_e, errors="coerce"); nn = pd.to_numeric(w.utm_n, errors="coerce")
    lo, la = t.transform(e.values, nn.values)
    dy = (la - att.loc[w.index, "dga_well_lat"].values) * 111.2
    dx = (lo - att.loc[w.index, "dga_well_lon"].values) * 111.2 * np.cos(np.radians(la))
    res[z] = pd.Series(np.hypot(dx, dy), index=w.index)
best = pd.concat(res, axis=1).min(axis=1)
zone = pd.concat(res, axis=1).idxmin(axis=1)
print(f"UTM vs lat/lon (mejor huso 18 o 19): mediana={best.median():.3f} km; >1 km: {int((best > 1).sum())}; >10 km: {int((best > 10).sum())}; huso 18 es el mejor en {int((zone == 18).sum())} pozos")
best.rename("km").loc[best > 1].to_csv(OUT / "vs_dga_utm_latlon_gt1km.csv")

# 4. Elevation: DGA vs DEM ---------------------------------------------------------------------
z0 = att.dga_well_elev == 0
d_el = (att.cr2sub_elev - att.dga_well_elev).abs()
print(f"Cota: DGA=0 en {int(z0.sum())} pozos; sin ellos, mediana |DEM-DGA|={d_el[~z0].median():.1f} m; >200 m: {int((d_el[~z0] > 200).sum())}")
d_el[(d_el > 200) & ~z0].rename("abs_diff_m").to_csv(OUT / "vs_dga_elev_gt200m.csv")

# 5. Basin: BNA polygon vs DGA code prefix ----------------------------------------------------
pref = pd.Series([int(str(i).zfill(8)[:3]) for i in att.index], index=att.index)
agree = pref == att.cr2sub_in_basin_bna
print(f"Cuenca BNA (polígono) vs prefijo del código DGA: coinciden {int(agree.sum())} de {len(att)}")
att.loc[~agree, ["dga_well_name", "dga_well_basin", "cr2sub_in_basin_bna"]].assign(prefijo_codigo=pref[~agree]).to_csv(OUT / "vs_dga_bna_prefix_mismatch.csv")

# Figure ----------------------------------------------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(13, 10))
ax = axs[0, 0]
sa, sb = a.stack(), b.stack(); k = sa.index.intersection(sb.dropna().index)
ax.scatter(sa[k], sb[k], s=1, c=BLUE, alpha=.3)
lim = [min(sa[k].min(), sb[k].min()), max(sa[k].max(), sb[k].max())]; ax.plot(lim, lim, "k--", lw=.6)
ax.set_xlabel("recalculado desde el export DGA [m]"); ax.set_ylabel("gwl_mon.csv v1.2 [m]")
ax.set_title(f"Serie mensual: {len(k):,} pozo-mes; distintos: {mism}; solo en uno de los dos: {only_src + only_csv}")
ax = axs[0, 1]
ax.hist(np.clip(best, 1e-3, None), bins=np.logspace(-3, 3, 50), color=BLUE); ax.set_xscale("log")
ax.set_xlabel("distancia entre UTM y lat/lon reportados [km]"); ax.set_ylabel("pozos")
ax.set_title(f"Coordenadas DGA: UTM vs lat/lon (>1 km: {int((best > 1).sum())} pozos)")
ax = axs[1, 0]
ax.scatter(att.dga_well_elev[~z0], att.cr2sub_elev[~z0], s=6, c=BLUE, alpha=.6, label=f"cota DGA > 0 ({int((~z0).sum())})")
ax.scatter(att.dga_well_elev[z0], att.cr2sub_elev[z0], s=6, c=ORANGE, alpha=.6, label=f"cota DGA = 0 ({int(z0.sum())})")
lim = [0, max(att.dga_well_elev.max(), att.cr2sub_elev.max())]; ax.plot(lim, lim, "k--", lw=.6)
ax.set_xlabel("cota DGA [m]"); ax.set_ylabel("cota DEM en la coordenada [m]"); ax.legend(); ax.set_title("Cota: DGA vs DEM (v1.2)")
ax = axs[1, 1]
ax.bar(["coincide", "no coincide"], [int(agree.sum()), int((~agree).sum())], color=[BLUE, ORANGE])
for i, v in enumerate([int(agree.sum()), int((~agree).sum())]): ax.text(i, v, str(v), ha="center", va="bottom")
ax.set_title("Cuenca BNA del polígono vs prefijo (3 dígitos) del código DGA")
fig.suptitle("G. CR2SUB v1.2 vs la fuente DGA (AM006T0034616_1.xlsx)")
fig.tight_layout(); fig.savefig(OUT / "fig_g_vs_dga.png", dpi=110)
