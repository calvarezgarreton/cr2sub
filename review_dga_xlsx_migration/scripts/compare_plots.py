#!/usr/bin/env python3
"""Verification figures v1.1 vs v1.2 for every attribute group (not only CLSoilMaps).

Reads ../versions/{v1.1,v1.2}/ and writes fig_*.png (+ a console summary) into
../comparison_v1.1_vs_v1.2/. Needs pandas, numpy, matplotlib.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "comparison_v1.1_vs_v1.2"
V1, V2 = ROOT / "versions/v1.1", ROOT / "versions/v1.2"
old = pd.read_csv(V1 / "cr2sub_v1.1_attributes.csv").set_index("cr2sub_id")
new = pd.read_csv(V2 / "cr2sub_v1.2_attributes.csv").set_index("cr2sub_id")
common = old.index.intersection(new.index)
o, n = old.loc[common], new.loc[common]
BLUE, ORANGE, GREY = "#1f4e79", "#e07b00", "#999999"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": .25})


def one_to_one(ax, a, b, title, mask=None, log=False, unit=""):
    ok = a.notna() & b.notna()
    same = ok & ((a - b).abs() < 1e-6)
    diff = ok & ~same
    ax.scatter(a[same], b[same], s=5, c=GREY, alpha=.5, label=f"igual ({same.sum()})")
    ax.scatter(a[diff], b[diff], s=10, c=ORANGE, alpha=.8, label=f"cambia ({diff.sum()})")
    both = pd.concat([a[ok], b[ok]])
    if len(both):
        lo, hi = both.min(), both.max()
        ax.plot([lo, hi], [lo, hi], "k--", lw=.7)
    if log:
        ax.set_xscale("symlog"); ax.set_yscale("symlog")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(f"v1.1 {unit}"); ax.set_ylabel(f"v1.2 {unit}")
    ax.legend(fontsize=7, loc="upper left")


# A. Groundwater statistics -------------------------------------------------------
fig, axs = plt.subplots(2, 3, figsize=(15, 9))
for ax, (c, t, u) in zip(axs.ravel(), [
        ("cr2sub_mean_gwl", "Nivel medio (raw)", "[m]"), ("cr2sub_sd_gwl", "Desv. estándar (raw)", "[m]"),
        ("cr2sub_cv_gwl", "CV (raw)", ""), ("cr2sub_clean_mean_gwl", "Nivel medio (clean)", "[m]"),
        ("cr2sub_clean_sd_gwl", "Desv. estándar (clean)", "[m]"), ("cr2sub_clean_cv_gwl", "CV (clean)", "")]):
    one_to_one(ax, o[c], n[c], t, unit=u)
fig.suptitle("A. Estadísticos del nivel freático: v1.1 vs v1.2 (pozos comunes); naranja = cambia por más meses de datos")
fig.tight_layout(); fig.savefig(OUT / "fig_a_gwl_stats.png", dpi=110); plt.close(fig)

# B. Series: observations per well, wells per year, outliers removed -------------
def load(v, stem):
    return pd.read_csv(V1 / f"cr2sub_v1.1_{stem}.csv" if v == 1 else V2 / f"cr2sub_v1.2_{stem}.csv",
                       index_col="date", parse_dates=True)
g1, g2 = load(1, "gwl_mon"), load(2, "gwl_mon")
c1, c2 = load(1, "gwl_mon_clean"), load(2, "gwl_mon_clean")
nobs1, nobs2 = g1.notna().sum(), g2.notna().sum()
cn1, cn2 = c1.notna().sum(), c2.notna().sum()
fig = plt.figure(figsize=(15, 8.5))
gs = fig.add_gridspec(2, 6)
axs = {(0, 0): fig.add_subplot(gs[0, 0:2]), (0, 1): fig.add_subplot(gs[0, 2:4]), (0, 2): fig.add_subplot(gs[0, 4:6]),
       (1, 1): fig.add_subplot(gs[1, 0:3]), (1, 2): fig.add_subplot(gs[1, 3:6])}
ax = axs[0, 0]
w = nobs1.index.intersection(nobs2.index)
ax.scatter(nobs1[w], nobs2[w], s=6, c=BLUE, alpha=.6)
mx = max(nobs1.max(), nobs2.max()); ax.plot([0, mx], [0, mx], "k--", lw=.7)
ax.set_xlabel("meses con dato, v1.1"); ax.set_ylabel("meses con dato, v1.2"); ax.set_title("Meses con dato por pozo")
ax = axs[0, 1]
d = (nobs2[w] - nobs1[w]); ax.hist(d[d > 0], bins=40, color=ORANGE)
ax.set_title(f"Meses agregados por pozo ({int((d > 0).sum())} pozos ganan datos, {int((d < 0).sum())} pierden)")
ax.set_xlabel("meses agregados"); ax.set_ylabel("pozos"); ax.set_yscale("log")
ax = axs[0, 2]
ax.plot(g1.notna().sum(axis=1).rolling(12, center=True).mean(), color=ORANGE, lw=1.2, label="v1.1")
ax.plot(g2.notna().sum(axis=1).rolling(12, center=True).mean(), color=BLUE, lw=1.2, label="v1.2")
ax.set_title("Pozos con dato por mes (media móvil 12 m)"); ax.legend(); ax.set_ylabel("pozos")
ax = axs[1, 1]
tot1 = (1 - cn1.sum() / nobs1.sum()) * 100; tot2 = (1 - cn2.sum() / nobs2.sum()) * 100
ax.bar(["v1.1", "v1.2"], [tot1, tot2], color=[ORANGE, BLUE]); ax.set_ylabel("% de datos eliminados")
ax.set_title("Outliers removidos, total"); [ax.text(i, v, f"{v:.1f} %", ha="center", va="bottom") for i, v in enumerate([tot1, tot2])]
ax = axs[1, 2]
yr1 = g1.notna().sum(axis=1).groupby(g1.index.year).sum(); yr2 = g2.notna().sum(axis=1).groupby(g2.index.year).sum()
ax.bar(yr2.index - .2, yr2.values, width=.4, color=BLUE, label="v1.2"); ax.bar(yr1.index + .2, yr1.values, width=.4, color=ORANGE, label="v1.1")
ax.set_title("Observaciones mensuales por año"); ax.legend()
fig.suptitle("B. Series de nivel freático: cantidad de datos y filtro de outliers")
fig.tight_layout(); fig.savefig(OUT / "fig_b_series.png", dpi=110); plt.close(fig)
print(f"B: outliers removidos v1.1={tot1:.2f}% v1.2={tot2:.2f}%; pozos que ganan meses={int((d>0).sum())}, mediana ganada={d[d>0].median()}, max={d.max()}")

# C. Location and elevation ---------------------------------------------------------
mv = common[((o.dga_well_lat - n.dga_well_lat).abs() > 0.01) | ((o.dga_well_lon - n.dga_well_lon).abs() > 0.01)]
fig, axs = plt.subplots(1, 3, figsize=(16, 6))
ax = axs[0]
ax.scatter(new.cr2sub_lon, new.cr2sub_lat, s=3, c=GREY, label="pozos")
for w_ in mv:
    ax.annotate("", xy=(n.loc[w_, "cr2sub_lon"], n.loc[w_, "cr2sub_lat"]), xytext=(o.loc[w_, "cr2sub_lon"], o.loc[w_, "cr2sub_lat"]),
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.2))
ax.scatter(o.loc[mv, "cr2sub_lon"], o.loc[mv, "cr2sub_lat"], s=25, c="k", marker="x", label="v1.1")
ax.scatter(n.loc[mv, "cr2sub_lon"], n.loc[mv, "cr2sub_lat"], s=25, c=ORANGE, label="v1.2")
ax.set_title(f"Pozos que cambian de ubicación ({len(mv)}, >0,01°)"); ax.legend(); ax.set_xlabel("lon"); ax.set_ylabel("lat")
for ax, df, name, col in ((axs[1], o, "v1.1", ORANGE), (axs[2], n, "v1.2", BLUE)):
    ax.scatter(df.dga_well_elev, df.cr2sub_elev, s=6, c=col, alpha=.6)
    lim = [min(df.dga_well_elev.min(), df.cr2sub_elev.min()), max(df.dga_well_elev.max(), df.cr2sub_elev.max())]
    ax.plot(lim, lim, "k--", lw=.7); ax.set_xlabel("cota reportada por DGA [m]"); ax.set_ylabel("cota FABDEM en la coordenada [m]")
    dd = (df.cr2sub_elev - df.dga_well_elev).abs(); ax.set_title(f"{name}: mediana |dif| = {dd.median():.0f} m; >200 m: {(dd > 200).sum()} pozos", fontsize=9)
fig.suptitle("C. Verificación de ubicación: cota reportada vs cota del DEM en la coordenada")
fig.tight_layout(); fig.savefig(OUT / "fig_c_location.png", dpi=110); plt.close(fig)
d1 = (o.cr2sub_elev - o.dga_well_elev).abs(); d2 = (n.cr2sub_elev - n.dga_well_elev).abs()
print(f"C: pozos movidos={len(mv)}; |cota DEM-DGA| >200 m: v1.1={int((d1>200).sum())} v1.2={int((d2>200).sum())}; mediana v1.1={d1.median():.1f} v1.2={d2.median():.1f}")
print("   pozos movidos, dif cota DEM-DGA v1.1 -> v1.2:", {int(w_): (round(float(d1[w_])), round(float(d2[w_]))) for w_ in mv})

# D. CAMELS v1 vs v2 -------------------------------------------------------------------
cols = [("cr2sub_camels_pr_yr", "Precipitación anual [mm]"), ("cr2sub_camels_pet_yr", "ETP anual [mm]"),
        ("cr2sub_camels_aridity", "Aridez (ETP/P)"), ("cr2sub_camels_sf", "Fracción nival"),
        ("cr2sub_camels_elev", "Elevación media [m]"), ("cr2sub_camels_slp", "Pendiente media [°]")]
fig, axs = plt.subplots(2, 3, figsize=(15, 9))
for ax, (c, t) in zip(axs.ravel(), cols):
    one_to_one(ax, o[c], n[c], t)
fig.suptitle("D. Atributos de la cuenca CAMELS: v1 (516 cuencas) vs v2 (860 cuencas), pozos comunes")
fig.tight_layout(); fig.savefig(OUT / "fig_d_camels.png", dpi=110); plt.close(fig)

# E. BNA changed wells + internal consistency ---------------------------------------
fig, axs = plt.subplots(2, 3, figsize=(15, 9))
for ax, (c, t) in zip(axs[0], [("cr2sub_bna_pr_yr", "BNA: precipitación [mm]"), ("cr2sub_bna_elev", "BNA: elevación media [m]"),
                               ("cr2sub_bna_slp", "BNA: pendiente media [°]")]):
    one_to_one(ax, o[c], n[c], t)
for ax, (pre, t) in zip(axs[1][:2], [("camels", "CAMELS"), ("bna", "BNA")]):
    dd = (n[f"cr2sub_{pre}_aridity"] - n[f"cr2sub_{pre}_pet_yr"] / n[f"cr2sub_{pre}_pr_yr"]).abs().dropna()
    ax.hist(dd, bins=30, color=BLUE); ax.set_yscale("log")
    ax.set_title(f"{t} v1.2: |aridez - ETP/P| (máx {dd.max():.3f}, n={len(dd)})"); ax.set_xlabel("diferencia")
ax = axs[1, 2]
ok = n.cr2sub_camels_pr_yr.notna()
sc = ax.scatter(n.cr2sub_lat[ok], n.cr2sub_camels_pr_yr[ok], s=5, c=BLUE, alpha=.5, label="CAMELS v2")
ax.scatter(n.cr2sub_lat, n.cr2sub_bna_pr_yr, s=5, c=ORANGE, alpha=.4, label="BNA")
ax.set_xlabel("latitud"); ax.set_ylabel("precipitación anual [mm]"); ax.legend(); ax.set_title("Precipitación vs latitud (gradiente N-S esperado)")
fig.suptitle("E. Cuencas BNA (9 pozos cambian) y consistencia interna de los atributos climáticos")
fig.tight_layout(); fig.savefig(OUT / "fig_e_bna_consistency.png", dpi=110); plt.close(fig)
for pre in ("camels", "bna"):
    dd = (n[f"cr2sub_{pre}_aridity"] - n[f"cr2sub_{pre}_pet_yr"] / n[f"cr2sub_{pre}_pr_yr"]).abs().dropna()
    print(f"E: {pre} aridez vs PET/P max diff={dd.max():.4f}; sf fuera de [0,1]: {int(((n[f'cr2sub_{pre}_sf']<0)|(n[f'cr2sub_{pre}_sf']>1)).sum())}")

# F. Soil map of values --------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(16, 7), sharey=True)
for ax, (c, t) in zip(axs, [("cr2sub_clsoilmap_awc_0_100cm", "AWC 0-100 cm [mm]"), ("cr2sub_clsoilmap_clay_0_100cm", "Arcilla 0-100 cm [%]"),
                            ("cr2sub_clsoilmap_ksat_0_100cm", "ksat 0-100 cm [cm/día]")]):
    ax.scatter(new.cr2sub_lon[new[c].isna()], new.cr2sub_lat[new[c].isna()], s=12, c="red", marker="x", label=f"sin dato ({int(new[c].isna().sum())})")
    s = ax.scatter(new.cr2sub_lon, new.cr2sub_lat, s=6, c=new[c], cmap="viridis", vmax=new[c].quantile(.98))
    plt.colorbar(s, ax=ax, shrink=.7); ax.set_title(t); ax.set_xlabel("lon"); ax.legend(loc="lower left", fontsize=7)
axs[0].set_ylabel("lat")
fig.suptitle("F. CLSoilMaps en v1.2: distribución espacial y pozos sin dato")
fig.tight_layout(); fig.savefig(OUT / "fig_f_soil_map.png", dpi=110); plt.close(fig)
print("figuras escritas en", OUT)
