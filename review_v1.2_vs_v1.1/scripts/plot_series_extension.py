#!/usr/bin/env python3
"""Plot a time-series chart for one well showing the historical record (old
input) against the points the new source (AM006T0034616_1.xlsx) adds on top
of it. Picks the well with the richest history among those with the most
extension-period coverage, unless WELL_CODE is overridden below.

Writes a PNG into ../comparison_v1.1_vs_v1.2/.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from datetime import datetime

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
V11_DIR = ROOT / "versions" / "v1.1"  # input viejo: reportes dga_xls_*
V13_DIR = ROOT / "versions" / "v1.2"  # input nuevo: AM006T0034616_1.xlsx
EVIDENCE_DIR = ROOT / "comparison_v1.1_vs_v1.2"

WELL_CODE = "5423019"

# Validated 2-slot categorical palette (dataviz skill reference palette,
# fixed order: slot 1 blue, slot 2 orange). Passed scripts/validate_palette.js
# --mode light: lightness band, chroma floor, CVD sep (dE 24.7), normal-vision
# sep (dE 33.6), contrast vs surface all PASS.
COLOR_OLD = "#2a78d6"   # slot 1 - existing record (both sources agree here)
COLOR_NEW = "#eb6834"   # slot 2 - points the new source adds
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#d8d7d2"


def load_series(path, well):
    with path.open() as f:
        r = csv.reader(f)
        header = next(r)
        idx = header.index(well)
        dates, values = [], []
        for row in r:
            v = row[idx]
            dates.append(datetime.strptime(row[0], "%Y-%m-%d"))
            values.append(None if v in ("", "NA") else float(v))
    return dates, values


def main():
    dates_old, vals_old = load_series(V11_DIR / "cr2sub_v1.1_gwl_mon.csv", WELL_CODE)
    dates_new, vals_new = load_series(V13_DIR / "cr2sub_v1.2_gwl_mon.csv", WELL_CODE)

    old_by_date = {d: v for d, v in zip(dates_old, vals_old) if v is not None}
    new_by_date = {d: v for d, v in zip(dates_new, vals_new) if v is not None}
    added_set = {d for d in new_by_date if d not in old_by_date}

    # Build one continuous monthly axis (new output's date range, which is a
    # superset) so real gaps break the line instead of being skipped and
    # drawn as a misleading straight connector.
    all_dates = sorted(set(dates_old) | set(dates_new))
    existing_series = [old_by_date.get(d, np.nan) for d in all_dates]
    added_series = [new_by_date[d] if d in added_set else np.nan for d in all_dates]

    n_existing = sum(1 for v in existing_series if not np.isnan(v))
    n_added = sum(1 for v in added_series if not np.isnan(v))

    # Bridge the join point so the added segment reads as a continuation of
    # the existing line rather than a disconnected island.
    added_bridged = list(added_series)
    last_existing_idx = max((i for i, v in enumerate(existing_series) if not np.isnan(v)), default=None)
    first_added_idx = next((i for i, v in enumerate(added_series) if not np.isnan(v)), None)
    if last_existing_idx is not None and first_added_idx is not None and last_existing_idx < first_added_idx:
        added_bridged[last_existing_idx] = existing_series[last_existing_idx]

    fig, ax = plt.subplots(figsize=(9, 4.2), dpi=200)
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")

    ax.plot(all_dates, existing_series, color=COLOR_OLD, linewidth=2,
             solid_capstyle="round", label=f"Registro existente ({n_existing} meses)",
             zorder=2)
    ax.plot(all_dates, added_bridged, color=COLOR_NEW, linewidth=2,
             solid_capstyle="round", marker="o", markersize=5,
             label=f"Agregado por la actualización ({n_added} meses)",
             zorder=3)

    ax.set_title(f"Pozo DGA {WELL_CODE} — nivel freático mensual",
                 color=TEXT_PRIMARY, fontsize=13, fontweight="bold", loc="left")
    ax.set_ylabel("Nivel (m)", color=TEXT_SECONDARY, fontsize=10)

    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9)
    ax.xaxis.set_major_locator(mdates.YearLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

    legend = ax.legend(loc="lower left", frameon=False, fontsize=9,
                        labelcolor=TEXT_PRIMARY)

    fig.tight_layout()
    out_path = EVIDENCE_DIR / f"series_extension_{WELL_CODE}.png"
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    existing_dates = [d for d, v in zip(all_dates, existing_series) if not np.isnan(v)]
    added_dates = sorted(added_set)
    print(f"Escrito {out_path}")
    print(f"Registro existente: {existing_dates[0].date()} a {existing_dates[-1].date()} "
          f"({n_existing} meses)")
    print(f"Agregado: {added_dates[0].date()} a {added_dates[-1].date()} "
          f"({n_added} meses)")


if __name__ == "__main__":
    main()
