#!/usr/bin/env python3
"""Compare cr2sub outputs built from the old input/DGA_GWL_obs_descarga/dga_xls_*
report folders against outputs built from the single AM006T0034616_1.xlsx export.

Reads the two output snapshots kept in ../versions/{v1.1,v1.3}
and writes evidence files (CSVs + console summary) into ../comparison_v1.1_vs_v1.3/.

Usage:
    python3 compare_outputs.py
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
V11_DIR = ROOT / "versions" / "v1.1"  # input viejo: reportes dga_xls_*
V13_DIR = ROOT / "versions" / "v1.3"  # input nuevo: AM006T0034616_1.xlsx
EVIDENCE_DIR = ROOT / "comparison_v1.1_vs_v1.3"


def load_csv(path):
    with path.open() as f:
        r = csv.reader(f)
        header = next(r)
        rows = list(r)
    return header, rows


def load_indexed(path):
    header, rows = load_csv(path)
    return header, {row[0]: row for row in rows}


def is_na(v):
    return v == "" or v == "NA"


def compare_timeseries(name):
    """Compare a date-indexed, well-column CSV (gwl_mon.csv / gwl_mon_clean.csv)."""
    h_old, d_old = load_indexed(V11_DIR / name)
    h_new, d_new = load_indexed(V13_DIR / name)

    idx_old = {c: i for i, c in enumerate(h_old)}
    idx_new = {c: i for i, c in enumerate(h_new)}

    cols_old, cols_new = set(h_old[1:]), set(h_new[1:])
    common_wells = sorted(cols_old & cols_new)
    common_dates = sorted(set(d_old) & set(d_new))

    gap_filled, lost, conflicts = [], [], []
    for dt in common_dates:
        ro, rn = d_old[dt], d_new[dt]
        for w in common_wells:
            vo, vn = ro[idx_old[w]], rn[idx_new[w]]
            if vo == vn:
                continue
            o_na, n_na = is_na(vo), is_na(vn)
            if o_na and not n_na:
                gap_filled.append((w, dt, vn))
            elif n_na and not o_na:
                lost.append((w, dt, vo))
            else:
                try:
                    if abs(float(vo) - float(vn)) > 1e-3:
                        conflicts.append((w, dt, vo, vn))
                except ValueError:
                    conflicts.append((w, dt, vo, vn))

    result = {
        "name": name,
        "wells_old": len(cols_old),
        "wells_new": len(cols_new),
        "wells_only_old": sorted(cols_old - cols_new),
        "date_range_old": (min(d_old), max(d_old)),
        "date_range_new": (min(d_new), max(d_new)),
        "n_dates_old": len(d_old),
        "n_dates_new": len(d_new),
        "gap_filled": gap_filled,
        "lost": lost,
        "conflicts": conflicts,
    }
    return result


def compare_attributes():
    h_old, d_old = load_indexed(V11_DIR / "cr2sub_v1.1_attributes.csv")
    h_new, d_new = load_indexed(V13_DIR / "cr2sub_v1.1_attributes.csv")

    idx_old = {c: i for i, c in enumerate(h_old)}
    idx_new = {c: i for i, c in enumerate(h_new)}
    common_ids = sorted(set(d_old) & set(d_new))
    missing_ids = sorted(set(d_old) - set(d_new))

    metadata_cols = [
        "dga_well_code", "dga_well_name", "dga_well_basin", "dga_well_subbasin",
        "cr2sub_lat", "cr2sub_lon", "dga_well_elev",
    ]
    diffs = []
    for cid in common_ids:
        for col in ["cr2sub_lat", "cr2sub_lon", "dga_well_elev"]:
            vo, vn = d_old[cid][idx_old[col]], d_new[cid][idx_new[col]]
            try:
                fo, fn = float(vo), float(vn)
            except ValueError:
                continue
            if abs(fo - fn) > 0.01:
                diffs.append((cid, d_old[cid][idx_old["dga_well_name"]], col, vo, vn))

    missing_rows = [d_old[cid] for cid in missing_ids]
    return {
        "header_old": h_old,
        "wells_old": len(d_old),
        "wells_new": len(d_new),
        "missing_ids": missing_ids,
        "missing_rows": missing_rows,
        "idx_old": idx_old,
        "metadata_diffs": diffs,
        "metadata_cols": metadata_cols,
    }


def pick_extension_examples(ts_result, n=3):
    """Pick a few wells that show data in the newly-extended tail (post Dec-2024)."""
    h_new, d_new = load_indexed(V13_DIR / ts_result["name"])
    idx_new = {c: i for i, c in enumerate(h_new)}
    extension_dates = sorted(dt for dt in d_new if dt > "2024-12")
    examples = []
    for w in h_new[1:]:
        found = []
        for dt in extension_dates:
            v = d_new[dt][idx_new[w]]
            if not is_na(v):
                found.append((dt, v))
        if found:
            examples.append((w, found))
        if len(examples) >= n:
            break
    return extension_dates, examples


def main():
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

    mon = compare_timeseries("cr2sub_v1.1_gwl_mon.csv")
    mon_clean = compare_timeseries("cr2sub_v1.1_gwl_mon_clean.csv")
    attrs = compare_attributes()
    ext_dates, ext_examples = pick_extension_examples(mon)

    # --- write evidence CSVs ---
    with (EVIDENCE_DIR / "metadata_diffs.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cr2sub_id", "dga_well_name", "column", "old_value", "new_value"])
        w.writerows(attrs["metadata_diffs"])

    with (EVIDENCE_DIR / "wells_missing_in_new_source.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(attrs["header_old"])
        w.writerows(attrs["missing_rows"])

    for label, result in [("gwl_mon", mon), ("gwl_mon_clean", mon_clean)]:
        with (EVIDENCE_DIR / f"{label}_gap_filled.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["dga_well_code", "fecha", "new_value"])
            w.writerows(result["gap_filled"])
        with (EVIDENCE_DIR / f"{label}_lost.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["dga_well_code", "fecha", "old_value"])
            w.writerows(result["lost"])
        with (EVIDENCE_DIR / f"{label}_conflicts.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["dga_well_code", "fecha", "old_value", "new_value"])
            w.writerows(result["conflicts"])

    with (EVIDENCE_DIR / "series_extension_examples.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dga_well_code", "fecha", "nivel"])
        for well, found in ext_examples:
            for dt, v in found:
                w.writerow([well, dt, v])

    # --- console summary ---
    for result in (mon, mon_clean):
        print(f"== {result['name']} ==")
        print(f"  wells: old={result['wells_old']} new={result['wells_new']} "
              f"only_in_old={result['wells_only_old']}")
        print(f"  date range old: {result['date_range_old']}  (n={result['n_dates_old']})")
        print(f"  date range new: {result['date_range_new']}  (n={result['n_dates_new']})")
        print(f"  gap_filled={len(result['gap_filled'])} "
              f"lost={len(result['lost'])} conflicts={len(result['conflicts'])}")
        print()

    print("== attributes ==")
    print(f"  wells: old={attrs['wells_old']} new={attrs['wells_new']} "
          f"missing_in_new={attrs['missing_ids']}")
    print(f"  metadata diffs (>0.01): {len(attrs['metadata_diffs'])}")
    print()
    print(f"Evidence written to {EVIDENCE_DIR}")


if __name__ == "__main__":
    main()
