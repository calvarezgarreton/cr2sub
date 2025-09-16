#!/usr/bin/env python
# coding: utf-8

# In[1]:

# functions_outliers.py
# -----------------------------------------------------------------------------
# We remove outliers with a three-stage, NaN-aware procedure tailored to groundwater 
# dynamics: (1) suppress “V/Λ” spikes on consecutive months using a P99 threshold on 
# consecutive jumps; (2) iteratively remove physically implausible drops/rises relative 
# to a local moving mean, gated by large residuals and (for rises) requiring future support; 
# and (3) perform a conservative distributional cleanup using skew-aware one-sided trimming 
# and two-sided robust-z tails, both guarded by residual gates computed from exclusion or 
# inlier-interpolated means. Thresholds are estimated from each well’s empirical distributions, 
# and all decisions avoid imputing gaps except when explicitly interpolating inliers to stabilize 
# the final reference.

# C. Alvarez-Garreton. 12 August 2025.
# -----------------------------------------------------------------------------

import numpy as np
import pandas as pd

__all__ = ["mean_at_dates", "run_pass", "final_tail_filters", "remove_isolated_islands"]


def mean_at_dates(series: pd.Series, dates: pd.Index, win_months: int) -> pd.Series:
    """
    NaN-aware symmetric moving mean centered at given dates.

    Builds a complete monthly grid from the observed span, grabs a ±win_months
    window around each date, and returns the mean ignoring NaNs.

    Parameters
    ----------
    series : pd.Series
        Groundwater series with a DateTimeIndex (monthly is ideal). May contain NaNs.
    dates : pd.Index
        Dates at which to compute the local mean. Must be within series span.
    win_months : int
        Half window in months for the symmetric mean.

    Returns
    -------
    pd.Series
        Local mean at `dates`. If a window is all-NaN, the result is NaN.
    """
    dates = pd.DatetimeIndex(dates)
    if len(dates) == 0 or series.dropna().empty:
        return pd.Series([], dtype=float)

    obs_ = series.dropna()
    full_idx = pd.date_range(obs_.index.min(), obs_.index.max(), freq="MS")
    s_full = series.reindex(full_idx)  # keep NaNs

    out = []
    for d in dates:
        left = d - pd.DateOffset(months=win_months)
        right = d + pd.DateOffset(months=win_months)
        w = s_full.loc[(s_full.index >= left) & (s_full.index <= right)].values
        out.append(np.nan if np.isnan(w).all() else np.nanmean(w))
    return pd.Series(out, index=dates)


def run_pass(s_in: pd.Series, win_months: int, perc_thr: float):
    """
    One filtering pass: Rule 0 (spike/dip) + Rule 1 (phys-aware jumps).

    Rule 0:
      Remove “spike/dip middles” where (i-1→i) and (i→i+1) are both large AND
      opposite-signed (V/Λ). Threshold = P99 of |consecutive raw diffs|.

    Rule 1:
      On the series after Rule 0, compute local means and residuals. For points
      whose |residual| > P(perc_thr) of |residuals|, compare the observation
      to the mean of past kept observations within ±win_months:
        - If a large DROP (< -P(perc_thr) of |Δ mean|), remove.
        - If a large RISE (>  P(perc_thr) of |Δ mean|), also require “future
          support”: the next ±win_months kept-mean must be close; otherwise remove.

    Parameters
    ----------
    s_in : pd.Series
        Input groundwater series (DateTimeIndex, monthly recommended).
    win_months : int
        Half-window in months for local means.
    perc_thr : float
        Percentile in (0,1), e.g. 0.99.

    Returns
    -------
    s_out : pd.Series
        Series after Rule 0 and Rule 1 (NaNs where removed).
    rule0_dates : pd.Index
        Dates removed by Rule 0.
    rule1_dates : pd.Index
        Dates removed by Rule 1.
    smooth_series : pd.Series
        Local mean at post-Rule0 observed dates (for plotting).
    """
    obs = s_in.dropna()
    if len(obs) == 0:
        return s_in.copy(), pd.Index([]), pd.Index([]), pd.Series([], dtype=float)

    # Month counter and values for “consecutive” logic (avoid gaps)
    m = (obs.index.year * 12 + obs.index.month).to_numpy()
    v = obs.to_numpy(float)
    n = len(obs)
    is_consec = (np.diff(m) == 1)

    # -------- Rule 0: remove spike/dip middles on consecutive months ----------
    dv = np.diff(v)  # signed diffs on observed-consecutive pairs
    if is_consec.any():
        thr_consec = np.quantile(np.abs(dv[is_consec]), 0.99)
        prev_diff = np.full(n, np.nan)
        next_diff = np.full(n, np.nan)
        idx_consec = np.where(is_consec)[0]
        prev_diff[idx_consec + 1] = dv[idx_consec]  # (i-1→i)
        next_diff[idx_consec]     = dv[idx_consec]  # (i→i+1)
        rule0_idx = np.where(
            ((prev_diff <= -thr_consec) & (next_diff >=  thr_consec)) |
            ((prev_diff >=  thr_consec) & (next_diff <= -thr_consec))
        )[0]
    else:
        rule0_idx = np.array([], dtype=int)

    s1 = s_in.copy()
    rule0_dates = obs.index[rule0_idx]
    s1.loc[rule0_dates] = np.nan  # drop only the middle spike/dip point

    # -------- Rule 1: physically-aware jumps, residual-gated ------------------
    obs2 = s1.dropna()
    if len(obs2) == 0:
        return s1, rule0_dates, pd.Index([]), pd.Series([], dtype=float)

    m2 = (obs2.index.year * 12 + obs2.index.month).to_numpy()
    v2 = obs2.to_numpy(float)
    n2 = len(obs2)

    means = mean_at_dates(s1, obs2.index, win_months).to_numpy(float)
    residuals = v2 - means
    thr_res = np.quantile(np.abs(residuals), perc_thr) if n else np.inf
    
    is_consec = (np.diff(m2) == 1)
    d_mu = np.abs(np.diff(means)[is_consec]) if is_consec.any() else np.array([])
    thr_jump = np.quantile(d_mu, perc_thr) if d_mu.size else np.inf
    
    keep = np.ones(n2, dtype=bool)

    for i in range(n2):
        if np.abs(residuals[i]) <= thr_res:
            continue  # only consider points with large residuals

        # Past kept observations within the window
        L = np.searchsorted(m2, m2[i] - win_months, side="left")
        pv = v2[L:i][keep[L:i]]
        if pv.size == 0:
            continue

        diff_prev = v2[i] - pv.mean()

        if diff_prev < -thr_jump:  # large DROP vs past
            keep[i] = False

        elif diff_prev > thr_jump:  # large RISE vs past → require future support
            R = np.searchsorted(m2, m2[i] + win_months, side="right")
            fv = v2[i+1:R][keep[i+1:R]]
            # If no future kept obs, or rise sits far above future mean, drop
            if fv.size == 0 or (v2[i] - fv.mean()) > thr_jump:
                keep[i] = False

    rule1_dates = obs2.index[~keep]
    s_out = s1.copy()
    s_out.loc[rule1_dates] = np.nan

    smooth_series = pd.Series(means, index=obs2.index)
    return s_out, rule0_dates, rule1_dates, smooth_series


def final_tail_filters(
    s_in: pd.Series,
    win_months: int,
    perc_thr: float,
    tail_p: float = 0.01,
    skew_thr: float = 0.30,
    min_n: int = 50,
    skew_tail_p: float = 0.01,
):
    """
    Final distributional cleanup with two conservative, gated steps.

    (A) One-sided (Bowley) skew tail:
        - Trigger only if |skew_Q| > skew_thr.
        - Select high (or low) tail via P(1 - skew_tail_p) or P(skew_tail_p).
        - Exclude these candidates, compute residual threshold from the remaining
          data’s moving mean, and remove candidates whose residual is still large.

    (B) Two-sided extreme tails:
        - Start with P(tail_p)/P(1-tail_p) cut and robust z (>4.5) around the median.
        - Build an inlier-only series (Plo–Phi), interpolate to stabilize the mean,
          and remove only candidates whose residual vs this inlier-mean exceeds
          P(perc_thr) of inlier residuals.

    Returns
    -------
    s_out : pd.Series
        Series after applying (A) and (B).
    tail_idx1 : pd.Index
        Dates removed by skew-aware step (A).
    tail_idx2 : pd.Index
        Dates removed by two-sided step (B).
    """
    gw_filt2 = s_in.copy()

    # ----------------------- (A) Skew-aware one-sided tail ----------------------
    vals = gw_filt2.dropna().to_numpy(float)
    tail_idx1 = pd.Index([])

    if vals.size >= min_n:
        q1, med, q3 = np.quantile(vals, [0.25, 0.50, 0.75])
        iqr = q3 - q1
        if iqr > 0:
            skew_q = (q3 + q1 - 2 * med) / iqr  # Bowley skewness
            if skew_q > skew_thr:
                hi = np.quantile(vals, 1 - skew_tail_p)
                cand = gw_filt2.index[gw_filt2 > hi]
            elif skew_q < -skew_thr:
                lo = np.quantile(vals, skew_tail_p)
                cand = gw_filt2.index[gw_filt2 < lo]
            else:
                cand = pd.Index([])

            if len(cand):
                # Exclude candidates, compute residual threshold from non-candidates
                s_excl = gw_filt2.copy()
                s_excl.loc[cand] = np.nan
                obs_ex = s_excl.dropna()

                if len(obs_ex):
                    mu_ex_obs = mean_at_dates(s_excl, obs_ex.index, win_months)
                    thr_res_ex = np.quantile(
                        np.abs((obs_ex - mu_ex_obs).to_numpy(float)),
                        perc_thr,
                    )
                else:
                    thr_res_ex = np.inf

                # Mean at candidate dates WITHOUT candidates
                mu_ex_cand = mean_at_dates(s_excl, cand, win_months)
                resid_cand = (gw_filt2.loc[cand] - mu_ex_cand).abs()
                tail_idx1 = resid_cand.index[resid_cand > thr_res_ex]

    gw_final0 = gw_filt2.copy()
    if len(tail_idx1):
        gw_final0.loc[tail_idx1] = np.nan

    # -------- (B) Two-sided Plo/Phi + robust-z + inlier-mean residual gate -------
    vals2 = gw_final0.dropna().to_numpy(float)
    tail_idx2 = pd.Index([])

    if vals2.size:
        p_lo_q, med2, p_hi_q = np.quantile(vals2, [tail_p, 0.50, 1 - tail_p])
        cand0 = gw_final0.index[(gw_final0 < p_lo_q) | (gw_final0 > p_hi_q)]

        if len(cand0):
            mad = np.median(np.abs(vals2 - med2))
            scale = 1.4826 * mad if mad > 0 else (np.quantile(vals2, 0.75) - np.quantile(vals2, 0.25)) / 1.349
            if not np.isfinite(scale) or scale <= 0:
                s = np.nanstd(vals2)
                scale = s if s > 0 else 1.0

            z = (gw_final0.loc[cand0] - med2).abs() / scale
            cand = z.index[z > 4]  # conservative

            if len(cand):
                # Inlier-only series (between Plo and Phi), then interpolate
                mask_inliers = (gw_final0 >= p_lo_q) & (gw_final0 <= p_hi_q)
                s_inliers = gw_final0.where(mask_inliers)
                s_interp = s_inliers.interpolate(method="time", limit_direction="both")

                # Local mean for any date on the interpolated inliers
                def mean_at_dates_interp(series: pd.Series, dates: pd.Index, win_m: int) -> pd.Series:
                    dates = pd.DatetimeIndex(dates)
                    if len(dates) == 0:
                        return pd.Series([], dtype=float)
                    full_idx = pd.date_range(series.index.min(), series.index.max(), freq="MS")
                    s_full = series.reindex(full_idx)
                    out = []
                    for d in dates:
                        left = d - pd.DateOffset(months=win_m)
                        right = d + pd.DateOffset(months=win_m)
                        w = s_full.loc[(s_full.index >= left) & (s_full.index <= right)].values
                        out.append(np.nan if np.isnan(w).all() else np.nanmean(w))
                    return pd.Series(out, index=dates)

                mu_cand = mean_at_dates_interp(s_interp, cand, win_months)
                mu_all = mean_at_dates_interp(s_interp, s_interp.dropna().index, win_months)
                thr_resI = np.quantile((s_interp.dropna() - mu_all).abs().to_numpy(float), perc_thr)

                resid_c = (gw_final0.loc[cand] - mu_cand).abs()
                tail_idx2 = resid_c.index[resid_c > thr_resI]

    s_out = gw_final0.copy()
    if len(tail_idx2):
        s_out.loc[tail_idx2] = np.nan

    return s_out, tail_idx1, tail_idx2


def remove_isolated_islands(s_in, min_nan_gap=6, max_island_len=2):
    """
    Drop short runs of observations (length <= max_island_len, in months)
    that are isolated by long NaN runs (>= min_nan_gap) on both sides,
    or at the edges (only the available side).

    Returns: (filtered_series, removed_dates_index)
    """
    s = s_in.copy()
    if s.dropna().empty:
        return s, pd.Index([])

    # Monthly grid so lengths == months
    start = s.index.min().replace(day=1)
    end   = s.index.max().replace(day=1)
    idx   = pd.date_range(start, end, freq="MS")
    sf    = s.reindex(idx)

    is_obs = sf.notna().to_numpy()

    # Find starts/ends of TRUE runs (observations)
    b = np.r_[0, is_obs.astype(int), 0]   # cast to int so diff is -1/0/1
    d = np.diff(b)
    starts = np.where(d == 1)[0]          # indices into sf/index
    ends   = np.where(d == -1)[0] - 1

    if len(starts) == 0:
        return s, pd.Index([])

    N = len(sf)
    drop_pos = []

    for i in range(len(starts)):
        L_obs = ends[i] - starts[i] + 1
        if L_obs > max_island_len:
            continue

        # NaN-gap to the left/right (in months)
        left_gap  = (starts[i] - (ends[i-1] + 1)) if i > 0 else starts[i]
        right_gap = (starts[i+1] - (ends[i] + 1)) if i < len(starts)-1 else (N - 1 - ends[i])

        # Interior: need both gaps long; Edges: only the existing side
        if ((0 < i < len(starts)-1 and left_gap >= min_nan_gap and right_gap >= min_nan_gap) or
            (i == 0 and right_gap >= min_nan_gap) or
            (i == len(starts)-1 and left_gap >= min_nan_gap)):
            drop_pos.extend(range(starts[i], ends[i] + 1))

    if not drop_pos:
        return s, pd.Index([])

    removed = sf.index[drop_pos].intersection(s.index)
    s.loc[removed] = np.nan
    return s, removed





