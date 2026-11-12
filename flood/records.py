"""Checks of the daily rainfall series and the flood-event table, and the relation between the event columns."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

SEVERITY_ORDER = ["Low", "Moderate", "High"]
MONSOON_MONTHS = [6, 7, 8, 9]


def load_rainfall(path):
    """Daily rainfall in mm, one column per district, indexed by date and without gaps."""
    rain = pd.read_csv(path)
    rain["Date"] = pd.to_datetime(rain["Date"])
    rain = rain.set_index("Date").sort_index().astype(float)
    expected = pd.date_range(rain.index[0], rain.index[-1], freq="D")
    if len(rain) != len(expected) or not rain.index.equals(expected):
        raise ValueError("the rainfall series must have exactly one row per day")
    return rain


def load_events(path):
    """Flood events with the columns region, date, severity, area_km2 and rainfall_mm."""
    events = pd.read_csv(path).rename(columns={
        "Region": "region", "Date": "date", "Flood Severity": "severity", "Affected Area (sq km)": "area_km2",
        "Rainfall (mm)": "rainfall_mm"})
    missing = {"region", "date", "severity", "area_km2", "rainfall_mm"} - set(events.columns)
    if missing:
        raise ValueError(f"missing columns: {sorted(missing)}")
    unknown = sorted(set(events["severity"]) - set(SEVERITY_ORDER))
    if unknown:
        raise ValueError(f"unknown severity levels: {unknown}")
    events["date"] = pd.to_datetime(events["date"])
    return events


def monthly_summary(rain, wet_day_mm=1.0):
    """Monthly total, wet days and largest daily value per district, for a series of one year."""
    months = rain.groupby(rain.index.month)
    return {"total_mm": months.sum(), "wet_days": months.agg(lambda values: int((values >= wet_day_mm).sum())),
            "largest_mm": months.max()}


def rainfall_checks(rain):
    """Facts that show whether the series can be gauge records."""
    monthly = rain.groupby(rain.index.month).sum()
    wet = monthly[(monthly > 0).all(axis=1)]
    # monthly totals of each district as multiples of its smallest wet month
    multiple = wet / wet.min(axis=0)
    monsoon = rain[rain.index.month.isin(MONSOON_MONTHS)]
    rank = monsoon.corr(method="spearman").to_numpy()
    pairs = rank[np.triu_indices_from(rank, k=1)]
    both = (monsoon > 0).all(axis=1).mean()
    return {
        "annual_total_mm": rain.sum().round(1).to_dict(),
        "dry_months": [int(month) for month in monthly.index[(monthly == 0).all(axis=1)]],
        "monthly_multiples": {name: [float(value) for value in multiple[name].round(3)] for name in multiple},
        "largest_departure_from_whole_multiples": float((multiple - multiple.round()).abs().to_numpy().max()),
        "largest_difference_between_districts": float((multiple.max(axis=1) - multiple.min(axis=1)).max()),
        "monsoon_rank_correlation_between_districts": [float(pairs.min()), float(pairs.max())],
        "monsoon_days_wet_in_all_districts": float(both),
        "monsoon_days_wet_in_all_districts_if_independent": float((monsoon > 0).mean().prod()),
    }


def event_rainfall(events, rain, windows=(1, 3, 7)):
    """Rainfall of the series in each event's district over windows of days that end on the event date."""
    unknown = sorted(set(events["region"]) - set(rain.columns))
    if unknown:
        raise ValueError(f"no rainfall series for: {unknown}")
    table = events[["region", "date", "rainfall_mm"]].copy()
    for days in windows:
        table[f"series_{days}d_mm"] = [
            float(rain.loc[date - pd.Timedelta(days=days - 1):date, region].sum())
            for date, region in zip(events["date"], events["region"])]
    return table


def severity_ranges(events):
    """Smallest and largest area and rainfall of each severity class, and whether neighbouring classes overlap."""
    table = events.groupby("severity")[["area_km2", "rainfall_mm"]].agg(["min", "max"]).reindex(SEVERITY_ORDER)
    overlap = {column: bool((table[column]["max"].to_numpy()[:-1] >= table[column]["min"].to_numpy()[1:]).any())
               for column in ("area_km2", "rainfall_mm")}
    return table, overlap


def area_rainfall_relation(events):
    """Straight-line fit of affected area on rainfall with its leave-one-out prediction error."""
    x = events["rainfall_mm"].to_numpy(float)
    y = events["area_km2"].to_numpy(float)
    n = len(x)
    fit = stats.linregress(x, y)
    pearson = stats.pearsonr(x, y)
    low, high = pearson.confidence_interval(0.95)
    predicted, mean_only = np.empty(n), np.empty(n)
    for i in range(n):
        others = np.arange(n) != i
        part = stats.linregress(x[others], y[others])
        predicted[i] = part.intercept + part.slope * x[i]
        mean_only[i] = y[others].mean()
    return {
        "events": int(n), "slope_km2_per_mm": float(fit.slope), "intercept_km2": float(fit.intercept),
        "rainfall_at_zero_area_mm": float(-fit.intercept / fit.slope),
        "r": float(pearson.statistic), "r_low": float(low), "r_high": float(high), "p": float(pearson.pvalue),
        "spearman": float(stats.spearmanr(x, y).statistic),
        "rmse_fit_km2": float(np.sqrt(np.mean((fit.intercept + fit.slope * x - y) ** 2))),
        "rmse_leave_one_out_km2": float(np.sqrt(np.mean((predicted - y) ** 2))),
        "rmse_leave_one_out_mean_only_km2": float(np.sqrt(np.mean((mean_only - y) ** 2))),
        "leave_one_out_prediction_km2": predicted.tolist(),
    }
