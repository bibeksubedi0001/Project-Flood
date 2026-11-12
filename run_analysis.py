"""Low ground along the rivers of one DEM scene, and checks of the rainfall and flood-event tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, LightSource, ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.patheffects import withStroke
from rasterio.warp import transform as reproject
from scipy import ndimage

from flood import (MONSOON_MONTHS, SEVERITY_ORDER, area_rainfall_relation, block_mean, event_rainfall,
                   height_above_drainage, load_events, load_rainfall, monthly_summary, rainfall_checks, read_dem,
                   resample, resampling_test, route_flow, severity_ranges, slope, stream_level, value_at_drainage,
                   write_raster)

# Liberation Serif has the metrics of Times New Roman and stands in where that font is not installed
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"]})

ROOT = Path(__file__).resolve().parent
FIGURES, RESULTS = ROOT / "figures", ROOT / "results"
SCENE = "AP_27418_PLR_F0540_RT1"
DEM_FILE = ROOT / SCENE / SCENE / f"{SCENE}.dem.tif"
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

THRESHOLDS = (5, 25, 100)                  # a cell is a stream when this area in km² drains to it inside the scene
REFERENCE = 25                             # the threshold of the maps
HEIGHTS = (5, 10, 20, 50)                  # m above the stream
LEVELS = ("filled", "fitted", "carved")    # three readings of the level of a stream, see flood.stream_level
CELL_FACTORS = (1, 2, 4)                   # cell sizes, in cells of the file, on which the main result is repeated
# Stretches of the edge of the scene, in km of easting or northing, across the valley floors through which the
# two rivers enter. No flow leaves there, and the course of each river is a stream from the edge on.
ENTRIES = {"Dudh Koshi": ("north", 470.5, 472.5), "Sun Koshi": ("west", 2999.8, 3001.9)}
TOWNS = {"Okhaldhunga": (86.50, 27.32), "Diktel": (86.80, 27.21)}       # approximate, for orientation only
WINDOW_KM = (12.0, 9.0)                    # width and height of the enlarged map
PIT_DEPTH = 10                             # m below the level of its stream from which ground counts as lying in a pit
RESAMPLING_MARGIN = 100                    # cells left out at the edge of the file, where data are missing

HEIGHT_COLOURS = ("#08306b", "#2171b5", "#6baed6", "#c6dbef")
LEVEL_COLOURS = {"filled": "tab:blue", "fitted": "black", "carved": "tab:red"}
SEVERITY_COLOURS = {"Low": "tab:green", "Moderate": "tab:orange", "High": "tab:red"}
OUTLINE = [withStroke(linewidth=2.5, foreground="white")]


def finish(fig, name):
    fig.savefig(FIGURES / name, dpi=150, bbox_inches="tight")
    plt.close(fig)


def place(dem, index):
    """Position of a cell given by its flat index, for the result tables."""
    row, col = divmod(int(index), dem.z.shape[1])
    east, north = dem.position(row, col)
    lon, lat = reproject(dem.crs, "EPSG:4326", [east], [north])
    return {"easting_km": east / 1000, "northing_km": north / 1000, "longitude": float(lon[0]), "latitude": float(lat[0])}


def spread(values):
    return {"mean": float(values.mean()), "median": float(np.median(values)),
            "percentile_5": float(np.percentile(values, 5)), "percentile_95": float(np.percentile(values, 95)),
            "rms": float(np.sqrt(np.mean(values**2))), "largest": float(values[np.argmax(np.abs(values))])}


def number(value, digits=1):
    """A number for the result tables; a value that rounds to zero is written without a sign."""
    return f"{round(value, digits) + 0.0:.{digits}f}"


def signed(text):
    """The same text with a true minus sign, for the figures."""
    return text.replace("-", "\N{MINUS SIGN}")


def extent_km(dem, rows=None, cols=None):
    """Left, right, bottom and top of a block of cells in km, in the order imshow expects."""
    rows = rows or (0, dem.z.shape[0])
    cols = cols or (0, dem.z.shape[1])
    left, top = dem.transform * (cols[0], rows[0])
    right, bottom = dem.transform * (cols[1], rows[1])
    return [left / 1000, right / 1000, bottom / 1000, top / 1000]


def hillshade(z, cell):
    return LightSource(azdeg=315, altdeg=45).hillshade(np.where(np.isnan(z), np.nanmin(z), z), dx=cell, dy=cell)


def block_min(values, factor):
    """Smallest value in each block of factor x factor cells, ignoring missing cells."""
    rows, cols = (np.array(values.shape) // factor) * factor
    blocks = np.where(np.isnan(values), np.inf, values)[:rows, :cols].reshape(rows // factor, factor, cols // factor, factor)
    smallest = blocks.min(axis=(1, 3))
    return np.where(np.isinf(smallest), np.nan, smallest)


def height_classes(height):
    """Index of the height class of each cell, masked above the last class and where there is no value."""
    classes = np.digitize(np.nan_to_num(height, nan=np.inf), HEIGHTS, right=True)
    return np.ma.masked_equal(classes, len(HEIGHTS))


def class_legend(extra=()):
    labels = [f"up to {HEIGHTS[0]} m"] + [f"{low} to {high} m" for low, high in zip(HEIGHTS[:-1], HEIGHTS[1:])]
    return [Patch(color=colour, label=label) for colour, label in zip(HEIGHT_COLOURS, labels)] + list(extra)


def river_entries(dem):
    """The edge cells in the stretches of ENTRIES and, in each stretch, the lowest cell: there the river enters."""
    valid = np.isfinite(dem.z)
    edge = valid & ~ndimage.binary_erosion(valid, structure=np.ones((3, 3)))
    rows, cols = np.indices(dem.z.shape)
    east, north = dem.transform * (cols + 0.5, rows + 0.5)
    left, bottom, right, top = dem.bounds()
    near = {"west": east < left + 1000, "east": east > right - 1000, "north": north > top - 1000, "south": north < bottom + 1000}
    closed, entry = np.zeros(dem.z.shape, bool), {}
    for name, (side, start, end) in ENTRIES.items():
        along = north if side in ("west", "east") else east
        stretch = edge & near[side] & (along >= 1000 * start) & (along <= 1000 * end)
        closed |= stretch
        cells = np.flatnonzero(stretch.ravel())
        entry[name] = int(cells[np.argmin(dem.z.ravel()[cells])])
    return closed, entry


def prepare(dem):
    """Flow routing with the river entries closed, the contributing area in km² and the cells of the two rivers."""
    closed, entry = river_entries(dem)
    routing = route_flow(dem.z, closed)
    area = routing.cells * dem.cell_area_km2
    river = np.zeros(dem.z.size, bool)
    for cell in entry.values():
        river[routing.flow_path(cell)] = True
    return routing, area, river.reshape(dem.z.shape), entry


def describe_scene(native):
    valid = np.isfinite(native.z)
    values = native.z[valid]
    west, south, east, north = native.bounds_lonlat()
    left, bottom, right, top = native.bounds()
    return {
        "file": DEM_FILE.name, "epsg": native.crs.to_epsg(), "cell_m": native.cell,
        "rows": int(native.z.shape[0]), "columns": int(native.z.shape[1]),
        "width_km": (right - left) / 1000, "height_km": (top - bottom) / 1000,
        "area_km2": float(valid.sum() * native.cell_area_km2), "share_without_data": float(1 - valid.mean()),
        "west": west, "south": south, "east": east, "north": north,
        "lowest_m": float(values.min()), "highest_m": float(values.max()), "median_m": float(np.median(values)),
        "whole_metres": bool(np.all(values == np.round(values))),
    }


def analyse_resampling(native):
    margin = RESAMPLING_MARGIN
    window = native.z[margin:-margin, margin:-margin].astype(np.float32)
    return {"cells": int(window.size), "from_cell_m": 2 * native.cell, "to_cell_m": native.cell,
            "mean_difference_between_neighbours_m": float(np.abs(np.diff(window, axis=1)).mean(dtype=np.float64)),
            "methods": resampling_test(window, native.cell)}


def analyse_terrain(dem, routing, area):
    valid = np.isfinite(dem.z)
    degrees = np.degrees(np.arctan(slope(dem.z, dem.cell)))
    depth = routing.filled - dem.z
    _, depressions = ndimage.label(np.nan_to_num(depth) > 1.0, structure=np.ones((3, 3)))
    deepest = int(np.nanargmax(depth))
    outlets = routing.outlets()
    return {
        "cell_m": dem.cell, "cells": int(valid.sum()), "area_km2": float(valid.sum() * dem.cell_area_km2),
        "slope_mean_deg": float(np.nanmean(degrees)), "slope_median_deg": float(np.nanmedian(degrees)),
        "share_flatter_than_5_deg": float(np.nanmean(degrees < 5)),
        "share_steeper_than_30_deg": float(np.nanmean(degrees > 30)),
        "share_raised_more_than_1_m_by_filling": float((depth[valid] > 1.0).mean()),
        "area_raised_more_than_1_m_km2": float((depth[valid] > 1.0).sum() * dem.cell_area_km2),
        "depressions_deeper_than_1_m": int(depressions),
        "largest_fill_depth_m": float(depth.ravel()[deepest]), "largest_fill_at": place(dem, deepest),
        "outlets": int(len(outlets)),
        "largest_outlets": [dict(place(dem, index), area_km2=float(area.ravel()[index]),
                                 elevation_m=float(dem.z.ravel()[index])) for index in outlets[:5]],
    }


def analyse_low_ground(dem, routing, area, river, thresholds=THRESHOLDS):
    """Streams by each threshold, their level by each reading, and the ground within given heights of them."""
    valid = np.isfinite(dem.z)
    steps = np.arange(0.0, HEIGHTS[-1] + 0.5)
    table, curves, kept = {}, {}, {}
    for threshold in thresholds:
        stream = (np.nan_to_num(area) >= threshold) | river
        length = routing.step_length(np.flatnonzero(stream.ravel())).sum() * dem.cell / 1000
        entry = {"stream_length_km": float(length), "stream_cells_km2": float(stream.sum() * dem.cell_area_km2)}
        for level in LEVELS:
            levels = stream_level(routing, stream, level)
            height = height_above_drainage(routing, stream, level)
            ordered = np.sort(height[height <= HEIGHTS[-1]])
            curves[(threshold, level)] = np.searchsorted(ordered, steps, side="right") * dem.cell_area_km2
            entry[level] = {"area_km2": {f"{h} m": float((height <= h).sum() * dem.cell_area_km2) for h in HEIGHTS},
                            "level_minus_dem_m": spread(levels[stream] - dem.z[stream])}
            if threshold == REFERENCE:
                kept[level] = {"levels": levels, "height": height}
        entry["share_reaching_no_stream"] = float(np.isnan(height[valid]).mean())
        table[f"{threshold} km²"] = entry
    return table, steps, curves, kept


def analyse_rivers(dem, routing, area, entry, kept):
    """The course of each river from where it enters the scene to where it leaves it."""
    rivers = []
    for name, cell in entry.items():
        path = routing.flow_path(cell)
        step = routing.step_length(path) * dem.cell / 1000
        distance = np.concatenate([[0.0], np.cumsum(step[:-1])])        # km downstream of the entry
        z = dem.z.ravel()[path]
        inside = area.ravel()[path]
        difference = {level: kept[level]["levels"].ravel()[path] - z for level in LEVELS}
        deepest = int(np.argmax(difference["filled"]))
        rivers.append({"name": name, "path": path, "distance": distance, "z": z, "difference": difference, "facts": {
            "name": name, "side": ENTRIES[name][0],
            "enters": place(dem, path[0]), "enters_elevation_m": float(z[0]),
            "leaves": place(dem, path[-1]), "leaves_elevation_m": float(z[-1]),
            "length_km": float(distance[-1]), "mean_fall_m_per_km": float((z[0] - z[-1]) / distance[-1]),
            "area_inside_scene_km2": float(inside[-1]),
            "km_from_entry_until_area": {f"{threshold} km²": float(distance[np.argmax(inside >= threshold)])
                                         for threshold in THRESHOLDS if (inside >= threshold).any()},
            "level_minus_dem_m": {level: spread(difference[level]) for level in LEVELS},
            "filled_minus_carved_median_m": float(np.median(difference["filled"] - difference["carved"])),
            "largest_fill_m": float(difference["filled"][deepest]), "largest_fill_km_from_entry": float(distance[deepest]),
        }})
    return rivers


def analyse_pits(dem, height):
    """Ground far below the level of its stream: pits of the elevation model, which the map shows as low ground."""
    below = np.nan_to_num(height, nan=0.0) < -PIT_DEPTH
    labels, count = ndimage.label(below, structure=np.ones((3, 3)))
    sizes = ndimage.sum(below, labels, np.arange(1, count + 1)) * dem.cell_area_km2
    largest = int(np.argmax(sizes)) + 1
    east, north = dem.position(*ndimage.center_of_mass(labels == largest))
    return {"depth_m": PIT_DEPTH, "area_km2": float(below.sum() * dem.cell_area_km2), "patches": int(count),
            "largest_patch_km2": float(sizes[largest - 1]), "largest_patch_easting_km": east / 1000,
            "largest_patch_northing_km": north / 1000,
            "area_below_the_stream_level_km2": float((height < 0).sum() * dem.cell_area_km2)}


def analyse_by_river(dem, routing, area, river, rivers, height):
    """The low ground of the reference threshold, split by the stream that the ground drains to."""
    stream = (np.nan_to_num(area) >= REFERENCE) | river
    code = np.zeros(dem.z.size)
    for position, one in enumerate(rivers, 1):
        code[one["path"]] = position
    owner = value_at_drainage(routing, stream, code)
    total = routing.step_length(np.flatnonzero(stream.ravel())).sum() * dem.cell / 1000
    lengths = [float(one["distance"][-1]) for one in rivers]
    parts = [(one["name"], position, length) for position, (one, length) in enumerate(zip(rivers, lengths), 1)]
    parts.append(("other streams", 0, float(total - sum(lengths))))
    table = {}
    for name, position, length in parts:
        areas = {f"{h} m": float(((height <= h) & (owner == position)).sum() * dem.cell_area_km2) for h in HEIGHTS}
        table[name] = {"length_km": length, "area_km2": areas,
                       "mean_width_m": {key: 1000 * value / length for key, value in areas.items()}}
    return table


def analyse_cell_size(native, done):
    """The low ground of the reference threshold on cells of other sizes; `done` holds the table of one factor."""
    table = {}
    for factor in sorted(set(CELL_FACTORS) | set(done)):
        if factor in done:
            row = done[factor]
        else:
            dem = native.coarsened(factor)
            routing, area, river, _ = prepare(dem)
            row = analyse_low_ground(dem, routing, area, river, (REFERENCE,))[0][f"{REFERENCE} km²"]
        table[f"{factor * native.cell:g} m"] = {
            "stream_length_km": row["stream_length_km"], "stream_cells_km2": row["stream_cells_km2"],
            **{level: row[level]["area_km2"] for level in LEVELS}}
    return table


def analyse_tables():
    rain = load_rainfall(ROOT / "Rainfall_Data.csv")
    events = load_events(ROOT / "Realistic_Flood_Occurrence_Data.csv")
    table = event_rainfall(events, rain)
    ranges, overlap = severity_ranges(events)
    windows = [column for column in table.columns if column.startswith("series_")]
    facts = {
        "rainfall": dict(rainfall_checks(rain), days=int(len(rain)), start=f"{rain.index[0]:%Y-%m-%d}",
                         end=f"{rain.index[-1]:%Y-%m-%d}", largest_daily_mm=rain.max().to_dict()),
        "events": {
            "events": int(len(events)), "per_region": events["region"].value_counts().to_dict(),
            "per_severity": events["severity"].value_counts().reindex(SEVERITY_ORDER).to_dict(),
            "stated_rainfall_above_series": {column: int((table["rainfall_mm"] > table[column]).sum()) for column in windows},
            "events_without_rain_in_series": {column: int((table[column] == 0).sum()) for column in windows},
            "correlation_of_stated_rainfall_with_series": {
                column: float(np.corrcoef(table["rainfall_mm"], table[column])[0, 1]) for column in windows},
            "severity_ranges": {level: {column: [float(ranges.loc[level, (column, "min")]),
                                                 float(ranges.loc[level, (column, "max")])]
                                        for column in ("area_km2", "rainfall_mm")} for level in SEVERITY_ORDER},
            "neighbouring_classes_overlap": overlap,
            "area_against_rainfall": area_rainfall_relation(events),
        },
    }
    return rain, events, table, facts


def densest_window(dem, height):
    """Rows and columns of the block of WINDOW_KM that holds the most ground within 10 m of the streams."""
    size = (int(1000 * WINDOW_KM[1] / dem.cell), int(1000 * WINDOW_KM[0] / dem.cell))
    density = ndimage.uniform_filter((np.nan_to_num(height, nan=np.inf) <= 10).astype(np.float32), size, mode="constant")
    row, col = np.unravel_index(np.argmax(density), density.shape)
    row = int(np.clip(row - size[0] // 2, 0, dem.z.shape[0] - size[0]))
    col = int(np.clip(col - size[1] // 2, 0, dem.z.shape[1] - size[1]))
    return (row, row + size[0]), (col, col + size[1])


def plot_terrain(dem, area, river, rivers, height, window):
    shown = block_mean(dem.z, 2)
    shade = hillshade(shown, 2 * dem.cell)
    extent = extent_km(dem)
    colours = ListedColormap(plt.cm.terrain(np.linspace(0.25, 1.0, 256)))
    scale = Normalize(0, 5000)
    relief = LightSource().blend_soft_light(colours(scale(np.nan_to_num(shown)))[..., :3], shade[..., None])
    relief = np.dstack([relief, np.isfinite(shown).astype(float)])

    fig, axes = plt.subplots(1, 2, figsize=(11, 8.8), sharey=True)
    axes[0].imshow(relief, extent=extent, interpolation="bilinear")
    rows, cols = np.nonzero((np.nan_to_num(area) >= REFERENCE) | river)
    east, north = dem.transform * (cols + 0.5, rows + 0.5)
    size = 0.3 + 3.0 * np.log10(np.maximum(area[rows, cols], REFERENCE) / REFERENCE)
    axes[0].scatter(east / 1000, north / 1000, s=size, color="#08519c", linewidths=0)
    for one in rivers:
        row, col = divmod(int(one["path"][int(0.55 * len(one["path"]))]), dem.z.shape[1])
        x, y = dem.position(row, col)
        axes[0].annotate(one["name"], (x / 1000, y / 1000), xytext=(8, 6), textcoords="offset points", fontsize=10,
                         fontstyle="italic", color="#08306b", path_effects=OUTLINE)
    for name, (lon, lat) in TOWNS.items():
        x, y = reproject("EPSG:4326", dem.crs, [lon], [lat])
        axes[0].plot(x[0] / 1000, y[0] / 1000, "s", color="black", markeredgecolor="white", markersize=5)
        axes[0].annotate(name, (x[0] / 1000, y[0] / 1000), xytext=(6, -3), textcoords="offset points", fontsize=9,
                         path_effects=OUTLINE)
    axes[0].set(title="Elevation and streams", xlabel="easting (km), UTM zone 45N", ylabel="northing (km)")
    fig.colorbar(plt.cm.ScalarMappable(norm=scale, cmap=colours), ax=axes[0], fraction=0.046, pad=0.02, extend="max",
                 label="elevation (m)")

    axes[1].imshow(shade, cmap="gray", vmin=-0.4, vmax=1.0, extent=extent, interpolation="bilinear")
    axes[1].imshow(height_classes(block_min(height, 4)), cmap=ListedColormap(HEIGHT_COLOURS),
                   norm=BoundaryNorm(np.arange(len(HEIGHTS) + 1) - 0.5, len(HEIGHTS)), extent=extent, interpolation="nearest")
    left, right, bottom, top = extent_km(dem, *window)
    axes[1].add_patch(Rectangle((left, bottom), right - left, top - bottom, fill=False, edgecolor="tab:red", linewidth=1.2))
    axes[1].annotate("enlarged in low_ground.png", (left, top), xytext=(0, 4), textcoords="offset points", fontsize=8,
                     color="tab:red", path_effects=OUTLINE)
    axes[1].legend(handles=class_legend(), title="height above the stream", fontsize=8, title_fontsize=8, loc="upper right")
    axes[1].set(title=f"Ground within {HEIGHTS[-1]} m above the streams\n(lowest class in each {4 * dem.cell:g} m block)",
                xlabel="easting (km), UTM zone 45N")
    fig.tight_layout()
    finish(fig, "terrain.png")


def plot_low_ground(dem, area, river, height, window, steps, curves):
    rows, cols = window
    block = (slice(*rows), slice(*cols))
    extent = extent_km(dem, rows, cols)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
    left.imshow(hillshade(dem.z[block], dem.cell), cmap="gray", vmin=-0.4, vmax=1.0, extent=extent)
    left.imshow(height_classes(height[block]), cmap=ListedColormap(HEIGHT_COLOURS),
                norm=BoundaryNorm(np.arange(len(HEIGHTS) + 1) - 0.5, len(HEIGHTS)), extent=extent, interpolation="nearest")
    stream_rows, stream_cols = np.nonzero(((np.nan_to_num(area) >= REFERENCE) | river)[block])
    x, y = dem.transform * (stream_cols + cols[0] + 0.5, stream_rows + rows[0] + 0.5)
    left.scatter(x / 1000, y / 1000, s=1.2, color="black", linewidths=0)
    left.legend(handles=class_legend([Line2D([], [], color="black", linewidth=1.5, label="stream")]),
                title="height above the stream", fontsize=8, title_fontsize=8, loc="upper right")
    left.set(xlabel="easting (km), UTM zone 45N", ylabel="northing (km)",
             title=f"Height above the fitted level of the streams from {REFERENCE} km²")

    for level in LEVELS:
        right.plot(steps, curves[(REFERENCE, level)], color=LEVEL_COLOURS[level],
                   label=f"{level} level, streams from {REFERENCE} km²")
    others = [threshold for threshold in THRESHOLDS if threshold != REFERENCE]
    for threshold, style in zip(others, ("--", ":", "-.")):
        right.plot(steps, curves[(threshold, "fitted")], color="0.45", linestyle=style,
                   label=f"fitted level, streams from {threshold} km²")
    right.set(xlabel="height above the stream (m)", ylabel="area at or below that height (km²)", xlim=(0, HEIGHTS[-1]),
              title="Low ground in the whole scene")
    right.set_ylim(bottom=0)
    right.grid(alpha=0.3)
    right.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    finish(fig, "low_ground.png")


def plot_rivers(rivers):
    fig, axes = plt.subplots(2, len(rivers), figsize=(12.5, 7.4), sharex="col", gridspec_kw={"height_ratios": [1, 1.1]})
    for column, river in enumerate(rivers):
        top, bottom = axes[0, column], axes[1, column]
        top.plot(river["distance"], river["z"], color="0.3", linewidth=0.9)
        top.set(ylabel="elevation in the DEM (m)",
                title=f"{river['name']}: {river['z'][0]:,.0f} m to {river['z'][-1]:,.0f} m in {river['distance'][-1]:.0f} km")
        top.grid(alpha=0.3)
        largest = max(np.abs(river["difference"][level]).max() for level in LEVELS)
        limit = min(60.0, 10.0 * np.ceil(1.1 * largest / 10.0))
        for level in ("carved", "filled", "fitted"):
            bottom.plot(river["distance"], river["difference"][level], color=LEVEL_COLOURS[level], linewidth=0.7, label=level)
        if largest > limit:
            bottom.text(0.99, 0.03, f"axis cut at ±{limit:.0f} m", transform=bottom.transAxes, ha="right", fontsize=8)
        bottom.set(ylim=(-limit, limit), ylabel="level given to the river minus the DEM (m)",
                   xlabel="distance downstream of the point where the river enters the scene (km)")
        bottom.grid(alpha=0.3)
        bottom.legend(fontsize=8, ncol=3, loc="upper center")
    fig.tight_layout()
    finish(fig, "river_levels.png")


def plot_resampling(native, dem, river, resampling):
    # a west-east line of cells through the river, half-way along its course in the scene
    row, col = divmod(int(river["path"][len(river["path"]) // 2]), dem.z.shape[1])
    factor = int(round(dem.cell / native.cell))
    row, col = factor * row, factor * col
    block = native.z[row - 24:row + 24, col - 40:col + 40]
    coarse = block_mean(block, 2)
    line = block.shape[0] // 2
    x = (np.arange(block.shape[1]) + 0.5) * native.cell
    styles = {"nearest": ("tab:orange", "-"), "bilinear": ("tab:green", "--"), "cubic": ("tab:blue", "-")}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), gridspec_kw={"width_ratios": [2.2, 1, 1]})
    axes[0].plot(x, block[line], "o", color="black", markersize=3, label=f"{native.cell:g} m cells of the file")
    for method, (colour, style) in styles.items():
        axes[0].plot(x, resample(coarse, 2, method)[line], color=colour, linestyle=style, linewidth=1.1,
                     drawstyle="steps-mid" if method == "nearest" else "default",
                     label=f"{method}, from {2 * native.cell:g} m means")
    axes[0].set(xlabel=f"distance along a west-east line across the {river['name']} (m)", ylabel="elevation (m)",
                title=f"Rebuilding the {native.cell:g} m cells from {2 * native.cell:g} m means")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    methods = list(resampling["methods"])
    for ax, key, label in ((axes[1], "elevation_rmse_m", "RMS error of the elevation (m)"),
                           (axes[2], "slope_rmse_deg", "RMS error of the slope (degrees)")):
        bars = ax.bar(methods, [resampling["methods"][method][key] for method in methods],
                      color=[styles[method][0] for method in methods])
        ax.bar_label(bars, fmt="%.2f", fontsize=9)
        ax.set(title=label)
        ax.margins(y=0.12)
        ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    finish(fig, "resampling.png")


def plot_rainfall(rain, facts):
    monthly = monthly_summary(rain)["total_mm"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.5, 4.4), gridspec_kw={"width_ratios": [1.5, 1]})
    width = 0.8 / len(monthly.columns)
    for position, name in enumerate(monthly.columns):
        left.bar(np.arange(12) + (position - (len(monthly.columns) - 1) / 2) * width, monthly[name], width, label=name)
    wet_months = [month for month in range(1, 13) if month not in facts["dry_months"]]
    for month, multiple in zip(wet_months, facts["monthly_multiples"][monthly.columns[0]]):
        left.text(month - 1, monthly.loc[month].max() + 8, f"×{multiple:.0f}", ha="center", fontsize=9)
    left.set_xticks(np.arange(12), MONTHS)
    left.set(ylabel="rainfall (mm)", ylim=(0, 1.1 * monthly.to_numpy().max()),
             title="Monthly rainfall; labels: multiple of the smallest monthly total above zero")
    left.grid(axis="y", alpha=0.3)
    left.legend(fontsize=8)

    monsoon = rain[rain.index.month.isin(MONSOON_MONTHS)]
    first, second = rain.columns[1], rain.columns[2]
    low, high = facts["monsoon_rank_correlation_between_districts"]
    right.scatter(monsoon[first], monsoon[second], s=14, color="tab:blue", alpha=0.6, linewidths=0)
    right.set(xlabel=f"{first} (mm per day)", ylabel=f"{second} (mm per day)",
              title=signed(f"June to September, day by day\n(rank correlation between districts {low:.2f} to {high:.2f})"))
    right.grid(alpha=0.3)
    fig.tight_layout()
    finish(fig, "rainfall_series.png")


def plot_events(events, table, relation, regions):
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.6))
    for position, region in enumerate(regions):         # the colours of the rainfall figure
        rows = table["region"] == region
        left.scatter(table.loc[rows, "series_7d_mm"], table.loc[rows, "rainfall_mm"], s=34, color=f"C{position}", label=region)
    limit = 1.05 * max(table["rainfall_mm"].max(), table["series_7d_mm"].max())
    left.plot([0, limit], [0, limit], color="0.4", linestyle="--", label="equal")
    left.set(xlim=(-8, limit), ylim=(0, limit), xlabel="rainfall of the series in the 7 days up to the event (mm)",
             ylabel="rainfall stated for the event (mm)", title="Event rainfall against the rainfall series")
    left.grid(alpha=0.3)
    left.legend(fontsize=8, loc="lower right")

    for level in SEVERITY_ORDER:
        rows = events["severity"] == level
        right.scatter(events.loc[rows, "rainfall_mm"], events.loc[rows, "area_km2"], s=34, color=SEVERITY_COLOURS[level],
                      label=f"{level} ({rows.sum()})")
    x = np.array([events["rainfall_mm"].min(), events["rainfall_mm"].max()])
    right.plot(x, relation["intercept_km2"] + relation["slope_km2_per_mm"] * x, color="black", linewidth=1,
               label=f"least squares, r = {relation['r']:.2f}")
    ranges = events.groupby("severity")["area_km2"].agg(["min", "max"]).reindex(SEVERITY_ORDER)
    for lower, upper in zip(ranges["max"].to_numpy()[:-1], ranges["min"].to_numpy()[1:]):
        right.axhline((lower + upper) / 2, color="0.5", linestyle=":", linewidth=1)
    right.set(xlabel="rainfall stated for the event (mm)", ylabel="affected area (km²)",
              title="Affected area, rainfall and severity class\n(dotted: between the areas of neighbouring classes)")
    right.grid(alpha=0.3)
    right.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    finish(fig, "flood_events.png")


def report(summary, rain, table):
    scene, sampling, terrain, low = summary["scene"], summary["resampling"], summary["terrain"], summary["low_ground"]
    rainfall, events, window = summary["rainfall"], summary["events"], summary["enlarged_map"]
    pits = summary["pits_under_the_fitted_level"]
    relation = events["area_against_rainfall"]
    lines = ["# Results", "", f"Command: `{summary['command']}`", "", "## The scene", "", "| | |", "|---|---|",
             f"| File | `{scene['file']}`, EPSG:{scene['epsg']} |",
             f"| Grid | {scene['columns']} x {scene['rows']} cells of {scene['cell_m']:g} m |",
             f"| Extent | {scene['width_km']:.1f} km by {scene['height_km']:.1f} km; {scene['west']:.2f} to {scene['east']:.2f} °E, "
             f"{scene['south']:.2f} to {scene['north']:.2f} °N |",
             f"| Area with data | {scene['area_km2']:,.0f} km² ({100 * scene['share_without_data']:.2f} % of the cells have no data) |",
             f"| Elevation | {scene['lowest_m']:.0f} to {scene['highest_m']:.0f} m, median {scene['median_m']:.0f} m; "
             f"{'stored in whole metres' if scene['whole_metres'] else 'not in whole metres'} |",
             "", f"## Rebuilding the {sampling['to_cell_m']:g} m cells from {sampling['from_cell_m']:g} m means", "",
             f"{sampling['cells']:,} cells. Neighbouring cells differ by {sampling['mean_difference_between_neighbours_m']:.2f} m "
             f"on average.", "", "| Method | RMS error of the elevation (m) | RMS error of the slope (degrees) |", "|---|---|---|"]
    lines += [f"| {method} | {row['elevation_rmse_m']:.2f} | {row['slope_rmse_deg']:.2f} |"
              for method, row in sampling["methods"].items()]
    deepest = terrain["largest_fill_at"]
    lines += ["", f"## Terrain on {terrain['cell_m']:g} m cells", "", "| | |", "|---|---|",
              f"| Area | {terrain['area_km2']:,.1f} km² in {terrain['cells']:,} cells |",
              f"| Slope | mean {terrain['slope_mean_deg']:.1f}°, median {terrain['slope_median_deg']:.1f}°; "
              f"{100 * terrain['share_flatter_than_5_deg']:.2f} % of the area is flatter than 5° and "
              f"{100 * terrain['share_steeper_than_30_deg']:.0f} % is steeper than 30° |",
              f"| Cells raised by more than 1 m when the depressions are filled | "
              f"{100 * terrain['share_raised_more_than_1_m_by_filling']:.2f} % ({terrain['area_raised_more_than_1_m_km2']:.1f} km² "
              f"in {terrain['depressions_deeper_than_1_m']:,} depressions) |",
              f"| Largest fill depth | {terrain['largest_fill_depth_m']:.1f} m at easting {deepest['easting_km']:.2f} km, "
              f"northing {deepest['northing_km']:.2f} km |",
              f"| Places where flow leaves the scene | {terrain['outlets']} |",
              "", "| Outlet | Area draining to it inside the scene (km²) | Easting (km) | Northing (km) | Longitude | Latitude | Elevation (m) |",
              "|---|---|---|---|---|---|---|"]
    lines += [f"| {position} | {row['area_km2']:,.1f} | {row['easting_km']:.2f} | {row['northing_km']:.2f} | {row['longitude']:.3f} | "
              f"{row['latitude']:.3f} | {row['elevation_m']:.0f} |" for position, row in enumerate(terrain["largest_outlets"], 1)]
    lines += ["", "## The two rivers", "",
              "| River | Enters at (easting, northing in km; elevation) | Leaves at | Length in the scene (km) | Mean fall (m/km) | "
              "Area draining to it inside the scene (km²) | Distance from the entry until "
              + ", ".join(f"{threshold}" for threshold in THRESHOLDS) + " km² drain to it inside the scene (km) |",
              "|---|---|---|---|---|---|---|"]
    for river in summary["rivers"]:
        enters, leaves = river["enters"], river["leaves"]
        lines.append(
            f"| {river['name']} | {river['side']} edge, {enters['easting_km']:.2f}, {enters['northing_km']:.2f}; "
            f"{river['enters_elevation_m']:,.0f} m | {leaves['easting_km']:.2f}, {leaves['northing_km']:.2f}; "
            f"{river['leaves_elevation_m']:,.0f} m | {river['length_km']:.1f} | {river['mean_fall_m_per_km']:.1f} | "
            f"{river['area_inside_scene_km2']:,.1f} | "
            + ", ".join(f"{value:.1f}" for value in river["km_from_entry_until_area"].values()) + " |")
    lines += ["", "Level given to the cells of each river minus their elevation in the DEM, in m:", "",
              "| River | Level | Median | 5th percentile | 95th percentile | RMS |", "|---|---|---|---|---|---|"]
    for river in summary["rivers"]:
        lines += [f"| {river['name']} | {level} | {number(row['median'])} | {number(row['percentile_5'])} | "
                  f"{number(row['percentile_95'])} | {number(row['rms'])} |" for level, row in river["level_minus_dem_m"].items()]
    lines += [""] + [f"- {river['name']}: the filled level lies {river['filled_minus_carved_median_m']:.1f} m above the carved one "
                     f"(median); the largest fill along the river is {river['largest_fill_m']:.1f} m, "
                     f"{river['largest_fill_km_from_entry']:.1f} km from the entry." for river in summary["rivers"]]
    lines += ["", "## Level of the streams", "", "Level of all stream cells minus their elevation in the DEM, in m:", "",
              "| Streams from | Level | Mean | Median | 5th percentile | 95th percentile | RMS | Largest |",
              "|---|---|---|---|---|---|---|---|"]
    for threshold, row in low.items():
        for level in LEVELS:
            d = row[level]["level_minus_dem_m"]
            lines.append(f"| {threshold} | {level} | {number(d['mean'], 2)} | {number(d['median'], 2)} | {number(d['percentile_5'])} | "
                         f"{number(d['percentile_95'])} | {number(d['rms'])} | {number(d['largest'])} |")
    lines += ["", "## Ground near the streams", "",
              "Area in km² at or below the given height above the level of the stream that the ground drains to. "
              "The stream cells themselves are included.", "",
              "| Streams from | Stream length (km) | Stream cells (km²) | Area that reaches no stream | Level | "
              + " | ".join(f"{h} m" for h in HEIGHTS) + " |", "|---|---|---|---|---|" + "---|" * len(HEIGHTS)]
    for threshold, row in low.items():
        for level in LEVELS:
            lines.append(f"| {threshold} | {row['stream_length_km']:.0f} | {row['stream_cells_km2']:.1f} | "
                         f"{100 * row['share_reaching_no_stream']:.1f} % | {level} | "
                         + " | ".join(f"{row[level]['area_km2'][f'{h} m']:.1f}" for h in HEIGHTS) + " |")
    lines += ["", f"The enlarged map covers easting {window['west_km']:.1f} to {window['east_km']:.1f} km and northing "
              f"{window['south_km']:.1f} to {window['north_km']:.1f} km: the block of {WINDOW_KM[0]:g} km by {WINDOW_KM[1]:g} km "
              f"with the most ground within 10 m of the streams ({window['area_within_10_m_km2']:.1f} km²).",
              "", f"Ground below the fitted level of its stream: {pits['area_below_the_stream_level_km2']:.1f} km². More than "
              f"{pits['depth_m']} m below it, in pits of the elevation model: {pits['area_km2']:.1f} km² in {pits['patches']} patches; "
              f"the largest, {pits['largest_patch_km2']:.1f} km², lies around easting {pits['largest_patch_easting_km']:.0f} km, "
              f"northing {pits['largest_patch_northing_km']:.0f} km. The maps show this ground in the lowest class.",
              "", f"### By river (fitted level, streams from {REFERENCE} km²)", "",
              "Area by the stream that the ground drains to, and that area divided by the length of the stream.", "",
              "| Stream | Length (km) | " + " | ".join(f"Area within {h} m (km²)" for h in HEIGHTS)
              + " | Mean width within 10 m (m) |", "|---|---|" + "---|" * (len(HEIGHTS) + 1)]
    lines += [f"| {name} | {row['length_km']:.1f} | " + " | ".join(f"{row['area_km2'][f'{h} m']:.1f}" for h in HEIGHTS)
              + f" | {row['mean_width_m']['10 m']:.0f} |" for name, row in summary["low_ground_by_river"].items()]
    lines += ["", f"### By cell size (streams from {REFERENCE} km²)", "",
              "| Cell | Stream length (km) | Stream cells (km²) | Level | " + " | ".join(f"{h} m" for h in HEIGHTS) + " |",
              "|---|---|---|---|" + "---|" * len(HEIGHTS)]
    for cell, row in summary["low_ground_by_cell_size"].items():
        lines += [f"| {cell} | {row['stream_length_km']:.0f} | {row['stream_cells_km2']:.1f} | {level} | "
                  + " | ".join(f"{row[level][f'{h} m']:.1f}" for h in HEIGHTS) + " |" for level in LEVELS]
    monthly = monthly_summary(rain)
    lines += ["", "## Rainfall table", "", f"{rainfall['days']} days, {rainfall['start']} to {rainfall['end']}.", "",
              "| Month | " + " | ".join(f"{name} (mm)" for name in rain.columns)
              + " | Multiple of the district's smallest monthly total above zero |",
              "|---|" + "---|" * (len(rain.columns) + 1)]
    wet_months = [month for month in range(1, 13) if month not in rainfall["dry_months"]]
    multiples = dict(zip(wet_months, zip(*[rainfall["monthly_multiples"][name] for name in rain.columns])))
    for month in range(1, 13):
        values = " | ".join(f"{monthly['total_mm'].loc[month, name]:.2f}" for name in rain.columns)
        multiple = ", ".join(f"{value:.3f}" for value in multiples[month]) if month in multiples else "-"
        lines.append(f"| {MONTHS[month - 1]} | {values} | {multiple} |")
    lines.append("| Year | " + " | ".join(f"{rainfall['annual_total_mm'][name]:,.1f}" for name in rain.columns) + " | |")
    low_rank, high_rank = rainfall["monsoon_rank_correlation_between_districts"]
    lines += ["", f"- Largest departure of a multiple from a whole number: "
              f"{rainfall['largest_departure_from_whole_multiples']:.4f}.",
              f"- Rank correlation of the daily values between districts, June to September: {low_rank:.2f} to {high_rank:.2f}.",
              f"- June to September days with rain in all districts: {100 * rainfall['monsoon_days_wet_in_all_districts']:.1f} % "
              f"({100 * rainfall['monsoon_days_wet_in_all_districts_if_independent']:.1f} % expected if the districts were "
              f"independent of each other).",
              "", "## Flood-event table", "",
              "| Region | Date | Stated rainfall (mm) | Series, event day (mm) | Series, 3 days (mm) | Series, 7 days (mm) |",
              "|---|---|---|---|---|---|"]
    lines += [f"| {row.region} | {row.date:%Y-%m-%d} | {row.rainfall_mm:.0f} | {row.series_1d_mm:.1f} | {row.series_3d_mm:.1f} | "
              f"{row.series_7d_mm:.1f} |" for row in table.itertuples()]
    above, correlation = events["stated_rainfall_above_series"], events["correlation_of_stated_rainfall_with_series"]
    dry = events["events_without_rain_in_series"]
    lines += ["", f"- Events whose stated rainfall is above the 7-day total of the series: {above['series_7d_mm']} of {events['events']}.",
              f"- Events for which the series has no rain: {dry['series_1d_mm']} on the event day, {dry['series_3d_mm']} in the "
              f"3 days and {dry['series_7d_mm']} in the 7 days up to it.",
              "- Correlation of the stated rainfall with the series: "
              + ", ".join(f"{correlation[f'series_{days}d_mm']:.2f} ({label})"
                          for days, label in ((1, "event day"), (3, "3 days"), (7, "7 days"))) + ".",
              "", "| Severity | Events | Affected area (km²) | Stated rainfall (mm) |", "|---|---|---|---|"]
    lines += [f"| {level} | {events['per_severity'][level]} | {row['area_km2'][0]:.0f} to {row['area_km2'][1]:.0f} | "
              f"{row['rainfall_mm'][0]:.0f} to {row['rainfall_mm'][1]:.0f} |" for level, row in events["severity_ranges"].items()]
    lines += ["", "Affected area against stated rainfall:", "", "| | |", "|---|---|",
              f"| Correlation | {relation['r']:.2f} (95 % interval {relation['r_low']:.2f} to {relation['r_high']:.2f}); "
              f"rank correlation {relation['spearman']:.2f} |",
              f"| Fitted line | {relation['slope_km2_per_mm']:.3f} km² per mm; zero area at {relation['rainfall_at_zero_area_mm']:.0f} mm |",
              f"| RMS error of the fit | {relation['rmse_fit_km2']:.1f} km² |",
              f"| RMS error for events left out of the fit | {relation['rmse_leave_one_out_km2']:.1f} km² "
              f"({relation['rmse_leave_one_out_mean_only_km2']:.1f} km² when the mean area of the other events is used) |", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factor", type=int, default=2, help="analyse on cells this many times larger than those of the file")
    parser.add_argument("--rasters", action="store_true",
                        help="also write the height above the streams and the contributing area as GeoTIFF")
    args = parser.parse_args()
    FIGURES.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)

    native = read_dem(DEM_FILE)
    dem = native.coarsened(args.factor)
    routing, area, river, entry = prepare(dem)
    low, steps, curves, kept = analyse_low_ground(dem, routing, area, river)
    height = kept["fitted"]["height"]
    rivers = analyse_rivers(dem, routing, area, entry, kept)
    resampling = analyse_resampling(native)
    rain, events, table, tables = analyse_tables()
    window = densest_window(dem, height)
    left, right, bottom, top = extent_km(dem, *window)
    block = height[slice(*window[0]), slice(*window[1])]

    command = "python run_analysis.py" + (f" --factor {args.factor}" if args.factor != parser.get_default("factor") else "")
    summary = {"command": command, "scene": describe_scene(native), "resampling": resampling,
               "terrain": analyse_terrain(dem, routing, area), "rivers": [one["facts"] for one in rivers], "low_ground": low,
               "low_ground_by_river": analyse_by_river(dem, routing, area, river, rivers, height),
               "low_ground_by_cell_size": analyse_cell_size(native, {args.factor: low[f"{REFERENCE} km²"]}),
               "pits_under_the_fitted_level": analyse_pits(dem, height),
               "enlarged_map": {"west_km": left, "east_km": right, "south_km": bottom, "north_km": top,
                                "area_within_10_m_km2": float((block <= 10).sum() * dem.cell_area_km2)},
               **tables}
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    (RESULTS / "summary.md").write_text(report(summary, rain, table), encoding="utf-8")
    if args.rasters:
        write_raster(RESULTS / "height_above_stream_m.tif", height, dem)
        write_raster(RESULTS / "contributing_area_km2.tif", area, dem)

    plot_terrain(dem, area, river, rivers, height, window)
    plot_low_ground(dem, area, river, height, window, steps, curves)
    plot_rivers(rivers)
    plot_resampling(native, dem, rivers[0], resampling)
    plot_rainfall(rain, tables["rainfall"])
    plot_events(events, table, tables["events"]["area_against_rainfall"], list(rain.columns))
    print(f"wrote {RESULTS / 'summary.md'} and the figures in {FIGURES}")


if __name__ == "__main__":
    main()
