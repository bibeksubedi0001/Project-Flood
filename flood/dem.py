"""Reading a DEM, changing its grid, and terrain attributes that need no flow routing."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import rasterio
from rasterio.warp import transform_bounds
from scipy import ndimage

RESAMPLING_ORDER = {"nearest": 0, "bilinear": 1, "cubic": 3}


@dataclass
class Dem:
    z: np.ndarray                 # elevation in metres, NaN outside the scene
    cell: float                   # cell size in metres
    transform: rasterio.Affine
    crs: object

    @property
    def cell_area_km2(self):
        return self.cell**2 / 1e6

    def bounds(self):
        """Left, bottom, right and top in the units of the coordinate system."""
        rows, cols = self.z.shape
        left, top = self.transform * (0, 0)
        right, bottom = self.transform * (cols, rows)
        return left, bottom, right, top

    def bounds_lonlat(self):
        return transform_bounds(self.crs, "EPSG:4326", *self.bounds())

    def coarsened(self, factor):
        """The same DEM with cells `factor` times larger, each the mean of the cells it covers."""
        if factor == 1:
            return self
        return Dem(block_mean(self.z, factor), self.cell * factor, self.transform * rasterio.Affine.scale(factor), self.crs)

    def position(self, row, col):
        """Easting and northing of the centre of a cell."""
        return self.transform * (col + 0.5, row + 0.5)


def block_mean(z, factor):
    """Mean over blocks of factor x factor cells; a block with a missing cell is missing."""
    rows, cols = (np.array(z.shape) // factor) * factor
    blocks = z[:rows, :cols].reshape(rows // factor, factor, cols // factor, factor)
    return blocks.mean(axis=(1, 3))


def read_dem(path):
    """Read band 1 of a raster as elevation in metres, with NaN for the no-data value."""
    with rasterio.open(path) as source:
        z = source.read(1).astype(np.float64)
        if source.nodata is not None:
            z[z == source.nodata] = np.nan
        if not np.isclose(source.res[0], source.res[1]):
            raise ValueError("the cells of the raster are not square")
        return Dem(z, float(source.res[0]), source.transform, source.crs)


def resample(z, factor, method="bilinear"):
    """Elevation on a grid with cells `factor` times smaller, by nearest, bilinear or cubic interpolation."""
    return ndimage.zoom(z, factor, order=RESAMPLING_ORDER[method], mode="nearest", grid_mode=True)


def resampling_test(z, cell, factor=2, border=8):
    """Rebuild a grid from its block means with each method and compare the result with the grid itself.

    Returns, per method, the RMS error of the elevation in metres and of the slope in degrees, leaving out
    `border` cells at the edges where the interpolation has no neighbours.
    """
    if np.isnan(z).any():
        raise ValueError("the grid must have no missing cells")
    coarse = block_mean(z, factor)
    truth = z[:coarse.shape[0] * factor, :coarse.shape[1] * factor]
    inner = (slice(border, -border), slice(border, -border))
    true_slope = np.degrees(np.arctan(slope(truth, cell)))[inner]
    errors = {}
    for method in RESAMPLING_ORDER:
        rebuilt = resample(coarse, factor, method)
        rebuilt_slope = np.degrees(np.arctan(slope(rebuilt, cell)))[inner]
        errors[method] = {
            "elevation_rmse_m": float(np.sqrt(np.mean((rebuilt[inner] - truth[inner]) ** 2, dtype=np.float64))),
            "slope_rmse_deg": float(np.sqrt(np.mean((rebuilt_slope - true_slope) ** 2, dtype=np.float64))),
        }
    return errors


def slope(z, cell):
    """Tangent of the slope angle by the method of Horn (1981); NaN next to missing cells."""
    p = np.pad(z, 1, constant_values=np.nan)
    a, b, c = p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:]
    d, f = p[1:-1, :-2], p[1:-1, 2:]
    g, h, i = p[2:, :-2], p[2:, 1:-1], p[2:, 2:]
    dz_dx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * cell)
    dz_dy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * cell)
    return np.hypot(dz_dx, dz_dy)


def write_raster(path, values, dem):
    """Write a grid on the cells of `dem` as a compressed single-band GeoTIFF."""
    profile = {"driver": "GTiff", "height": values.shape[0], "width": values.shape[1], "count": 1,
               "dtype": "float32", "crs": dem.crs, "transform": dem.transform, "nodata": np.nan,
               "compress": "deflate", "predictor": 3, "tiled": True}
    with rasterio.open(path, "w", **profile) as target:
        target.write(values.astype(np.float32), 1)
