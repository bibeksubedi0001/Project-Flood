"""Checks of the terrain functions on surfaces with known answers, and of the checks of the two tables."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy import ndimage, optimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from flood import (Dem, area_rainfall_relation, block_mean, event_rainfall, height_above_drainage, load_events,
                   load_rainfall, monthly_summary, rainfall_checks, read_dem, resample, resampling_test, route_flow,
                   severity_ranges, slope, stream_level, value_at_drainage, write_raster)
from flood.hydro import _monotone_fit


def valley(rows=30, cols=21, centre=10):
    """A V-shaped valley that falls gently towards the last row; the channel is the column `centre`."""
    y, x = np.mgrid[0:rows, 0:cols]
    return 1.0 * np.abs(x - centre) + 0.001 * (rows - 1 - y)


class Grid(unittest.TestCase):
    def test_block_mean(self):
        z = np.arange(24.0).reshape(4, 6)
        np.testing.assert_allclose(block_mean(z, 2), [[3.5, 5.5, 7.5], [15.5, 17.5, 19.5]])
        z[0, 0] = np.nan
        self.assertTrue(np.isnan(block_mean(z, 2)[0, 0]))
        self.assertEqual(block_mean(np.ones((5, 7)), 2).shape, (2, 3))

    def test_slope_of_a_plane_is_exact(self):
        y, x = np.mgrid[0:12, 0:15]
        tangent = slope(0.03 * 25 * x - 0.04 * 25 * y, 25.0)
        np.testing.assert_allclose(tangent[1:-1, 1:-1], 0.05, rtol=1e-12)
        self.assertTrue(np.isnan(tangent[0]).all())

    def test_resampling_keeps_constants_and_planes(self):
        self.assertEqual(resample(np.ones((6, 8)), 2, "cubic").shape, (12, 16))
        for method in ("nearest", "bilinear", "cubic"):
            np.testing.assert_allclose(resample(np.full((6, 8), 7.0), 2, method), 7.0)
        y, x = np.mgrid[0:60, 0:60]
        errors = resampling_test(2.0 * x + 3.0 * y, 1.0)
        self.assertAlmostEqual(errors["bilinear"]["elevation_rmse_m"], 0.0, places=9)
        self.assertAlmostEqual(errors["bilinear"]["slope_rmse_deg"], 0.0, places=9)
        self.assertLess(errors["cubic"]["elevation_rmse_m"], 0.01)
        # a fine cell lies half a cell from its block centre in both directions: errors of 2.5 and 0.5
        self.assertAlmostEqual(errors["nearest"]["elevation_rmse_m"], np.sqrt(3.25), places=9)

    def test_resampling_test_rejects_missing_cells(self):
        z = np.ones((40, 40))
        z[5, 5] = np.nan
        with self.assertRaises(ValueError):
            resampling_test(z, 1.0)

    def test_raster_round_trip_and_coarsening(self):
        values = np.arange(48.0).reshape(6, 8)
        values[0, 0] = np.nan
        dem = Dem(values, 25.0, rasterio.Affine(25.0, 0.0, 400000.0, 0.0, -25.0, 3000000.0),
                  rasterio.crs.CRS.from_epsg(32645))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "dem.tif"
            write_raster(path, dem.z, dem)
            again = read_dem(path)
        np.testing.assert_array_equal(again.z, dem.z)
        self.assertEqual(again.cell, 25.0)
        self.assertEqual(again.bounds(), (400000.0, 2999850.0, 400200.0, 3000000.0))
        self.assertEqual(dem.position(0, 0), (400012.5, 2999987.5))
        coarse = dem.coarsened(2)
        self.assertEqual(coarse.cell, 50.0)
        self.assertEqual(coarse.bounds(), (400000.0, 2999850.0, 400200.0, 3000000.0))
        np.testing.assert_array_equal(coarse.z, block_mean(values, 2))
        self.assertAlmostEqual(coarse.cell_area_km2, 0.0025)
        west, south, east, north = dem.bounds_lonlat()
        self.assertTrue(85.9 < west < east < 86.1 and 27.0 < south < north < 27.2)


class Routing(unittest.TestCase):
    def test_plane_drains_along_the_slope(self):
        z = np.tile(np.arange(8.0)[::-1], (6, 1))
        routing = route_flow(z)
        np.testing.assert_array_equal(routing.filled, z)
        np.testing.assert_array_equal(routing.cells, np.tile(np.arange(1.0, 9.0), (6, 1)))

    def test_every_cell_reaches_an_outlet_on_rough_ground_with_gaps(self):
        z = 100 * ndimage.gaussian_filter(np.random.default_rng(0).normal(size=(80, 90)), 3)
        z[10:14, 20:30] = np.nan
        z[0, :5] = np.nan
        valid = np.isfinite(z)
        routing = route_flow(z)
        self.assertEqual(routing.cells.ravel()[routing.outlets()].sum(), valid.sum())
        self.assertTrue(np.all(routing.filled[valid] >= z[valid]))
        self.assertGreater((routing.filled[valid] > z[valid]).sum(), 0)
        self.assertTrue(np.isnan(routing.cells[~valid]).all())
        draining = np.flatnonzero(routing.receiver >= 0)
        level = routing.filled.ravel()
        self.assertTrue(np.all(level[routing.receiver[draining]] <= level[draining]))
        self.assertTrue(valid.ravel()[routing.receiver[draining]].all())

    def test_pit_is_filled_to_its_rim(self):
        y, x = np.mgrid[0:21, 0:21]
        cone = np.hypot(x - 10, y - 10)
        routing = route_flow(cone)
        inside = cone < 10
        np.testing.assert_allclose(routing.filled[inside], 10.0)
        np.testing.assert_array_equal(routing.filled[~inside], cone[~inside])
        self.assertEqual(routing.cells.ravel()[routing.outlets()].sum(), cone.size)

    def test_flat_ground_drains_without_loops(self):
        z = np.zeros((15, 17))
        z[:, 0] = -1.0
        routing = route_flow(z)
        self.assertEqual(routing.cells.ravel()[routing.outlets()].sum(), z.size)

    def test_valley_accumulation_and_height_above_the_channel(self):
        rows, cols, centre = 30, 21, 10
        routing = route_flow(valley(rows, cols, centre))
        np.testing.assert_array_equal(routing.cells[:, centre], cols * np.arange(1, rows + 1))
        hand = height_above_drainage(routing, routing.cells >= cols)
        np.testing.assert_allclose(hand, np.abs(np.arange(cols) - centre) * np.ones((rows, 1)), atol=1e-9)

    def test_flow_path_runs_down_the_side_and_along_the_channel(self):
        rows, cols, centre = 30, 21, 10
        routing = route_flow(valley(rows, cols, centre))
        outlet = routing.outlets()[0]
        self.assertEqual(divmod(int(outlet), cols), (rows - 1, centre))
        path = routing.flow_path(5 * cols + 3)
        np.testing.assert_array_equal(path[:8], 5 * cols + np.arange(3, 11))
        np.testing.assert_array_equal(path[8:], np.arange(6, rows) * cols + centre)
        self.assertEqual(path[-1], outlet)
        np.testing.assert_array_equal(routing.step_length(path), [1.0] * (len(path) - 1) + [0.0])

    def test_closed_edge_turns_the_flow_inwards(self):
        z = np.tile(np.arange(8.0), (6, 1))              # falls to the west edge
        open_edge = route_flow(z)
        self.assertTrue(np.all(open_edge.outlets() % 8 == 0))
        closed = np.zeros(z.shape, bool)
        closed[:, 0] = True
        routing = route_flow(z, closed)
        np.testing.assert_array_equal(np.sort(routing.outlets()), [1, 5 * 8 + 1])
        np.testing.assert_array_equal(routing.filled[:, 0], 1.0)
        self.assertEqual(routing.cells.ravel()[routing.outlets()].sum(), z.size)
        with self.assertRaises(ValueError):
            route_flow(z, np.ones(z.shape, bool))

    def test_height_is_missing_where_no_stream_is_reached(self):
        routing = route_flow(np.tile(np.arange(8.0)[::-1], (6, 1)))
        stream = np.zeros((6, 8), bool)
        stream[2, 5:] = True
        hand = height_above_drainage(routing, stream)
        np.testing.assert_array_equal(hand[2], [5, 4, 3, 2, 1, 0, 0, 0])
        self.assertTrue(np.isnan(np.delete(hand, 2, axis=0)).all())

    def test_every_cell_takes_the_value_of_the_stream_cell_it_drains_to(self):
        rows, cols, centre = 30, 21, 10
        routing = route_flow(valley(rows, cols, centre))
        stream = routing.cells >= cols                   # the channel
        row_number = np.tile(np.arange(rows)[:, None], (1, cols))
        # every cell runs along its row to the channel
        np.testing.assert_array_equal(value_at_drainage(routing, stream, row_number), row_number)
        upper = stream & (row_number < 10)
        found = value_at_drainage(routing, upper, row_number)
        np.testing.assert_array_equal(found[:10], row_number[:10])
        self.assertTrue(np.isnan(found[10:]).all())

    def test_ground_without_barriers_is_left_as_it_is(self):
        z = valley()
        routing = route_flow(z)
        np.testing.assert_array_equal(routing.carved, z)
        stream = routing.cells >= z.shape[1]
        for level in ("filled", "fitted", "carved"):
            levels = stream_level(routing, stream, level)
            np.testing.assert_allclose(levels[stream], z[stream], rtol=0, atol=1e-12)
            self.assertTrue(np.isnan(levels[~stream]).all())
            np.testing.assert_allclose(height_above_drainage(routing, stream, level),
                                       height_above_drainage(routing, stream, "filled"), rtol=0, atol=1e-12)
        with self.assertRaises(ValueError):
            height_above_drainage(routing, stream, "smoothed")

    def test_fitted_level_pools_a_pit_with_the_cells_downstream(self):
        # one stream of four cells, listed from its top: 5, 1, 4, 3; the outlet is the last
        z = np.array([5.0, 1.0, 4.0, 3.0])
        receiver = np.array([1, 2, 3, -1])
        level = _monotone_fit(np.array([3, 2, 1, 0]), receiver, z, np.ones(4, bool))
        np.testing.assert_allclose(level, [5.0, 8 / 3, 8 / 3, 8 / 3])
        # cells outside the stream take no part
        level = _monotone_fit(np.array([3, 2, 1, 0]), receiver, z, np.array([False, True, True, True]))
        np.testing.assert_allclose(level[1:], 8 / 3)
        self.assertTrue(np.isnan(level[0]))

    def test_fitted_level_is_the_least_squares_solution_on_a_tree(self):
        for seed in range(6):
            rng = np.random.default_rng(seed)
            count = 16
            receiver = np.array([-1] + [int(rng.integers(0, k)) for k in range(1, count)])
            z = rng.normal(size=count).round(2)
            level = _monotone_fit(np.arange(count), receiver, z, np.ones(count, bool))
            self.assertTrue(np.all(level[1:] >= level[receiver[1:]] - 1e-12))
            constraints = [{"type": "ineq", "fun": lambda f, k=k: f[k] - f[receiver[k]]} for k in range(1, count)]
            best = optimize.minimize(lambda f: np.sum((f - z) ** 2), np.full(count, z.mean()), constraints=constraints,
                                     method="SLSQP", options={"ftol": 1e-12, "maxiter": 500})
            np.testing.assert_allclose(level, best.x, atol=1e-4)
            self.assertLessEqual(np.sum((level - z) ** 2), best.fun + 1e-9)

    def test_carving_cuts_a_barrier_where_filling_makes_a_lake(self):
        rows, cols, centre = 30, 21, 10
        z = valley(rows, cols, centre)
        z[:2] += 10.0                                    # a head wall: the valley drains only past the dam
        z[15, centre - 2:centre + 3] += 5.0              # the dam; the lowest way round it is at 3.014
        routing = route_flow(z)
        np.testing.assert_allclose(routing.filled[2:15, centre - 2:centre + 3], 3.014)
        self.assertTrue(np.all(routing.carved <= z))
        self.assertGreaterEqual((z - routing.carved > 1.0).sum(), 1)
        draining = np.flatnonzero(routing.receiver >= 0)
        carved = routing.carved.ravel()
        self.assertTrue(np.all(carved[routing.receiver[draining]] <= carved[draining]))
        stream = routing.cells >= 5 * cols
        fitted = stream_level(routing, stream, "fitted").ravel()
        on_stream = np.flatnonzero(stream.ravel() & (routing.receiver >= 0))
        self.assertTrue(np.all(fitted[routing.receiver[on_stream]] <= fitted[on_stream] + 1e-12))
        self.assertTrue(np.all(fitted[stream.ravel()] <= routing.filled[stream] + 1e-12))
        self.assertTrue(np.all(fitted[stream.ravel()] >= routing.carved[stream] - 1e-12))
        on_filled, on_fitted, on_carved = (height_above_drainage(routing, stream, level)
                                           for level in ("filled", "fitted", "carved"))
        self.assertGreaterEqual(np.nanmin(on_carved), 0.0)
        self.assertLess(np.nanmin(on_filled), -2.9)      # the floor of the lake lies below the level of its outlet
        for height in (on_filled, on_fitted, on_carved):
            self.assertTrue(np.all(height[stream] == 0.0))
        reached = np.isfinite(on_fitted)
        self.assertTrue(np.all(on_filled[reached] <= on_fitted[reached] + 1e-12))
        self.assertTrue(np.all(on_fitted[reached] <= on_carved[reached] + 1e-12))
        # the lake is low ground against the filled level; against the carved one only the valley floor is low
        self.assertGreater((on_filled[2:15] <= 0.5).sum(), 2 * (on_carved[2:15] <= 0.5).sum())
        np.testing.assert_allclose(on_carved[2:10, centre], 0.0, atol=0.02)
        np.testing.assert_allclose(on_carved[2:10, centre + 2], 2.0, atol=0.02)


class Tables(unittest.TestCase):
    def test_rainfall_file_must_have_every_day(self):
        days = pd.date_range("2024-01-01", "2024-01-10").delete(4)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "rain.csv"
            pd.DataFrame({"Date": days.strftime("%m/%d/%Y"), "A": 1.0}).to_csv(path, index=False)
            with self.assertRaises(ValueError):
                load_rainfall(path)

    def test_rainfall_checks_find_whole_multiples_and_dry_months(self):
        days = pd.date_range("2024-01-01", "2024-12-31")
        multiple = np.array([0, 0, 1, 2, 3, 4, 6, 8, 4, 2, 1, 0])[days.month - 1]
        rain = pd.DataFrame({"A": 80.0 * multiple / days.days_in_month, "B": 70.0 * multiple / days.days_in_month},
                            index=days)
        checks = rainfall_checks(rain)
        self.assertEqual(checks["dry_months"], [1, 2, 12])
        self.assertEqual(checks["annual_total_mm"], {"A": 2480.0, "B": 2170.0})
        np.testing.assert_allclose(checks["monthly_multiples"]["B"], [1, 2, 3, 4, 6, 8, 4, 2, 1])
        self.assertLess(checks["largest_departure_from_whole_multiples"], 1e-9)
        summary = monthly_summary(rain)
        self.assertAlmostEqual(summary["total_mm"].loc[8, "A"], 640.0)
        self.assertEqual(summary["wet_days"].loc[8, "A"], 31)
        self.assertEqual(summary["wet_days"].loc[1, "A"], 0)

    def test_event_rainfall_sums_the_days_up_to_the_event(self):
        rain = pd.DataFrame({"A": np.arange(1.0, 11.0)}, index=pd.date_range("2024-01-01", "2024-01-10"))
        events = pd.DataFrame({"region": ["A"], "date": pd.to_datetime(["2024-01-05"]), "rainfall_mm": [40.0]})
        table = event_rainfall(events, rain, windows=(1, 3))
        self.assertEqual(table.loc[0, "series_1d_mm"], 5.0)
        self.assertEqual(table.loc[0, "series_3d_mm"], 12.0)
        with self.assertRaises(ValueError):
            event_rainfall(events.assign(region="B"), rain)

    def test_severity_ranges_report_overlap(self):
        events = pd.DataFrame({"severity": ["Low", "Low", "Moderate", "Moderate", "High", "High"],
                               "area_km2": [1, 5, 6, 9, 10, 20], "rainfall_mm": [10, 50, 40, 70, 80, 90]})
        table, overlap = severity_ranges(events)
        self.assertEqual(list(table.index), ["Low", "Moderate", "High"])
        self.assertEqual(overlap, {"area_km2": False, "rainfall_mm": True})

    def test_line_is_recovered_and_left_out_events_are_predicted_worse(self):
        x = 10.0 * np.arange(1, 13)
        y = 0.5 * x - 2 + 0.3 * np.tile([1.0, -1.0], 6)
        relation = area_rainfall_relation(pd.DataFrame({"rainfall_mm": x, "area_km2": y}))
        self.assertAlmostEqual(relation["slope_km2_per_mm"], 0.5, delta=0.01)
        self.assertAlmostEqual(relation["rainfall_at_zero_area_mm"], 4.0, delta=1.0)
        self.assertGreater(relation["r"], 0.999)
        self.assertGreater(relation["rmse_leave_one_out_km2"], relation["rmse_fit_km2"])
        self.assertLess(relation["rmse_leave_one_out_km2"], relation["rmse_leave_one_out_mean_only_km2"])

    def test_tables_of_the_repository(self):
        rain = load_rainfall(ROOT / "Rainfall_Data.csv")
        events = load_events(ROOT / "Realistic_Flood_Occurrence_Data.csv")
        checks = rainfall_checks(rain)
        self.assertEqual(checks["dry_months"], [1, 2, 12])
        self.assertLess(checks["largest_departure_from_whole_multiples"], 0.002)
        self.assertLess(checks["monsoon_rank_correlation_between_districts"][1], 0.0)
        table = event_rainfall(events, rain)
        self.assertEqual(int((table["rainfall_mm"] > table["series_7d_mm"]).sum()), 16)
        self.assertEqual(severity_ranges(events)[1], {"area_km2": False, "rainfall_mm": True})
        self.assertAlmostEqual(area_rainfall_relation(events)["r"], 0.769, places=3)

    def test_bad_event_tables_are_rejected(self):
        events = pd.read_csv(ROOT / "Realistic_Flood_Occurrence_Data.csv")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "events.csv"
            events.drop(columns="Region").to_csv(path, index=False)
            with self.assertRaises(ValueError):
                load_events(path)
            events.assign(**{"Flood Severity": "Severe"}).to_csv(path, index=False)
            with self.assertRaises(ValueError):
                load_events(path)


if __name__ == "__main__":
    unittest.main()
