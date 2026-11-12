# Results

Command: `python run_analysis.py`

## The scene

| | |
|---|---|
| File | `AP_27418_PLR_F0540_RT1.dem.tif`, EPSG:32645 |
| Grid | 3281 x 5729 cells of 12.5 m |
| Extent | 41.0 km by 71.6 km; 86.43 to 86.85 °E, 27.06 to 27.71 °N |
| Area with data | 2,921 km² (0.54 % of the cells have no data) |
| Elevation | 189 to 5806 m, median 1746 m; stored in whole metres |

## Rebuilding the 12.5 m cells from 25 m means

17,034,849 cells. Neighbouring cells differ by 4.41 m on average.

| Method | RMS error of the elevation (m) | RMS error of the slope (degrees) |
|---|---|---|
| nearest | 3.98 | 2.00 |
| bilinear | 1.09 | 1.56 |
| cubic | 0.55 | 0.92 |

## Terrain on 25 m cells

| | |
|---|---|
| Area | 2,920.0 km² in 4,672,013 cells |
| Slope | mean 27.8°, median 27.4°; 0.97 % of the area is flatter than 5° and 40 % is steeper than 30° |
| Cells raised by more than 1 m when the depressions are filled | 0.89 % (26.0 km² in 2,848 depressions) |
| Largest fill depth | 232.8 m at easting 470.92 km, northing 3059.84 km |
| Places where flow leaves the scene | 175 |

| Outlet | Area draining to it inside the scene (km²) | Easting (km) | Northing (km) | Longitude | Latitude | Elevation (m) |
|---|---|---|---|---|---|---|
| 1 | 1,886.2 | 443.92 | 3003.41 | 86.434 | 27.152 | 293 |
| 2 | 619.3 | 471.57 | 2993.59 | 86.713 | 27.064 | 202 |
| 3 | 125.5 | 443.92 | 3053.69 | 86.432 | 27.606 | 1668 |
| 4 | 64.4 | 443.92 | 3015.81 | 86.433 | 27.264 | 513 |
| 5 | 63.3 | 443.92 | 3026.34 | 86.433 | 27.359 | 879 |

## The two rivers

| River | Enters at (easting, northing in km; elevation) | Leaves at | Length in the scene (km) | Mean fall (m/km) | Area draining to it inside the scene (km²) | Distance from the entry until 5, 25, 100 km² drain to it inside the scene (km) |
|---|---|---|---|---|---|---|
| Dudh Koshi | north edge, 471.52, 3065.09; 2,430 m | 443.92, 3003.41; 293 m | 99.6 | 21.4 | 1,886.2 | 0.8, 3.8, 11.7 |
| Sun Koshi | west edge, 443.92, 3001.16; 275 m | 471.57, 2993.59; 202 m | 41.1 | 1.8 | 619.3 | 1.9, 7.1, 11.9 |

Level given to the cells of each river minus their elevation in the DEM, in m:

| River | Level | Median | 5th percentile | 95th percentile | RMS |
|---|---|---|---|---|---|
| Dudh Koshi | filled | 2.2 | 0.0 | 17.6 | 11.1 |
| Dudh Koshi | fitted | 0.0 | -7.5 | 7.1 | 7.9 |
| Dudh Koshi | carved | -3.5 | -29.5 | 0.0 | 23.2 |
| Sun Koshi | filled | 2.5 | 0.0 | 9.5 | 4.6 |
| Sun Koshi | fitted | -0.2 | -3.9 | 4.5 | 2.7 |
| Sun Koshi | carved | -7.8 | -20.0 | -0.8 | 10.0 |

- Dudh Koshi: the filled level lies 7.5 m above the carved one (median); the largest fill along the river is 184.8 m, 6.4 km from the entry.
- Sun Koshi: the filled level lies 10.5 m above the carved one (median); the largest fill along the river is 23.0 m, 30.8 km from the entry.

## Level of the streams

Level of all stream cells minus their elevation in the DEM, in m:

| Streams from | Level | Mean | Median | 5th percentile | 95th percentile | RMS | Largest |
|---|---|---|---|---|---|---|---|
| 5 km² | filled | 2.03 | 0.00 | 0.0 | 9.2 | 8.3 | 216.5 |
| 5 km² | fitted | 0.00 | 0.00 | -3.2 | 3.1 | 5.0 | 130.1 |
| 5 km² | carved | -2.70 | 0.00 | -12.8 | 0.0 | 10.3 | -216.5 |
| 25 km² | filled | 2.98 | 0.75 | 0.0 | 12.2 | 7.4 | 184.8 |
| 25 km² | fitted | 0.00 | 0.00 | -4.6 | 4.5 | 4.7 | 149.6 |
| 25 km² | carved | -4.66 | -1.00 | -19.0 | 0.0 | 13.0 | -184.8 |
| 100 km² | filled | 3.99 | 1.75 | 0.0 | 15.0 | 9.2 | 184.8 |
| 100 km² | fitted | 0.00 | 0.00 | -5.7 | 5.8 | 6.0 | 149.6 |
| 100 km² | carved | -6.88 | -2.75 | -23.8 | 0.0 | 16.9 | -184.8 |

## Ground near the streams

Area in km² at or below the given height above the level of the stream that the ground drains to. The stream cells themselves are included.

| Streams from | Stream length (km) | Stream cells (km²) | Area that reaches no stream | Level | 5 m | 10 m | 20 m | 50 m |
|---|---|---|---|---|---|---|---|---|
| 5 km² | 865 | 18.2 | 2.1 % | filled | 56.5 | 75.4 | 106.5 | 187.2 |
| 5 km² | 865 | 18.2 | 2.1 % | fitted | 49.4 | 69.1 | 101.0 | 182.3 |
| 5 km² | 865 | 18.2 | 2.1 % | carved | 39.2 | 59.9 | 93.4 | 176.1 |
| 25 km² | 408 | 8.5 | 5.5 % | filled | 38.4 | 48.6 | 64.9 | 105.8 |
| 25 km² | 408 | 8.5 | 5.5 % | fitted | 32.7 | 43.7 | 60.6 | 102.0 |
| 25 km² | 408 | 8.5 | 5.5 % | carved | 23.4 | 35.4 | 53.8 | 96.4 |
| 100 km² | 230 | 4.7 | 9.9 % | filled | 28.7 | 35.2 | 45.5 | 70.3 |
| 100 km² | 230 | 4.7 | 9.9 % | fitted | 24.1 | 31.2 | 42.0 | 67.2 |
| 100 km² | 230 | 4.7 | 9.9 % | carved | 15.6 | 23.8 | 35.9 | 62.0 |

The enlarged map covers easting 443.9 to 455.9 km and northing 2997.8 to 3006.8 km: the block of 12 km by 9 km with the most ground within 10 m of the streams (6.9 km²).

Ground below the fitted level of its stream: 12.2 km². More than 10 m below it, in pits of the elevation model: 3.3 km² in 235 patches; the largest, 1.0 km², lies around easting 471 km, northing 3060 km. The maps show this ground in the lowest class.

### By river (fitted level, streams from 25 km²)

Area by the stream that the ground drains to, and that area divided by the length of the stream.

| Stream | Length (km) | Area within 5 m (km²) | Area within 10 m (km²) | Area within 20 m (km²) | Area within 50 m (km²) | Mean width within 10 m (m) |
|---|---|---|---|---|---|---|
| Dudh Koshi | 99.6 | 10.4 | 12.7 | 16.6 | 25.4 | 128 |
| Sun Koshi | 41.1 | 7.4 | 9.6 | 12.5 | 18.0 | 233 |
| other streams | 266.9 | 15.0 | 21.4 | 31.6 | 58.5 | 80 |

### By cell size (streams from 25 km²)

| Cell | Stream length (km) | Stream cells (km²) | Level | 5 m | 10 m | 20 m | 50 m |
|---|---|---|---|---|---|---|---|
| 12.5 m | 411 | 4.2 | filled | 38.7 | 48.8 | 64.9 | 105.7 |
| 12.5 m | 411 | 4.2 | fitted | 33.7 | 44.4 | 61.3 | 102.4 |
| 12.5 m | 411 | 4.2 | carved | 23.1 | 35.4 | 54.1 | 96.7 |
| 25 m | 408 | 8.5 | filled | 38.4 | 48.6 | 64.9 | 105.8 |
| 25 m | 408 | 8.5 | fitted | 32.7 | 43.7 | 60.6 | 102.0 |
| 25 m | 408 | 8.5 | carved | 23.4 | 35.4 | 53.8 | 96.4 |
| 50 m | 402 | 16.9 | filled | 39.8 | 49.1 | 66.0 | 107.7 |
| 50 m | 402 | 16.9 | fitted | 34.0 | 43.8 | 61.4 | 103.4 |
| 50 m | 402 | 16.9 | carved | 27.5 | 37.6 | 56.1 | 98.7 |

## Rainfall table

366 days, 2024-01-01 to 2024-12-31.

| Month | Sankhuwasabha (mm) | Khotang (mm) | Okhaldhunga (mm) | Multiple of the district's smallest monthly total above zero |
|---|---|---|---|---|
| Jan | 0.00 | 0.00 | 0.00 | - |
| Feb | 0.00 | 0.00 | 0.00 | - |
| Mar | 80.01 | 75.01 | 70.01 | 1.000, 1.000, 1.000 |
| Apr | 160.00 | 150.02 | 139.98 | 2.000, 2.000, 2.000 |
| May | 240.00 | 225.01 | 209.99 | 3.000, 3.000, 3.000 |
| Jun | 320.00 | 300.00 | 280.00 | 4.001, 4.000, 4.001 |
| Jul | 479.99 | 449.99 | 419.99 | 6.001, 6.000, 6.001 |
| Aug | 639.98 | 599.99 | 560.01 | 8.001, 8.000, 8.001 |
| Sep | 319.98 | 300.01 | 280.01 | 4.000, 4.000, 4.001 |
| Oct | 160.00 | 149.98 | 139.99 | 2.000, 2.000, 2.000 |
| Nov | 79.99 | 75.00 | 69.99 | 1.000, 1.000, 1.000 |
| Dec | 0.00 | 0.00 | 0.00 | - |
| Year | 2,480.0 | 2,325.0 | 2,170.0 | |

- Largest departure of a multiple from a whole number: 0.0013.
- Rank correlation of the daily values between districts, June to September: -0.16 to -0.06.
- June to September days with rain in all districts: 2.5 % (4.8 % expected if the districts were independent of each other).

## Flood-event table

| Region | Date | Stated rainfall (mm) | Series, event day (mm) | Series, 3 days (mm) | Series, 7 days (mm) |
|---|---|---|---|---|---|
| Sankhuwasabha | 2024-07-15 | 319 | 0.0 | 0.0 | 133.6 |
| Khotang | 2024-08-10 | 189 | 0.0 | 0.0 | 33.0 |
| Okhaldhunga | 2024-09-05 | 138 | 0.0 | 0.0 | 0.0 |
| Sankhuwasabha | 2024-06-20 | 240 | 0.0 | 25.8 | 62.6 |
| Khotang | 2024-07-30 | 229 | 0.0 | 0.0 | 0.0 |
| Okhaldhunga | 2024-08-25 | 74 | 0.0 | 39.8 | 99.2 |
| Sankhuwasabha | 2024-09-10 | 100 | 0.0 | 100.7 | 100.7 |
| Khotang | 2024-10-05 | 167 | 0.0 | 0.0 | 63.5 |
| Okhaldhunga | 2024-07-25 | 237 | 0.0 | 111.4 | 160.5 |
| Sankhuwasabha | 2024-08-15 | 307 | 41.5 | 114.8 | 228.9 |
| Khotang | 2024-09-20 | 179 | 0.0 | 0.0 | 37.7 |
| Okhaldhunga | 2024-07-05 | 105 | 0.0 | 69.4 | 115.0 |
| Sankhuwasabha | 2024-08-20 | 218 | 26.8 | 26.8 | 90.0 |
| Khotang | 2024-09-25 | 235 | 0.0 | 0.0 | 8.0 |
| Okhaldhunga | 2024-07-10 | 74 | 0.0 | 0.0 | 0.0 |
| Sankhuwasabha | 2024-08-30 | 260 | 79.1 | 157.5 | 247.7 |
| Khotang | 2024-09-15 | 249 | 8.0 | 18.2 | 111.6 |
| Okhaldhunga | 2024-08-01 | 273 | 0.0 | 0.0 | 12.0 |
| Sankhuwasabha | 2024-08-18 | 180 | 0.0 | 0.0 | 114.8 |
| Khotang | 2024-09-12 | 63 | 1.9 | 93.4 | 131.9 |

- Events whose stated rainfall is above the 7-day total of the series: 16 of 20.
- Events for which the series has no rain: 15 on the event day, 10 in the 3 days and 3 in the 7 days up to it.
- Correlation of the stated rainfall with the series: 0.37 (event day), 0.00 (3 days), 0.27 (7 days).

| Severity | Events | Affected area (km²) | Stated rainfall (mm) |
|---|---|---|---|
| Low | 6 | 6 to 19 | 63 to 138 |
| Moderate | 7 | 20 to 46 | 167 to 249 |
| High | 7 | 50 to 96 | 218 to 319 |

Affected area against stated rainfall:

| | |
|---|---|
| Correlation | 0.77 (95 % interval 0.50 to 0.90); rank correlation 0.82 |
| Fitted line | 0.281 km² per mm; zero area at 48 mm |
| RMS error of the fit | 17.8 km² |
| RMS error for events left out of the fit | 19.6 km² (29.3 km² when the mean area of the other events is used) |
