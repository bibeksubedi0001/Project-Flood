import numpy as np
import rasterio
from scipy.interpolate import RectBivariateSpline

def resample_spline(input_path, output_path, new_resolution):
    # Open the source DEM file
    with rasterio.open(input_path) as src:
        # Read DEM data and metadata
        dem_data = src.read(1, masked=True)  # Read the first band
        src_transform = src.transform

        # Get the original grid
        x = np.arange(0, dem_data.shape[1]) * src_transform[1] + src_transform[0]
        y = np.arange(0, dem_data.shape[0]) * src_transform[5] + src_transform[3]

        # Create the new grid for resampling
        x_new = np.arange(x[0], x[-1] + new_resolution, new_resolution)
        y_new = np.arange(y[0], y[-1] + new_resolution, new_resolution)

        # Perform spline interpolation
        spline = RectBivariateSpline(y, x, dem_data.filled(np.nan))
        z_new = spline(y_new, x_new)

        # Update metadata for the new raster
        new_transform = rasterio.transform.from_origin(
            x_new[0], y_new[0], new_resolution, new_resolution
        )
        profile = src.profile
        profile.update({
            "transform": new_transform,
            "width": len(x_new),
            "height": len(y_new),
            "dtype": "float32",
            "compress": "lzw",
        })

        # Save the resampled DEM
        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(z_new.astype("float32"), 1)

    print(f"Resampled DEM saved at: {output_path}")

# Define input and output paths
input_dem_path = r"D:\Flood_Prediction\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1.dem.tif"
output_dem_path = r"D:\Flood_Prediction\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1\Resample_Spline.dem.tif"

# Set new resolution (6.25m for both x and y)
new_resolution = 6.25

# Call the function
resample_spline(input_dem_path, output_dem_path, new_resolution)
4