import rasterio
from rasterio.enums import Resampling

def resample_nearest_neighbor(input_path, output_path, new_resolution):
    # Open the source DEM file
    with rasterio.open(input_path) as src:
        # Calculate the scaling factor
        scale_x = src.res[0] / new_resolution
        scale_y = src.res[1] / new_resolution

        # Calculate new dimensions
        new_width = int(src.width * scale_x)
        new_height = int(src.height * scale_y)

        # Update transform for the new resolution
        new_transform = src.transform * src.transform.scale(
            1 / scale_x, 1 / scale_y
        )

        # Update metadata for the new raster
        profile = src.profile
        profile.update({
            "transform": new_transform,
            "width": new_width,
            "height": new_height,
            "compress": "lzw",
        })

        # Perform resampling
        with rasterio.open(output_path, "w", **profile) as dst:
            for i in range(1, src.count + 1):  # For each band
                resampled_data = src.read(
                    i,
                    out_shape=(src.count, new_height, new_width),
                    resampling=Resampling.nearest  # Nearest Neighbor Resampling
                )
                dst.write(resampled_data, i)

    print(f"Resampled DEM saved at: {output_path}")

# Define input and output paths
input_dem_path = r"D:\Flood_Prediction\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1.dem.tif"
output_dem_path = r"D:\Flood_Prediction\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1\Resample_Nearest.dem.tif"

# Set new resolution (6.25m for both x and y)
new_resolution = 6.25

# Call the function
resample_nearest_neighbor(input_dem_path, output_dem_path, new_resolution)