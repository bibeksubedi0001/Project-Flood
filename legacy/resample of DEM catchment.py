import rasterio
from rasterio.enums import Resampling

def resample_dem(input_path, output_path, scale_factor):
    # Open the source DEM file
    with rasterio.open(input_path) as src:
        # Calculate new transform and dimensions
        transform = src.transform
        new_transform = transform * transform.scale(1 / scale_factor, 1 / scale_factor)
        new_width = int(src.width * scale_factor)
        new_height = int(src.height * scale_factor)

        # Create metadata for the resampled output
        profile = src.profile
        profile.update({
            'transform': new_transform,
            'width': new_width,
            'height': new_height
        })

        # Perform resampling
        with rasterio.open(output_path, 'w', **profile) as dst:
            for i in range(1, src.count + 1):
                resampled_data = src.read(
                    i,
                    out_shape=(new_height, new_width),
                    resampling=Resampling.bilinear
                )
                dst.write(resampled_data, i)

# Define input and output paths
input_dem_path = r'D:\Flood_Prediction\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1.dem.tif'
output_dem_path = r'D:\Flood_Prediction\AP_27418_PLR_F0540_RT1\AP_27418_PLR_F0540_RT1\Resample.dem.tif'
scale_factor = 12.5 / 6.25  # Scaling to achieve 6.25m resolution

# Call the function
resample_dem(input_dem_path, output_dem_path, scale_factor)
