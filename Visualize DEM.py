import rasterio
import matplotlib.pyplot as plt
import numpy as np

# Corrected DEM file path
dem_path = r"D:/Flood_Prediction/AP_27418_PLR_F0540_RT1/AP_27418_PLR_F0540_RT1/AP_27418_PLR_F0540_RT1.dem.tif.ovr"

# Open the DEM file
try:
    with rasterio.open(dem_path) as dem:
        # Read the DEM data
        dem_data = dem.read(1)  # Read the first band
        print("DEM Metadata:")
        print(dem.meta)
        print("\nBounds:")
        print(dem.bounds)
        print("\nCRS (Coordinate Reference System):")
        print(dem.crs)

        # Handle nodata values
        dem_data = np.where(dem_data == dem.nodata, np.nan, dem_data)

        # Plot DEM using imshow
        plt.figure(figsize=(10, 6))
        im = plt.imshow(dem_data, cmap="terrain", extent=dem.bounds, origin="upper")
        plt.colorbar(im, label="Elevation (m)")
        plt.title("Digital Elevation Model (DEM)")
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        plt.show()

except FileNotFoundError:
    print(f"Error: DEM file not found at {dem_path}. Ensure the file exists.")
