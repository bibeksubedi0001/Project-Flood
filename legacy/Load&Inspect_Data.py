import pandas as pd

# Step 1: Load the Rainfall Data
rainfall_path = "D:\Flood_Prediction\Rainfall_Data.csv"
try:
    rainfall_df = pd.read_csv(rainfall_path)
    print("Rainfall Data Loaded Successfully:\n")
    print(rainfall_df.head())
    print("\nRainfall Data Info:\n")
    print(rainfall_df.info())
except FileNotFoundError:
    print(f"Error: Rainfall data not found at {rainfall_path}. Ensure the file exists.")

# Step 2: Load the Flood Occurrence Data
flood_path = "D:\Flood_Prediction\Realistic_Flood_Occurrence_Data.csv"
try:
    flood_df = pd.read_csv(flood_path)
    print("\nFlood Occurrence Data Loaded Successfully:\n")
    print(flood_df.head())
    print("\nFlood Occurrence Data Info:\n")
    print(flood_df.info())
except FileNotFoundError:
    print(f"Error: Flood occurrence data not found at {flood_path}. Ensure the file exists.")

# Step 3: Check for Missing Values
print("\nMissing Values in Rainfall Data:\n")
print(rainfall_df.isnull().sum() if 'rainfall_df' in locals() else "Rainfall data not loaded.")

print("\nMissing Values in Flood Occurrence Data:\n")
print(flood_df.isnull().sum() if 'flood_df' in locals() else "Flood occurrence data not loaded.")

# Step 4: Summarize Observations
if 'rainfall_df' in locals():
    print("\nRainfall Summary:\n")
    print(rainfall_df.describe())

if 'flood_df' in locals():
    print("\nFlood Occurrence Summary:\n")
    print(flood_df.describe())
