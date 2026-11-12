import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Load the Rainfall Data
rainfall_path = "D:\Flood_Prediction\Rainfall_Data.csv"
try:
    rainfall_df = pd.read_csv(rainfall_path)
    # Ensure the Date column is in datetime format
    rainfall_df['Date'] = pd.to_datetime(rainfall_df['Date'], errors='coerce')
    # Check for invalid dates and drop them
    if rainfall_df['Date'].isnull().sum() > 0:
        print("Warning: Invalid dates detected. Cleaning the Date column.")
        rainfall_df = rainfall_df.dropna(subset=['Date'])
    # Create the Month column
    rainfall_df['Month'] = rainfall_df['Date'].dt.month
    print("Rainfall Data Loaded Successfully:\n")
    print(rainfall_df.head())
except FileNotFoundError:
    print(f"Error: Rainfall data not found at {rainfall_path}. Ensure the file exists.")

# Step 2: Load the Flood Occurrence Data
flood_path = "D:\Flood_Prediction\Realistic_Flood_Occurrence_Data.csv"
try:
    flood_df = pd.read_csv(flood_path)
    print("\nFlood Occurrence Data Loaded Successfully:\n")
    print(flood_df.head())
except FileNotFoundError:
    print(f"Error: Flood occurrence data not found at {flood_path}. Ensure the file exists.")

# Step 3: Handle Missing Values
if 'rainfall_df' in locals():
    print("\nMissing Values in Rainfall Data Before Cleaning:")
    print(rainfall_df.isnull().sum())

    # Replace missing values with the monthly mean for each region
    rainfall_df['Sankhuwasabha'] = rainfall_df['Sankhuwasabha'].fillna(
        rainfall_df.groupby('Month')['Sankhuwasabha'].transform('mean')
    )
    rainfall_df['Khotang'] = rainfall_df['Khotang'].fillna(
        rainfall_df.groupby('Month')['Khotang'].transform('mean')
    )
    rainfall_df['Okhaldhunga'] = rainfall_df['Okhaldhunga'].fillna(
        rainfall_df.groupby('Month')['Okhaldhunga'].transform('mean')
    )

    print("\nMissing Values in Rainfall Data After Cleaning:")
    print(rainfall_df.isnull().sum())

if 'flood_df' in locals():
    print("\nMissing Values in Flood Data Before Cleaning:")
    print(flood_df.isnull().sum())

    # Replace missing values with the median for numeric columns
    flood_df['Rainfall (mm)'] = flood_df['Rainfall (mm)'].fillna(flood_df['Rainfall (mm)'].median())
    flood_df['Affected Area (sq km)'] = flood_df['Affected Area (sq km)'].fillna(
        flood_df['Affected Area (sq km)'].median()
    )

    print("\nMissing Values in Flood Data After Cleaning:")
    print(flood_df.isnull().sum())

# Step 4: Feature Engineering
if 'rainfall_df' in locals():
    # Create total monthly rainfall for each region
    monthly_rainfall = rainfall_df.groupby('Month')[['Sankhuwasabha', 'Khotang', 'Okhaldhunga']].sum().reset_index()
    print("\nTotal Monthly Rainfall:")
    print(monthly_rainfall)

    # Create a binary column for rainy days
    rainfall_df['Rainy Day (Sankhuwasabha)'] = rainfall_df['Sankhuwasabha'] > 0
    rainfall_df['Rainy Day (Khotang)'] = rainfall_df['Khotang'] > 0
    rainfall_df['Rainy Day (Okhaldhunga)'] = rainfall_df['Okhaldhunga'] > 0

if 'flood_df' in locals():
    # Encode Flood Severity into numeric values
    flood_df['Flood Severity Numeric'] = flood_df['Flood Severity'].map({'High': 3, 'Moderate': 2, 'Low': 1})

    # Calculate the ratio of Affected Area to Rainfall
    flood_df['Area to Rainfall Ratio'] = flood_df['Affected Area (sq km)'] / flood_df['Rainfall (mm)']
    print("\nFlood Data with New Features:")
    print(flood_df.head())

# Step 5: Save Cleaned Data
if 'rainfall_df' in locals():
    rainfall_df.to_csv("D:\Flood_Prediction\Rainfall_Data_update.csv", index=False)
    print("\nCleaned Rainfall Data Saved Successfully.")

if 'flood_df' in locals():
    flood_df.to_csv("D:\Flood_Prediction\Realistic_Flood_Occurrence_Data_final.csv", index=False)
    print("\nCleaned Flood Data Saved Successfully.")
