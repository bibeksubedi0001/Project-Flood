import matplotlib.pyplot as plt
import pandas as pd
# Step 1: Load the Rainfall Data
rainfall_path = "D:\Flood_Prediction\Rainfall_Data.csv"
rainfall_df = pd.read_csv(rainfall_path)
# Step 2: Load the Flood Occurrence Data
flood_path = "D:\Flood_Prediction\Realistic_Flood_Occurrence_Data.csv"
flood_df = pd.read_csv(flood_path)
# Flood severity counts
severity_counts = flood_df['Flood Severity'].value_counts()
print("Flood Severity Counts:")
print(severity_counts)

# Plot flood severity
severity_counts.plot(kind='bar', figsize=(6, 4), color=['red', 'orange', 'green'])
plt.title("Flood Severity Distribution")
plt.xlabel("Severity Level")
plt.ylabel("Number of Events")
plt.show()

# Correlation analysis
correlation = flood_df[['Rainfall (mm)', 'Affected Area (sq km)']].corr()
print("\nCorrelation Between Rainfall and Affected Area:")
print(correlation)

# Scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(flood_df['Rainfall (mm)'], flood_df['Affected Area (sq km)'], alpha=0.7)
plt.title("Rainfall vs. Affected Area")
plt.xlabel("Rainfall (mm)")
plt.ylabel("Affected Area (sq km)")
plt.grid()
plt.show()
