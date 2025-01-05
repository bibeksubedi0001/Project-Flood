import matplotlib.pyplot as plt
import pandas as pd

# Step 1: Load the Rainfall Data
rainfall_path = "D:\Flood_Prediction\Rainfall_Data.csv"
rainfall_df = pd.read_csv(rainfall_path)
# Step 2: Load the Flood Occurrence Data
flood_path = "D:\Flood_Prediction\Realistic_Flood_Occurrence_Data.csv"
flood_df = pd.read_csv(flood_path)
# Convert Date column to datetime if not already done
rainfall_df['Date'] = pd.to_datetime(rainfall_df['Date'])

# Extract the month from the Date column
rainfall_df['Month'] = rainfall_df['Date'].dt.month

# Group by Month and calculate the sum for numeric columns only
monthly_rainfall = rainfall_df.groupby('Month')[['Sankhuwasabha', 'Khotang', 'Okhaldhunga']].sum()

# Plot monthly trends
monthly_rainfall.plot(kind='bar', figsize=(10, 6))
plt.title("Monthly Rainfall Trends")
plt.xlabel("Month")
plt.ylabel("Total Rainfall (mm)")
plt.legend(["Sankhuwasabha", "Khotang", "Okhaldhunga"])
plt.show()
