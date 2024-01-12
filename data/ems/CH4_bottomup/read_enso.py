import pandas as pd

# Read the space-separated file
df = pd.read_csv('enso_meiv2.txt', sep='\s+', skiprows=1, header=None)

# Create an empty DataFrame with columns for year and each month
columns = ['Year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# List to store rows
rows = []

# Iterate over the rows 
for i in range(0, len(df)):
    year = df.iloc[i, 0]
    monthly_values = df.iloc[i, 1:].tolist()
    row = [year] + monthly_values
    rows.append(row)

# Convert the list of rows to a DataFrame
data = pd.DataFrame(rows, columns=columns)

# Set 'Year' column as index
data.set_index('Year', inplace=True)

# Compute annual average
data['annual_avg'] = data.drop('annual_avg', errors='ignore').mean(axis=1)  # Drop 'annual_avg' to avoid including it in the calculation, if it exists

# Define the range
start_year = 1984
end_year = 2022

# Select years within the range and display
subset = data.loc[start_year:end_year, 'annual_avg']
subset.to_csv('enso_index.csv')
print(data)
