import xarray as xr
import pandas as pd
import numpy as np

def read_strat_H2O(file_path, var_name='combinedh2oq', 
                                     lat_min=-23.5, lat_max=23.5, 
                                     level_min=90, level_max=80):
    """
    Reads stratospheric H2O data from a NetCDF file and calculates the area-weighted
    average over the tropics and within a specified altitude range, indexed by timestamp.

    Parameters:
    - file_path: str, the path to the NetCDF file.
    - var_name: str, the specific data variable name in the file.
    - lat_min: float, the minimum latitude for defining the tropics.
    - lat_max: float, the maximum latitude for defining the tropics.
    - level_min: float, the minimum pressure level (in hPa) for the altitude range.
    - level_max: float, the maximum pressure level (in hPa) for the altitude range.

    Returns:
    - A pandas DataFrame with the area-weighted average strat H2O data over the tropics
      and within the specified altitude range, indexed by timestamp.
    """
    # Load the dataset
    ds = xr.open_dataset(file_path)
    
    # Select the data variable and focus on the specified level and latitude range
    data_var = ds[var_name].sel(level=slice(level_min, level_max), lat=slice(lat_min, lat_max))
    
    # Calculate the cosine of latitude for area weighting
    lat_rad = np.deg2rad(data_var.lat)  # Convert latitudes to radians
    weights = np.cos(lat_rad)
    weights /= weights.mean()  # Normalize weights so they sum to 1 across the selected latitude range
    
    # Calculate the area-weighted average over the tropics and specified level range
    tropical_avg = data_var.weighted(weights).mean(dim=['lat', 'lon', 'level'])
    
    # Convert to a pandas DataFrame for easy use, indexed by time
    tropical_avg_df = tropical_avg.to_dataframe().dropna().reset_index()
    
    # Set the time as the DataFrame index
    tropical_avg_df['time'] = tropical_avg_df['time'].dt.to_period('M').dt.to_timestamp()
    tropical_avg_df.set_index('time', inplace=True)
    tropical_avg_df.rename(columns={var_name: 'strat_H2O'}, inplace=True)
    
    return tropical_avg_df


def calc_heating_rate_series(series):
    """
    Calculates the heating rate from a pandas Series of temperatures,
    ensuring the output Series is indexed by time.
    """
    delta_T = series.diff()  # Change in temperature
    delta_t = series.index.to_series().diff().dt.total_seconds() / (3600 * 24)  # Change in time in days
    heating_rate = delta_T / delta_t  # Heating rate in K/day
    return heating_rate.rename('BDO_ind')  # Rename the series

def read_reanalysis(file_path, temp_var='t'):
    ds = xr.open_dataset(file_path)
    weights = np.cos(np.deg2rad(ds.latitude))
    
    # Calculate the mean temperature between 70-80 hPa and directly convert to Series
    temp_70_80 = ds[temp_var].sel(level=70).weighted(weights).mean(dim=['latitude', 'longitude']).to_series()
    bdo_df = calc_heating_rate_series(temp_70_80).to_frame()    
    # Calculate the mean temperature at 500 hPa
    temp_500 = ds[temp_var].sel(level=500).weighted(weights).mean(dim=['latitude', 'longitude']).to_series()
    # Convert to a pandas Series to ensure we have a time-indexed Series
    deltaT_df = temp_500.to_frame()
    print(deltaT_df.head())
    deltaT_df.rename(columns={temp_var: 'DeltaT'}, inplace=True)

    return bdo_df, deltaT_df



def read_qbo_indexes(filename, pressure_level):
    """
    Reads a file containing zonal wind data and returns a time series of wind speeds
    at a specified pressure level, accounting for empty lines and headers.

    Parameters:
    - filename: str, the name of the file containing the data.
    - pressure_level: str, the pressure level (hPa) of interest as a string.

    Returns:
    - A pandas DataFrame containing the time series of wind speeds at the specified pressure level.
    """
    with open(filename, 'r') as file:
        data = file.readlines()
    
    year = None
    monthly_data = []
    for line in data:
        line = line.strip()
        if not line:  # Skip empty lines
            continue
        if line.isdigit():
            year = line  # Found a year line
        else:
            parts = line.split()
            if parts[0] == 'hPa':  # Skip the header line
                continue
            if parts[0] == pressure_level:
                for month_idx, value in enumerate(parts[1:], 1):
                    monthly_data.append({
                        'Year': int(year),
                        'Month': month_idx,  # Keep as numeric month for now
                        'WindSpeed': float(value) / 10  # Convert to m/s
                    })
    
    df = pd.DataFrame(monthly_data)
    
    # Convert 'Month' to datetime for better handling and readability
    df['Date'] = pd.to_datetime(df[['Year', 'Month']].assign(DAY=1))
    
    # Optional: Drop the now redundant 'Year' and 'Month' columns
    df = df.drop(['Year', 'Month'], axis=1)
    
    # Set 'Date' as the index
    df.set_index('Date', inplace=True)
    df.rename(columns={'WindSpeed': 'QBO_ind'}, inplace=True)
    
    return df


def combine_dataframes(df_dict, time_col=None):
    """
    Combines separate DataFrames (each representing a different variable) into a single DataFrame.
    
    Parameters:
    - df_dict: dict, a dictionary where each key is a variable name and each value is a DataFrame 
               containing a timeseries for that variable.
    - time_col: str, optional, the name of the column to use as the DateTimeIndex if the DataFrames 
                do not already have a DateTimeIndex. If None, it is assumed the index is already a DateTimeIndex.
    
    Returns:
    - A single pandas DataFrame with each variable as a column, aligned by the common DateTimeIndex.
    """
    combined_df = pd.DataFrame()
    
    for var_name, df in df_dict.items():
        # If a time column is specified, set it as the index (assumes datetime format is convertible)
        if time_col is not None and time_col in df.columns:
            df = df.set_index(pd.to_datetime(df[time_col]))
            df.drop(columns=[time_col], inplace=True)
        
        # If combining for the first time, initialize combined_df with the current df
        if combined_df.empty:
            combined_df = df[[var_name]]  # Use double brackets to ensure DataFrame structure
        else:
            # Join the current df with the combined_df on the index
            combined_df = combined_df.join(df[var_name], how='outer')
    
    return combined_df


import statsmodels.api as sm

import pandas as pd
import statsmodels.api as sm

import pandas as pd
import statsmodels.api as sm

def regress_variables(df, dependent_var, independent_vars, lags={'QBO_ind': 3, 'BDO_ind': 1, 'DeltaT': 1}):
    """
    Performs multivariate linear regression on a specified dependent variable against
    a set of lagged independent variables.

    Parameters:
    - df: pandas DataFrame containing the columns for the dependent and independent variables,
          indexed by a DateTimeIndex.
    - dependent_var: str, the name of the dependent variable column in the DataFrame.
    - independent_vars: list of str, the names of the independent variable columns in the DataFrame.
    - lags: dict, a dictionary specifying the lag for each independent variable.
    
    Returns:
    - The regression results summary.
    """
    # Prepare DataFrame by lagging independent variables as specified
    df_lagged = df.copy()
    for var in independent_vars:
        if var in lags:
            df_lagged[f'{var}_lagged'] = df[var].shift(lags[var])
        else:
            df_lagged[f'{var}_lagged'] = df[var]  # No lag if not specified

    # Remove rows with any NaN values resulting from the lagging
    df_clean = df_lagged.dropna()

    # Define the dependent variable (Y) and independent variables (X) with lagged names
    Y = df_clean[dependent_var]
    X = df_clean[[f'{var}_lagged' for var in independent_vars]]
    X = sm.add_constant(X)  # Adds a constant term to the predictors

    # Perform the linear regression
    model = sm.OLS(Y, X).fit()

    # Return the regression results summary
    return model.summary()



def calculate_anomalies(df):
    """
    Calculates anomalies for each column in the DataFrame with respect to the monthly mean
    across all years included in the dataset.

    Parameters:
    - df: pandas DataFrame, the combined DataFrame with each variable as a column, 
          indexed by a DateTimeIndex.

    Returns:
    - A pandas DataFrame with the same structure as the input, where each value is replaced
      by the anomaly from the monthly mean for its respective month.
    """
    # Copy the DataFrame to avoid modifying the original data
    anomalies_df = df.copy()
    
    # Loop through each column to calculate anomalies
    for column in df.columns:
        # Calculate the monthly mean across all years
        monthly_mean = df[column].groupby(df.index.month).mean()
        
        # Subtract the monthly mean from each value to get the anomaly
        anomalies_df[column] = df[column] - df.index.month.map(monthly_mean)
    
    return anomalies_df


strat_h2o_file = '/Users/newtonnguyen/Documents/projects/strat-trop-box-model/strat_data/swoosh-v02.6-198401-202212-lonlatpress-20deg-5deg-L31.nc'
strat_h2o_df = read_strat_H2O(strat_h2o_file)

reanalysis_file = '/Users/newtonnguyen/Documents/projects/OHMethane/data/atm_dynamics/ERA5_temperature.nc'
bdo_df, deltaT_df = read_reanalysis(reanalysis_file)


qbo_filename = '/Users/newtonnguyen/Documents/projects/OHMethane/data/atm_dynamics/qbo.dat'
pressure_level = "50"  # The pressure level you're interested in, as a string
qbo_df = read_qbo_indexes(qbo_filename, pressure_level)

combined_data_dict = {'QBO_ind' : qbo_df, 'BDO_ind' : bdo_df, 'DeltaT' : deltaT_df, 'strat_H2O' : strat_h2o_df}
combined_df = combine_dataframes(combined_data_dict)
combined_df.dropna(inplace=True)
print(combined_df.head())
#combined_df = calculate_anomalies(combined_df)

dependent_var = 'strat_H2O'
independent_vars = ['QBO_ind', 'BDO_ind', 'DeltaT']
results = regress_variables(combined_df, dependent_var, independent_vars)

print(results)