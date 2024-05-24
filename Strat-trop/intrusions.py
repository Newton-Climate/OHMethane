"""
-- Newton Nguyen

Script for calculating tracer-tracer correlations
Classifies intrusion and non-intrusion events from HIPPO Data
Classification is based on ozone gradient

Functions are defined first. 
The bottom of the script calls the functions.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def read_hippo(filename, round_coordinates=True, round_to=1):
    """
    Function to read the HIPPO aircraft atmospheric chemistry observations.
    Observations are point samples every 10 secs.

    Note: obs timestamps are rounded down to the first of month
    option to round lon and lat to nearest degree or specified value.

    Parameters:
    * filename (str): path of HIPPO file
    * round_coordinates (bool): if True, then round to nearest lon and lat degree or specified value
    * round_to (float): the value to which lon and lat should be rounded (default is 1 degree)

    Returns:
    hippo_df (pd.DataFrame): the processed HIPPO files with the 
    correct columns (time, lon, lat, CH4, N2O, O3)
    """
    hippo_df = pd.read_csv(filename, delim_whitespace=True)
    # Create a timestamp column that is pd.DateTime from Year and DOY
    hippo_df['time'] = pd.to_datetime(hippo_df['Year'].astype(str) + '-' + hippo_df['DOY'].astype(str), format='%Y-%j')
    hippo_df['time'] = hippo_df['time'].dt.to_period('M').dt.to_timestamp()

    # Extract relevant columns
    subset_cols = ['H.no', 'time', 'GGLON', 'GGLAT', 'PSXC', 'CH4_QCLS', 'N2O_QCLS', 'O3_ppb', 'H2O_UWV', 'CO.X']
    hippo_df = hippo_df[subset_cols]
    hippo_df.dropna(subset=subset_cols, inplace=True)
    # Rename columns to more generic names
    hippo_df.rename(columns={'GGLON': 'lon', 'GGLAT': 'lat', 'PSXC': 'p', 'CH4_QCLS': 'CH4', 'N2O_QCLS': 'N2O', 'O3_ppb': 'O3', 'H2O_UWV': 'H2O', 'CO.X': 'CO'}, inplace=True)

    # --- #
    # Round lon and lat to specified value if round_coordinates is True
    # --- #
    if round_coordinates:
        hippo_df['lon'] = np.round(hippo_df['lon'] / round_to) * round_to
        hippo_df['lat'] = np.round(hippo_df['lat'] / round_to) * round_to
    # Round pressure to the nearest 100 hPa
    hippo_df['p'] = np.round(hippo_df['p'] / 100) * 100

    return hippo_df


def calculate_vertical_gradient(data):
    """
    Calculate the vertical gradient of O3 with respect to pressure.
    
    Parameters:
    * data (pd.DataFrame): DataFrame with pressure and O3 concentration

    Returns:
    pd.DataFrame: DataFrame with vertical gradient of O3 added in-place
    """
    # Sort by pressure within each (lon, lat) group
    data = data.sort_values(by=['lon', 'lat', 'p'])
    # Calculate the vertical gradient
    data['gradient_O3'] = data.groupby(['lon', 'lat'])['O3'].diff() / data.groupby(['lon', 'lat'])['p'].diff()

    return data


def detect_intrusions(data, threshold):
    """
    Detect potential stratospheric intrusion events based on the O3 gradient threshold.
    
    Parameters:
    * data (pd.DataFrame): DataFrame with calculated gradients
    * threshold (float): The threshold for detecting intrusions (should be negative)
    Units are ppb/HPa

    Returns:
    pd.DataFrame: DataFrame with intrusion events marked
    """
    data['is_intrusion'] = data['gradient_O3'].abs() > threshold

    return data


def filter_colocate_data(data):
    """
    Filter data to ensure intrusion and non-intrusion events are co-located.
    Algorithm overview:
    1. find unique spatial coordinates in lon, lat, and p
    2. for each uniq cordinate, find the data that is in that coordinate.
        Note that spatial coordinates have been rounded 
        so we will have multiple obs in one coordinate
    3. Take the first obs that is in the unique coordinate
        note that data has already been classified 
        into intrusion and non-intrusion events

    
    Parameters:
    * data (pd.DataFrame): DataFrame with intrusion events marked

    Returns:
    pd.DataFrame: DataFrame with co-located data
    """
    unique_locations = data[['lon', 'lat', 'p']].drop_duplicates()
    intrusion_data = []
    non_intrusion_data = []

    for loc in unique_locations.itertuples(index=False):
        # subset the entire dataset by unique location
        loc_data = data[(data['lon'] == loc.lon) & (data['lat'] == loc.lat) & (data['p'] == loc.p)]

# --- #
        # extract rows that have been classified as unique locations
        # save each unique location as intrusion or non-intrusion dataframes
        # there will be multiple rows here
        # due to sampling rate from aircraft being so high
        # --- #
        intrusion_loc_data = loc_data[loc_data['is_intrusion']]
        non_intrusion_loc_data = loc_data[~loc_data['is_intrusion']]
        
        # save the first row of observations in colocated data
        if not intrusion_loc_data.empty and not non_intrusion_loc_data.empty:
            intrusion_data.append(intrusion_loc_data)
            non_intrusion_data.append(non_intrusion_loc_data)

    intrusion_data = pd.concat(intrusion_data)
    non_intrusion_data = pd.concat(non_intrusion_data)
    
    return intrusion_data, non_intrusion_data


def filter_dataframe(df, criteria):
    """
    Filter the DataFrame based on the criteria specified in a dictionary.

    Parameters:
    * df (pd.DataFrame): The DataFrame to filter
    * criteria (dict): Dictionary specifying the filtering criteria
    key is the specific variable of interest
    value (either single or tupple) is the filter criteria

    Returns:
    pd.DataFrame: The filtered DataFrame
    """
    for key, value in criteria.items():
        if isinstance(value, tuple) and len(value) == 2:
            # Range filter
            df = df[(df[key] >= value[0]) & (df[key] <= value[1])]
        elif isinstance(value, (int, float)):
            # Single value filter
            df = df[df[key] == value]
        else:
            raise ValueError("Invalid filter value: must be a tuple (for range) or a single value")

    return df



def plot_correlations(intrusion_data, non_intrusion_data, tracers):
    """
    Plot correlations between chemical tracers for intrusion and non-intrusion events.
    
    Parameters:
    * intrusion_data (pd.DataFrame): DataFrame with intrusion event data
    * non_intrusion_data (pd.DataFrame): DataFrame with non-intrusion event data
    * tracers (list): List of chemical tracer columns to analyze
    """
    intrusion_corr = intrusion_data[tracers].corr()
    non_intrusion_corr = non_intrusion_data[tracers].corr()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sns.heatmap(intrusion_corr, ax=axes[0], annot=True, cmap='coolwarm')
    axes[0].set_title('Intrusion Event Correlations')

    sns.heatmap(non_intrusion_corr, ax=axes[1], annot=True, cmap='coolwarm')
    axes[1].set_title('Non-Intrusion Event Correlations')

    # plt.show()



def calculate_correlation(data, tracer_1, tracer_2):
    """
    Calculate the correlation between two tracers.
    
    Parameters:
    * data (pd.DataFrame): DataFrame with the tracer data
    * tracer_1 (str): The first tracer column name
    * tracer_2 (str): The second tracer column name

    Returns:
    float: The correlation coefficient between the two tracers
    """
    return data[[tracer_1, tracer_2]].corr().iloc[0, 1]

def plot_correlation_comparison(intrusion_data, non_intrusion_data, tracer_1, tracer_2):
    """
    Plot the comparison of correlations between two tracers during intrusion and non-intrusion events.
    
    Parameters:
    * intrusion_data (pd.DataFrame): DataFrame with intrusion event data
    * non_intrusion_data (pd.DataFrame): DataFrame with non-intrusion event data
    * tracer_1 (str): The first tracer column name
    * tracer_2 (str): The second tracer column name
    """
    # Calculate correlations
    correlation_intrusion = calculate_correlation(intrusion_data, tracer_1, tracer_2)
    correlation_non_intrusion = calculate_correlation(non_intrusion_data, tracer_1, tracer_2)

    print(f'Correlations for {tracer_1} and {tracer_2}: \n')
    print(f"Correlation during intrusion events: {correlation_intrusion}")
    print(f"Correlation not during intrusion events: {correlation_non_intrusion}")

    # Plot comparison
    plt.figure(figsize=(10, 5))
    plt.scatter(intrusion_data[tracer_1], intrusion_data[tracer_2], color='red', label='Intrusion Events')
    plt.scatter(non_intrusion_data[tracer_1], non_intrusion_data[tracer_2], color='blue', label='Non-Intrusion Periods')
    plt.xlabel(f'{tracer_1} (ppbv)')
    plt.ylabel(f'{tracer_2} (ppbv)')
    plt.legend()
    plt.title(f'Comparison of {tracer_1} versus {tracer_2} during stratospheric intrusions')
    plt.savefig(f'figs/intrusion_comparison_{tracer_1}_{tracer_2}.pdf')
    # plt.show()

filename = "../data/obs/aircraft/hippo_merged/HIPPO_all_missions_merge_10s_20121129.tbl"
threshold = -10 # ppb O3 / HPa
tracers = ['CH4', 'N2O', 'O3', 'H2O', 'CO']  # List of chemical tracers
# --- #
# Read and preprocess data
# --- #
hippo_df = read_hippo(filename, round_coordinates=True, round_to=0.5)
# Calculate vertical gradients
hippo_df = calculate_vertical_gradient(hippo_df)
# Detect intrusion events based on ozone vertical gradient
hippo_df = detect_intrusions(hippo_df, threshold)

# --- #
# Define criteria for filtering the dataframe
# --- #
# Define filtering criteria in dict
filter_criteria = {
    'p': (300, 600),
    # 'gradient_O3': (-300, None),  # Using None for open-ended range
    # 'H2O': (None, 100),
     'lat': (-60, 60)
}

# keep original data as copy
unfiltered_df = hippo_df.copy() 
# Filter intrusion data
hippo_df = filter_dataframe(hippo_df, filter_criteria)
# colocate data 
# and sort into intrusion and non-intrusion dataframes
intrusion_data, non_intrusion_data = filter_colocate_data(hippo_df)

# --- #
# Data display and plot figures
# --- #
tracer_1 = 'CH4'
tracer_2 = 'O3'
# plot correlations in scatter plot
plot_correlation_comparison(intrusion_data, non_intrusion_data, tracer_1, tracer_2)
# Plot correlations in heatmap
plot_correlations(intrusion_data, non_intrusion_data, tracers)
