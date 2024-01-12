### Newton Nguyen 
### Script for plotting tracor correlations 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from datetime import datetime, timedelta


def matlab2datetime(timestamps):
    """
    Converts MATLAB-style timestamps from days to DateTime objects.
    
    Args:
    timestamps (Series): Pandas Series of timestamps in days.
    
    Returns:
    Series: Converted timestamps as DateTime objects.
    """
    matlab_base_date = pd.Timestamp('0000-12-31')  # MATLAB's default base date
    return timestamps.apply(lambda x: matlab_base_date + timedelta(days=x))


def read_csv_to_dataframe(file_path, time_column='timestamp'):
    """
    Reads a CSV file into a pandas DataFrame 
    and converts a specified timestamp column from Matlab Datenum to Python DateTime.
    
    Args:
    file_path (str): Path to the CSV file.
    time_column (str): Name of the column containing the timestamp in nanoseconds.
    
    Returns:
    DataFrame: The CSV data as a pandas DataFrame with DateTime conversion.
    """
    df = pd.read_csv(file_path)
    df[time_column] = pd.to_datetime((df[time_column]))
    return df


def annual_average(df):
    """
    Computes the annual average of the given DataFrame, ignoring NaN values.

    Args:
    df (DataFrame): The pandas DataFrame containing the data.

    Returns:
    DataFrame: Annual average data as a pandas DataFrame.
    """
    # Ensure 'timestamp' is in datetime format and set as index
    #df['timestamp'] = pd.to_datetime(df['timestamp'])
    df.set_index('timestamp', inplace=True)

    # Take annual average
    annual_avg_df = df.resample('Y').mean()

    # Reset the index and make timestamp behave like a normal pd column 
    annual_avg_df.reset_index(inplace=True)

    return annual_avg_df




def calculate_growth_rate(df, column_name):
    """
    Calculates the growth rate (difference) of a specified concentration column,
    placing NaN at the beginning of the resulting difference vector.
    
    Args:
    df (DataFrame): The pandas DataFrame containing the data.
    column_name (str): The name of the column to calculate the growth rate for.
    
    Returns:
    ndarray: Growth rate of the specified column with NaN at the beginning.
    """
    differences = np.diff(df[column_name])
    return np.insert(differences, 0, np.nan)


def plot_tracor_growth_correlation(df, strat_column, trop_column, save_fig=False, fig_dir='figures'):
    """
    Plots the relationship between strat concentration and the growth rate of tropospheric concentration.
    
    Args:
    df (DataFrame): The pandas DataFrame containing the data.
    strat_column (str): The name of the strat concentration column.
    trop_column (str): The name of the tropospheric concentration column.
    """
    df['growth_rate'] = calculate_growth_rate(df, trop_column)
    sns.scatterplot(data=df, x=strat_column, y='growth_rate')
    plt.xlabel(strat_column)
    plt.ylabel(trop_column + ' growth-rate')
    plt.title(f'Correlation of {strat_column} and {trop_column} growth-rate')
    if save_fig:
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        fig_filename = f"{fig_dir}/{strat_column}_vs_{trop_column}_growth_rate.png"
        plt.savefig(fig_filename, dpi=300)
        print(f"Figure saved as {fig_filename}")

    

def print_correlation(df, strat_column, trop_column):
    """
    Prints the correlation coefficient between strat concentration and the growth rate of tropospheric concentration,
    ignoring NaN values.
    
    Args:
    df (DataFrame): The pandas DataFrame containing the data.
    strat_column (str): The name of the strat concentration column.
    trop_column (str): The name of the tropospheric concentration column.
    """
    df['growth_rate'] = calculate_growth_rate(df, trop_column)
    correlation = df[strat_column].corr(df['growth_rate'])
    print(f'Correlation between {strat_column} and Growth Rate of {trop_column}: {correlation}')



### main script here 
df = read_csv_to_dataframe('monthly_obs.csv')
df = annual_average(df)
#df = df.iloc[8:]
# Define pairs for plotting the correlation 
# the tupple plots (strat_obs, trop_tracor_growth_rate)
h2o_strat = 'h2o_global_strat'
#h2o_strat = 'h2o_tropical_strat'
pairs = [
    (h2o_strat, 'nh_ch4'), (h2o_strat, 'sh_ch4'),
    (h2o_strat, 'nh_ch4c13'), (h2o_strat, 'sh_ch4c13'),
    ('nh_n2o_strat', 'nh_n2o'), ('sh_n2o_strat', 'sh_n2o'),
    (h2o_strat, 'nh_n2o_strat'), (h2o_strat, 'sh_n2o_strat')
]

# Plot and print correlations
for strat_column, trop_column in pairs:
    print_correlation(df, strat_column, trop_column)
    plot_tracor_growth_correlation(df, strat_column, trop_column)

