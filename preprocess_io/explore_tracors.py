import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


PLOT_FLAG = True  # Set to True to see plots, False to skip plotting
MAX_LAGS = 10  # Maximum number of lags for cross-correlation

def plot_ccf(x, y, ax):
    # Align series by index
    x, y = x.align(y, join='inner')
    
    # Drop NaN values after aligning
    x, y = x.dropna(), y.dropna()

    # Check if either series is empty
    if len(x) == 0 or len(y) == 0:
        ax.set(title=f'Insufficient data for {y.name}', xlabel='Lag', ylabel='Correlation Coefficient')
        return

    lags = np.min([MAX_LAGS, len(x) - 1, len(y)-1])
    ccf_values = ax.xcorr(x, y, usevlines=True, maxlags=lags, normed=True, lw=2)
    ax.set(title=f'Cross-correlation with {y.name}', xlabel='Lag', ylabel='Correlation Coefficient')

    # Identify and print lag with maximum correlation
#    print(ccf_values[0])
#    print(ccf_values[1])
    maxlag = ccf_values[0][np.argmax(np.abs(ccf_values[1]))]
    maxcorrelation = ccf_values[1][np.argmax(np.abs(ccf_values[1]))]
    print(f"Max cross-correlation with {y.name} is {maxcorrelation:.2f} at lag {maxlag}")


data = pd.read_csv('tracers.csv', parse_dates=True, index_col='year')
#data.enso_index = np.abs(data.enso_index)

time_periods = {
#    "1990-2000": ("1990-01-01", "2000-12-31"),
    "2000-2007": ("2000-01-01", "2007-12-31"),
    "2007-2019": ("2007-01-01", "2019-12-31"),
    "2019-2022": ("2019-01-01", "2022-12-31"),
    "2000-2021": ("2000-01-01", "2021-12-31")
}

for period, (start_date, end_date) in time_periods.items():
    print(f"\nAnalysis for {period}\n" + "="*50)
    period_data = data.loc[start_date:end_date]

    # Calculating simple correlation
    corr = period_data.corr()
    print(f"Correlation with H2O_strat_ppm:\n{corr['H2O_strat_ppm']}\n")

    # Plotting cross-correlation if PLOT_FLAG is set to True
    if PLOT_FLAG:
        columns_to_analyze = [col for col in period_data.columns if col != 'H2O_strat_ppm']
        fig, axes = plt.subplots(nrows=len(columns_to_analyze), figsize=(10, 6*len(columns_to_analyze)))

        # If only one column, axes will not be an array. Convert it into an array for consistency
        if len(columns_to_analyze) == 1:
            axes = [axes]

        for idx, column in enumerate(columns_to_analyze):
            plot_ccf(period_data["H2O_strat_ppm"].dropna(), period_data[column].dropna(), axes[idx])

        plt.tight_layout()
#        plt.show()
