# agrigate the Strat N2O ECMWF daily reanalysis  files to yearly to save data 

import os
import xarray as xr
import glob

# Define the file directory and file pattern
data_dir = os.getcwd()

file_pattern = 'cams73_latest_n2o_conc_surface_inst_1996*.nc'

def process_and_save_one_year(year):
    # Create the file pattern for the specific year
    year_file_pattern = file_pattern.replace('1996', str(year))
    print(year_file_pattern)
    matching_files = glob.glob(year_file_pattern)
    print(matching_files)
    # Load the dataset for the year using dask
    ds = xr.open_mfdataset(os.path.join(data_dir, year_file_pattern), engine='netcdf4', chunks='auto')

    # Take the mean over the time dimension
    mean_ds = ds.mean(dim='time')

    # Save the yearly mean data to a new NetCDF file
    output_filename = 'combined/mean_n2o_conc_surface_{year}.nc'
    mean_ds.to_netcdf(os.path.join(data_dir, output_filename))

def main():
    start_year = 1984
    end_year = 2022  # Adjust this to your desired end year

    for year in range(start_year, end_year + 1):
        process_and_save_one_year(year)

if __name__ == "__main__":
    main()

