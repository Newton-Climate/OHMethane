import os
import xarray as xr
import numpy as np

def process_netCDF_files(directory="."):
    # Generate a list of netCDF files in the directory
    files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".nc")]

    # Load the dataset from the netCDF files
    dataset = xr.open_mfdataset(files, group="N2O PressureZM")

    # Compute the annual time average while skipping NaN values
    dataset = dataset.resample(time="AS").mean(skipna=True)

    # Split the dataset into northern and southern hemisphere
    tropics = dataset.sel(lat=slice(-30, 30), lev=slice(100, 45))
#    southern_hemisphere = dataset.sel(lat=slice(-90, 0), lev=slice(100, 45))

    # Spatially average over the hemispheres while skipping NaN values
    tropics_avg = tropics.mean(dim=["lat", "lev"], skipna=True)
#    sh_avg = southern_hemisphere.mean(dim=["lat", "lev"], skipna=True)

    # Concatenate the results along the time axis
    #result = xr.concat([nh_avg, sh_avg], dim="time")

    return tropics_avg

tropics = process_netCDF_files()