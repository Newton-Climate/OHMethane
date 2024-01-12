import os
import xarray as xr
import numpy as np

def calculate_observation_uncertainty(rms, std):
    # Square the RMS and standard deviation values
    rms_squared = rms ** 2
    std_squared = std ** 2

    # Calculate the observation uncertainty using the RSS method
    uncertainty = np.sqrt(rms_squared + std_squared)

    return uncertainty

def calculate_hemisphere_uncertainty(dataset, hemisphere_lat_range):
    # Average over specified pressure levels
    dataset = dataset.sel(lev=slice(100, 45)).mean(dim="lev", skipna=True)

    # Extract the RMS precision and standard deviation values for the specified hemisphere
    rms_values = dataset['rms_uncertainty'].sel(lat=hemisphere_lat_range).mean(dim="lat", skipna=True)
    std_values = dataset['std_dev'].sel(lat=hemisphere_lat_range).mean(dim="lat", skipna=True)

    # Calculate the observation uncertainty for the hemisphere dataset
    uncertainty = calculate_observation_uncertainty(rms_values, std_values)

    return uncertainty

def process_netCDF_files(directory="."):
    # Generate a list of netCDF files in the directory
    files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".nc")]

    # Load the dataset from the netCDF files
    dataset = xr.open_mfdataset(files, group="N2O PressureZM")

    # Compute the annual time average
    dataset = dataset.resample(time="MS").mean(skipna=True)

    # Split the dataset into northern and southern hemisphere
    northern_hemisphere = dataset.sel(lat=slice(0, 23.5), lev=slice(200, 45.0))
    southern_hemisphere = dataset.sel(lat=slice(-23.5, 0), lev=slice(200, 45))


    # Calculate the observation uncertainty for the northern and southern hemisphere datasets
    uncertainty_north = calculate_hemisphere_uncertainty(northern_hemisphere, slice(0, 23.5))
    uncertainty_south = calculate_hemisphere_uncertainty(southern_hemisphere, slice(-23.5, 0))

    # take average over each hemisphere (NH and SH)
    northern_hemisphere = northern_hemisphere.mean(dim=["lev", "lat"], skipna=True)
    print(uncertainty_north)
    southern_hemisphere = southern_hemisphere.mean(dim=["lev", "lat"], skipna=True)

    # Create a new dataset for each hemisphere, including the observations and uncertainties
    north_data = xr.Dataset({
        "observation": northern_hemisphere.value,
        "uncertainty": uncertainty_north
    }, coords={"time": northern_hemisphere.time})

    south_data = xr.Dataset({
        "observation": southern_hemisphere.value,
        "uncertainty": uncertainty_south
    }, coords={"time": southern_hemisphere.time})

    return north_data, south_data






def export_to_netcdf(nh_data, sh_data, output_file):
    # Create a new dataset for the output file
    output_dataset = xr.Dataset()

    # Add the variables and coordinates for the northern hemisphere
    for var_name in nh_data.variables:
        output_dataset[var_name + "_nh"] = nh_data[var_name]
    for coord_name in nh_data.coords:
        output_dataset[coord_name] = nh_data[coord_name]

    # Add the variables and coordinates for the southern hemisphere
    for var_name in sh_data.variables:
        output_dataset[var_name + "_sh"] = sh_data[var_name]
    for coord_name in sh_data.coords:
        output_dataset[coord_name] = sh_data[coord_name]

    # Export the output dataset to NetCDF file
    output_dataset.to_netcdf(output_file)




# Example usage
north_data, south_data = process_netCDF_files()
output_file = 'N2O_lowerstrat_obs.nc'  # The desired output NetCDF file path

export_to_netcdf(north_data, south_data, output_file)
print(north_data.observation.values)
