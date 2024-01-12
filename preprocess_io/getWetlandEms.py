import xarray as xr
import numpy as np

def load_data(directory):
    files = f"{directory}/*.nc"
    dataset = xr.open_mfdataset(files)
    return dataset

def calculate_annual_average(dataset, variable_name):
    annual_average = dataset[variable_name].resample(time='AS').mean(dim='time')
    return annual_average

def separate_hemispheres(annual_average):
    northern_hemisphere = annual_average.sel(lat=slice(0, 90))
    southern_hemisphere = annual_average.sel(lat=slice(-90, 0))
    return northern_hemisphere, southern_hemisphere

def convert_units(ch4_data, latitudes):
    # Conversion factor: mg/m^2/day to Tg/yr
    conversion_factor = 1e-15 * 365.25

    # Calculate grid cell areas using latitude
    area_m2 = calculate_grid_cell_area(latitudes)

    # Apply area weighting and conversion factor
    ch4_data_converted = ch4_data * conversion_factor * area_m2

    return ch4_data_converted

def calculate_grid_cell_area(latitude):
    # Earth's semi-major axis in meters
    semi_major_axis = 6378137.0

    # Approximate length of one degree at the equator in meters
    degree_length_at_equator = 111319.9

    # Convert latitude to radians
    lat_rad = np.deg2rad(latitude)

    # Calculate the width and height of the grid cell in degrees
    height_deg = 0.5
    width_deg = 0.5 * np.cos(lat_rad)

    # Convert the width and height from degrees to meters
    height_m = height_deg * degree_length_at_equator
    width_m = width_deg * degree_length_at_equator * np.cos(lat_rad)

    # Calculate the area of the grid cell in square meters
    area_m2 = height_m * width_m

    return area_m2

# Example usage
data_directory = '/Users/newtonnguyen/Downloads/MonthlyWetland_CH4_WetCHARTs_1915/data/'

# Load data
dataset = load_data(data_directory)

# Calculate annual average
variable_name = 'wetland_CH4_emissions'
annual_average = calculate_annual_average(dataset, variable_name)

# Separate hemispheres
northern_hemisphere, southern_hemisphere = separate_hemispheres(annual_average)

# Get latitude values
latitudes = annual_average['lat']

# Convert units with area weighting and take mean along the 'model' dimension
northern_hemisphere_converted = convert_units(northern_hemisphere, latitudes).mean(dim='model')
southern_hemisphere_converted = convert_units(southern_hemisphere, latitudes).mean(dim='model')

# Calculate total emissions in each hemisphere
total_emissions_northern = northern_hemisphere_converted.sum(dim=['lat', 'lon'])
total_emissions_southern = southern_hemisphere_converted.sum(dim=['lat', 'lon'])

print("Total emissions (Northern Hemisphere):", total_emissions_northern.values)
print("Total emissions (Southern Hemisphere):", total_emissions_southern.values)


# Save results to a new NetCDF file
output_file = 'wetland_ems.nc'

# Create a new dataset for the hemispheric emissions
output_dataset = xr.Dataset(
    {
        'nh_wetland_ems': total_emissions_northern,
        'sh_wetland_ems': total_emissions_southern
    },
    coords={
        'time': annual_average.time,
    }
)

# Save the dataset to a new NetCDF file
output_dataset.to_netcdf(output_file)

print("Results saved to:", output_file)
