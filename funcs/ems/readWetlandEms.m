function ems_struct = readWetlandEms(file_path, sYear, eYear)

    % Read the variable and coordinates from the NetCDF file
    nh_wetlands = ncread(file_path, 'nh_wetland_ems');
    sh_wetlands = ncread(file_path, 'sh_wetland_ems');
    t_data = ncread(file_path, 'time');

    % time axis of data 
  t_data = [t_data/365 + 2001] ;
    n_data = length(t_data);

    % find index of eYear and sYear in t_data
    data_start = find(t_data == sYear);
    data_stop = find(t_data == eYear);

% index for output array
    t_out = [sYear:eYear];
    out_start = find(t_out == sYear);
    out_stop = find(t_out == eYear);


    % Check if sYear or eYear is outside the range
    if sYear < t_data(1)
        data_start = 1;
        out_start = find(t_out == t_data(1));
    end
    if eYear > t_data(end)
        data_stop = length(t_data);
        out_stop = find(t_out == t_data(end));
    end

    ems_struct = struct();  % Initialize the output struct

    % Create arrays of NaN values
    num_years = eYear - sYear + 1;
    nan_data = NaN(num_years, 1);

        % assign the correct ds slice and timespan to ems_struct 
    ems_struct.nh_wetlands = nan_data;
    ems_struct.sh_wetlands = nan_data;

    % Fill in the correct data
        if data_start <= data_stop
    ems_struct.nh_wetlands(out_start:out_stop) = nh_wetlands(data_start:data_stop);
    ems_struct.sh_wetlands(out_start:out_stop) = sh_wetlands(data_start:data_stop);
end