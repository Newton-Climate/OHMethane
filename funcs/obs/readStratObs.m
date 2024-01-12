function obs = readStratObs(filename, St, tRes)
    % Open the NetCDF file
    ncid = netcdf.open(filename, 'NC_NOWRITE');

    
    % Get information about the file
    [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);

    % Initialize the data structure
    obs = struct();

    % Read the variable data
    for i = 0:numvars-1 % Adjust for 0-based indexing in netCDF API
        [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid, i);
        obs.(varname) = netcdf.getVar(ncid, i);
    end

    % Close the NetCDF file
    netcdf.close(ncid);

    % Convert time from Unix epoch to MATLAB datetime
    % Check if 'time' is a field in the struct
    if isfield(obs, 'time')
        referenceDate = datenum('1970-01-01'); % Change this to the actual reference date;
        obs.time = datetime(referenceDate + double(obs.time), 'ConvertFrom', 'datenum');    
    end

    St = datetime(St, 'ConvertFrom', 'datenum'); % Convert St to datetime

    % Select data within the specified timespan
    t_start = find(obs.time >= St(1), 1);
    t_stop = find(obs.time <= St(end), 1, 'last');

    % Filter the fields for the selected time span
    fieldNames = fieldnames(obs);
    for i = 1:numel(fieldNames)
        field = fieldNames{i};
        obs.(field) = obs.(field)(t_start:t_stop);
        % Check if annual average is needed
        if strcmp(tRes, 'year') && ~strcmp(field, 'time')
            [obs.(field), years] = computeAnnualAverage(obs.(field), obs.time(t_start:t_stop));
        end
    end
    obs.time = years
end

function [annualAverage, uniqueYears] = computeAnnualAverage(monthlyData, timeVector)
    % ... [Function remains the same as previously defined] ...
end
