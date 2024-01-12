function out = read_strat_obs(filename, St, tRes)
    % Open the NetCDF file
    variableNames = {'combinedanomfillanomh2oq', 'combinedeqfillh2oq', 'time'};
    ncid = netcdf.open(filename, 'NC_NOWRITE');
    % Initialize the data structure
    out = struct();
    fields = {'anomaly', 'obs', 'time'};

    % Read the variable data
    for i = 1:numel(variableNames)
        varid = netcdf.inqVarID(ncid, variableNames{i});
        out.(fields{i}) = netcdf.getVar(ncid, varid);
    end % for loop 

    % Close the NetCDF file
    netcdf.close(ncid);

    % Convert time from days since Jan 1, 1984 to a MATLAB datetime
    out.time = datetime(out.time + datenum(1984,1,1), 'ConvertFrom', 'datenum');
    St = datetime(St, 'ConvertFrom', 'datenum');

    % Select data within the specified timespan
    t_start = find(out.time >= St(1), 1);
    t_stop = find(out.time <= St(end)+16, 1, 'last')
 %  to include end of period

    % Filter the fields for the selected time span
    for i=1:length(fields)
        field = fields{i};
        out.(field) = out.(field)(t_start:t_stop);
        % Check if annual average is needed
        if strcmp(tRes, 'year') && ~strcmp(field, 'time')
            [out.(field), time] = computeAnnualAverage(out.(field), out.time(t_start:t_stop));
        end % conditional 
    end % for loop 
%    out.time = time;
end % function 

function [annualAverage, uniqueYears] = computeAnnualAverage(monthlyData, timeVector)
    % Convert time to years and months
    dates = datevec(timeVector);
    years = dates(:,1);

    % Initialize the annual average vector
    uniqueYears = unique(years);
    annualAverage = zeros(length(uniqueYears), 1);

    % Loop through each year and compute the average
    for i = 1:length(uniqueYears)
        year = uniqueYears(i);
        % Find indices for the current year
        indices = years == year;
        % Calculate the average for the current year
        annualAverage(i) = mean(monthlyData(indices));
    end % for loop
end % function 
