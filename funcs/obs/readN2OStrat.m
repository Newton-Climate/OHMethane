function obs = readN2OStrat(file_path, St, obs)
    % Read the NetCDF file
    data = ncinfo(file_path);  % Get information about the NetCDF file
    
    % Initialize the data structure to hold the variables
    data_vars = struct();
    
    % Access the variable names
    var_names = {data.Variables.Name};  % Cell array of variable names
    
    % Access the data variables
    for i = 1:numel(var_names)
        var_name = var_names{i};
        var_data = ncread(file_path, var_name);
        
        % Store the variable data in the data_vars structure
        data_vars.(var_name) = var_data;
    end
    
    % Extract time from the netCDF file
    n2o_time = datenum(data_vars.time) + datenum(2004,1,1);
    
    % Find the indices of the time range in the netCDF file
    n2o_ind1 = find(n2o_time >= St(1), 1, 'first');
    n2o_ind2 = find(n2o_time <= St(end), 1, 'last');

    % Prepare a series for NaN values
    nT = length(St);
    series = NaN * ones(nT, 1);
    
    % Initialize the indices for obs struct
    obs_ind1 = 1;
    obs_ind2 = nT;
    
    % Adjust the indices based on the available time range in n2o_time
    if n2o_time(n2o_ind1) > St(1)
        obs_ind1 = find(St == n2o_time(n2o_ind1), 1, 'first');
    end
    if n2o_time(n2o_ind2) < St(end)
        obs_ind2 = find(St == n2o_time(n2o_ind2), 1, 'last');
    end

    % Assign the data to the obs struct
    obs.nh_n2o_strat = series;
n2o_ind1
n2o_ind2
obs_ind1
obs_ind2
    obs.nh_n2o_strat(obs_ind1:obs_ind2) = data_vars.observation_nh(n2o_ind1:n2o_ind2) * 1e9;
    obs.nh_n2o_strat_err = series;
    obs.nh_n2o_strat_err(obs_ind1:obs_ind2) = data_vars.uncertainty_nh(n2o_ind1:n2o_ind2) * 1e9;
    obs.sh_n2o_strat = series;
    obs.sh_n2o_strat(obs_ind1:obs_ind2) = data_vars.observation_sh(n2o_ind1:n2o_ind2) * 1e9;
    obs.sh_n2o_strat_err = series;
    obs.sh_n2o_strat_err(obs_ind1:obs_ind2) = data_vars.uncertainty_sh(n2o_ind1:n2o_ind2) * 1e9;
end


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