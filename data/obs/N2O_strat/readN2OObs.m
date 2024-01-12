function out = readN2OStrat(file_path)
    % Read the NetCDF file
    data = ncinfo(file_path);  % Get information about the NetCDF file
    
    % Initialize the data structure to hold the variables and dimensions
    data_vars = struct();
    data_dims = struct();
    
    % Access the variable names and dimensions
    var_names = {data.Variables.Name};  % Cell array of variable names
    dim_names = {data.Dimensions.Name};  % Cell array of dimension names
    
    % Access the data variables
    for i = 1:numel(var_names)
        var_name = var_names{i};
        var_data = ncread(file_path, var_name);
        
        % Store the variable data in the data_vars structure
        data_vars.(var_name) = var_data;
    end
    
    % Access the coordinate variables
    for i = 1:numel(dim_names)
        dim_name = dim_names{i};
        dim_data = ncread(file_path, dim_name);
        
        % Store the coordinate data in the data_dims structure
        data_dims.(dim_name) = dim_data;
    end
    
    % Assign the data_vars and data_dims to the 'data' output variable
    data.data_vars = data_vars;
    data.data_dims = data_dims;

    % output the final struct 
    out.time = [2004:2022];
    out.nh_n2o_strat = data_vars.observation_nh;
    out.nh_n2o_strat_err = data_vars.uncertainty_nh;
    out.sh_n2o_strat = data_vars.observation_sh;
    out.sh_n2o_strat_err = data_vars.uncertainty_sh;

end
