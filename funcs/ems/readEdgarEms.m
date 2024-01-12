function ems_struct = readEdgarEms(filename, sYear, eYear)
    % filename: name of npy file as STR 
    % sYear: start year for model run
    % eYear: end year of data for model run 

    %%% output 
    % struct: struct of Edgar inventory between sYear and eYear 

    % read data file 
    % dimensions are (sector, time, hemisphere (either nh or sh))
    ds = readNPY(filename);
    [~, n_time, ~] = size(ds);

    % time axis of data 
    t = 1970:1970+n_time-1;
    % find index of eYear and sYear in t
    t_start = find(t == sYear);
    t_stop = find(t == eYear);

% index for output array
    t_out = [sYear:eYear];
    out_start = find(t_out == sYear);
    out_stop = find(t_out == eYear);


    % Check if sYear or eYear is outside the range
    if sYear < t(1)
        t_start = 1;
        out_start = find(t_out == t(1));
    end
    if eYear > t(end)
        t_stop = length(t);
        out_stop = find(t_out == t(end));
    end

    sectors = {'ff', 'waste', 'rice', 'animal', 'biofuel'};

    ems_struct = struct();  % Initialize the output struct

    % Create arrays of NaN values
    num_years = eYear - sYear + 1;
    nan_data = NaN(num_years, 1);

    for i = 1:length(sectors)
        sector = sectors{i};  % Retrieve the sector name

        % get fieldname for ems_struct 
        nh_fieldname = [sector '_nh'];
        sh_fieldname = [sector '_sh'];

        % assign the correct ds slice and timespan to ems_struct 
        ems_struct.(nh_fieldname) = nan_data;
        ems_struct.(sh_fieldname) = nan_data;

        % Fill in the correct data
        if t_start <= t_stop
            ems_struct.(nh_fieldname)(out_start:out_stop) = ds(i, t_start:t_stop, 2);
            ems_struct.(sh_fieldname)(out_start:out_stop) = ds(i, t_start:t_stop, 1);
        end
    end

    % add up biofuels and fossil fuels 
    ems_struct.energy_nh = ems_struct.ff_nh + ems_struct.biofuel_nh;
    ems_struct.energy_sh = ems_struct.ff_sh + ems_struct.biofuel_sh;
    ems_struct = rmfield(ems_struct, 'ff_nh');
    ems_struct = rmfield(ems_struct, 'ff_sh');
    ems_struct = rmfield(ems_struct, 'biofuel_nh');
    ems_struct = rmfield(ems_struct, 'biofuel_sh');

end