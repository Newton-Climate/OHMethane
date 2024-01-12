function ems_struct = readFireEms(file_path, sYear, eYear)

% Read the file
file = fopen(file_path);
lines = textscan(file, '%s', 'Delimiter', '\n');
lines = lines{1};
fclose(file);

% Initialize structures
nh_ems = struct();
sh_ems = struct();
years = [];

% Process the lines
for i = 1:length(lines)
    line = strtrim(lines{i});
    
    % Skip if the line is a comment, empty, only whitespace, or ends with ':'
    if isempty(line) || line(1) == '#' || all(isspace(line)) || line(end) == ':'
        continue
    end

    % Split line into parts
    parts = strsplit(line, ' ');


    if startsWith(parts{1}, 'Region') 
        % Get years
        years = [years cellfun(@str2num, parts(2:end-2))];
    else
        % Get emissions
        ems = cellfun(@str2num, parts(2:end-2));




        % Divide into NH and SH
        if ismember(parts{1}, {'BONA', 'TENA', 'CEAM', 'NHSA', 'EURO', 'MIDE', 'NHAF', 'BOAS', 'TEAS'})
            nh_ems.(parts{1}) = ems/1e2;
        else
            sh_ems.(parts{1}) = ems/1e2;
        end
    end

end

% Sum NH and SH emissions
nh_ems.Total = zeros(size(years));
sh_ems.Total = zeros(size(years));

fields = fieldnames(nh_ems);
for i = 1:numel(fields)
    nh_ems.Total = nh_ems.Total + nh_ems.(fields{i});
end

fields = fieldnames(sh_ems);
for i = 1:numel(fields)
    sh_ems.Total = sh_ems.Total + sh_ems.(fields{i});
end

    %%% Now, output the final data struct
    % time axis of data 
  t_data = years;
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
    ems_struct.nh_fires = nan_data;
    ems_struct.sh_fires = nan_data;

    % Fill in the correct data
        if data_start <= data_stop
    ems_struct.nh_fires(out_start:out_stop) = nh_ems.Total(data_start:data_stop);
    ems_struct.sh_fires(out_start:out_stop) = sh_ems.Total(data_start:data_stop);
end

end

