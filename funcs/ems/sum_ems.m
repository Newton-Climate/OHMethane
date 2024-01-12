function [nh_emissions_total, sh_emissions_total] = sumHemisphericEms(ems)
    %%% = calculate the total emissions per hemisphere (nh and sh)
    %%% = from bottom up inventories denoted by a struct
    %%%% = input:
    %%%% = ems: a struct that has fieldnames for each sector .e.g, ems.nh_wetland, ems.sh_wetland
    %%%% =  output: a vector timeseries of emissions for each hemisphere 
    %%%% = 

    nh_emissions_total = 0;
    sh_emissions_total = 0;
    
    % Check sectors in NH
    nh_sectors = fieldnames(ems);
    for i = 1:numel(nh_sectors)
        sector = nh_sectors{i};
        if endsWith(sector, '_nh')
            nh_emissions_total = nh_emissions_total + ems.(sector);
        end
    end
    
    % Check sectors in SH
    sh_sectors = fieldnames(ems);
    for i = 1:numel(sh_sectors)
        sector = sh_sectors{i};
        if endsWith(sector, '_sh')
            sh_emissions_total = sh_emissions_total + ems.(sector);
        end
    end
end
