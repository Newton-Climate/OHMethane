function [nh_total, sh_total] = hemisphericEms(ems)
    % initialize total emissions as zeros
    nh_total = 0;
    sh_total = 0;
    
    % get all the field names of the structure
    fields = fieldnames(ems);
    
    for i = 1:length(fields)
        field_name = fields{i};
        % check if field name ends with '_nh' or '_sh'
        if endsWith(field_name, '_nh')
            nh_total = nh_total + ems.(field_name);
        elseif endsWith(field_name, '_sh')
            sh_total = sh_total + ems.(field_name);
        end
    end
end
