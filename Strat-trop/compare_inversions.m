load('model_outputs/inversion_noStrat.mat')
ems_nostrat = ems_anal;
ch4_nostrat = ems_anal(:,1) + ems_anal(:,2);
oh_nostrat = (out_anal.nh_oh + out_anal.sh_oh) / 2; 
model_nostrat_MCF = out_anal;
obs_nostrat_MCF = obs;
oh_nostrat = 100 * (oh_nostrat - mean(oh_nostrat)) ./ mean(oh_nostrat);
ems_nostrat = makeEmsStruct(ems_nostrat);


load('model_outputs/inversion_noStratMCF.mat')
ems_nostrat_MCF = ems_anal;
ch4_nostrat_MCF = ems_anal(:,1) + ems_anal(:,2);
oh_nostrat_MCF = (out_anal.nh_oh + out_anal.sh_oh) / 2;
model_nostrat = out_anal
obs_nostrat = obs;
oh_nostrat_MCF = 100 * (oh_nostrat_MCF - mean(oh_nostrat_MCF)) ./ mean(oh_nostrat_MCF);
ems_nostrat_MCF = makeEmsStruct(ems_nostrat_MCF);

load('model_outputs/inversion_H2O_constraint.mat')
ems_strat = ems_anal;
ch4_strat = ems_anal(:,1) + ems_anal(:,2);
model_strat = out_anal;
obs_strat = obs;
oh_strat = (out_anal.nh_oh + out_anal.sh_oh) / 2;
oh_strat = 100 * (oh_strat - mean(oh_strat)) ./ mean(oh_strat);
ems_strat = makeEmsStruct(ems_strat);


export_to_csv()

% Create a tiled chart layout

% Your data setup
oh = [oh_strat, oh_nostrat_MCF, oh_nostrat];
ch4 = [ch4_strat, ch4_nostrat_MCF, ch4_nostrat];
t = [sYear:eYear];

figure(1);
tlo = tiledlayout(2, 1, "TileSpacing", "compact");

% Your data setup...
% Define plot aesthetics...

% Plotting CH4 emissions
nexttile;
hold on; % This ensures all plots are held on the same subplot
for i = 1:3
    plot(t, ch4(:,i), 'Color', colors{i}, 'LineStyle', lineStyles{i}, ...
         'LineWidth', lineWidth);
end
ylabel('CH4 emissions (Tg/yr)');
set(gca, 'FontSize', fontSize);
xlim([1992, 2022]);

% Plotting OH variability
nexttile;
hold on; % This ensures all plots are held on the same subplot
for i = 1:3
    plot(t, oh(:,i), 'Color', colors{i}, 'LineStyle', lineStyles{i}, ...
         'LineWidth', lineWidth);
end
xlabel('Year');
ylabel('[OH] variability (%)');
set(gca, 'FontSize', fontSize);
xlim([1992, 2022]);

% Adding central title
sgtitle(tlo, 'CH4 and [OH] Inversions ', 'FontSize', 14);

% Adding a shared legend for the entire layout
lgd = legend({'No strat', 'CH4 strat exchange / no MCF exchange', 'CH4 and MCF exchange'}, ...
             'Location', 'northoutside', 'Orientation', 'horizontal');
lgd.Layout.Tile = 'north'; % Places legend at the top, outside the subplots

% Save the figure
saveas(figure(1), 'model_outputs/ch4_oh_variability.pdf', 'pdf');


% Setup for the figure
% Create a tiled chart layout
figure(2);
tlo = tiledlayout(4, 1, "TileSpacing","compact"); 

% Assuming models and obs are already loaded with the required fields
models = [model_nostrat, model_nostrat_MCF, model_strat];

% Loop for plotting NH and SH fits

    % NH CH4 Fit Plot
    nexttile; % Automatically creates the next subplot in the layout
    hold on;
for i = 1:numel(models)
    plot(t, 100 * (obs.nh_ch4 - models(i).nh_ch4) ./ obs.nh_ch4, ...
         'Color', colors{i}, 'LineStyle', lineStyles{i}, ...
         'LineWidth', lineWidth);
    ylabel('NH CH4 fits (%)');
    xlim([1992, 2022]);
    set(gca, 'FontSize', fontSize);

end

    % NH MCF Fit Plot
    nexttile;
    hold on;
for i = 1:numel(models)

    plot(t, 100 * (obs.nh_mcf - models(i).nh_mcf) ./ obs.nh_mcf, ...
         'Color', colors{i}, 'LineStyle', lineStyles{i}, ...
         'LineWidth', lineWidth);
    ylabel('NH MCF fits (%)');
    xlim([1992, 2022]);
    set(gca, 'FontSize', fontSize);

end


    % SH CH4 Fit Plot

    nexttile;
    hold on;
for i = 1:numel(models)
    plot(t, 100 * (obs.sh_ch4 - models(i).sh_ch4) ./ obs.sh_ch4, ...
         'Color', colors{i}, 'LineStyle', lineStyles{i}, ...
         'LineWidth', lineWidth);
    ylabel('SH CH4 fits (%)');
    xlim([1992, 2022]);
    set(gca, 'FontSize', fontSize);

end

    % SH MCF Fit Plot
    nexttile;
    hold on;
for i = 1:numel(models)
    plot(t, 100 * (obs.sh_mcf - models(i).sh_mcf) ./ obs.sh_mcf, ...
         'Color', colors{i}, 'LineStyle', lineStyles{i}, ...
         'LineWidth', lineWidth);
    xlabel('Year');
    ylabel('SH MCF fits (%)');
    xlim([1992, 2022]);
    set(gca, 'FontSize', fontSize);

end

% Adding central title and grid
sgtitle(tlo, 'Box model fits', 'FontSize', 14);
tlo.TileSpacing = 'compact'; % Reduces spacing between plots
tlo.Padding = 'compact'; % Reduces padding around plots

% Adding Legend
% Adding Legend to the last subplot
lgd = legend({'No strat', 'CH4 strat exchange / no MCF exchange', 'CH4 and MCF exchange'}, ...
             'Location', 'northoutside', 'Orientation', 'horizontal');
lgd.Layout.Tile = 'north'; % Places legend at the top, outside the subplots

% Save the figure
saveas(figure(2), 'model_outputs/inversion_results_fits.pdf', 'pdf');



function ems_struct = makeEmsStruct(ems_array)

  % - NH CH4 emissions
% - SH CH4 emissions
% - NH CH4C13 composition
% - SH CH4C13 composition
% - NH MCF emissions
% - SH MCF emissions
% - NH N2O emissions
% - SH N2O emissions
% - NH C2H6 emissions
% - SH C2H6 emissions
% - NH OH emissions
% - SH OH emissions
% - NH CO emissions
% - SH CO emissions
% - Strat-trop exchange
% - NH arbitrary OH reaction rate
% - SH arbitrary OH reaction rate

ems_struct = struct;

ems_struct.nh_ch4 = ems_array(:,1);
ems_struct.sh_ch4 = ems_array(:,2);
ems_struct.nh_ch4c13 = ems_array(:,3);
ems_struct.sh_ch4c13 = ems_array(:,4);
ems_struct.nh_mcf = ems_array(:,5);
ems_struct.sh_mcf = ems_array(:,6);
ems_struct.nh_n2o = ems_array(:,7);
ems_struct.sh_n2o = ems_array(:,8);
ems_struct.nh_c2h6 = ems_array(:,9);
ems_struct.sh_c2h6 = ems_array(:,10);
ems_struct.nh_oh = ems_array(:,11);
ems_struct.sh_oh = ems_array(:,12);
ems_struct.nh_co = ems_array(:,13);
ems_struct.sh_co = ems_array(:,14);
ems_struct.tau = ems_array(:,15);
ems_struct.kx_nh = ems_array(:,16);
ems_struct.kx_sh = ems_array(:,17);
end


function export_to_csv()
    % Load your data
    load('model_outputs/inversion_noStrat.mat');
    create_csv(ems_anal, out_anal, obs, 'model_outputs/no_strat.csv', sYear, eYear);

    load('model_outputs/inversion_noStratMCF.mat');
    create_csv(ems_anal, out_anal, obs, 'model_outputs/no_strat_MCF.csv', sYear, eYear);

    load('model_outputs/inversion_H2O_constraint.mat');
    create_csv(ems_anal, out_anal, obs, 'model_outputs/strat.csv', sYear, eYear);
end

function create_csv(ems_anal, model_struct, obs_struct, filename, sYear, eYear)
    % Create EMS structure
    ems_struct = makeEmsStruct(ems_anal);
    
    % Prepare data table
    timestamp = (sYear:eYear)';
    data_table = table(timestamp);
    
    % Append EMS data
    ems_fields = fieldnames(ems_struct);
    for i = 1:numel(ems_fields)
        data_table.(strcat(ems_fields{i}, '_ems')) = ems_struct.(ems_fields{i});
    end
    
    % Append MODEL and OBS data (similarly as for EMS)
    % Append EMS data
    obs_fields = fieldnames(obs_struct);
    for i = 1:numel(obs_fields)
        data_table.(strcat(obs_fields{i}, '_obs')) = obs_struct.(obs_fields{i});
    end

    % Append EMS data
    model_fields = fieldnames(model_struct);
    for i = 1:numel(model_fields)
        data_table.(strcat(model_fields{i}, '_model')) = model_struct.(model_fields{i});
    end


    % ...

    % Export to CSV
    writetable(data_table, filename);
end
