H2O_obs = read_H2O('global_strat_obs_monthly.nc')
H2O_strat = read_strat_obs('global_strat_obs.nc', sYear, eYear);
H2O_obs = read_H2O('tropic_strat_obs_yearly.nc')

% Extracting data based on latitude bands
%tropical_idx = find(latitudes >= -15 & latitudes <= 15);
%northern_polar_idx = find(latitudes >= 50 & latitudes <= 90);
%southern_polar_idx = find(latitudes >= -90 & latitudes <= -50);

% Initialize matrices for storing gradients
gradient_northern_time = zeros(size(H2O_obs.tim));
gradient_southern_time = zeros(size(H2O_obs.tim));

% Loop through each time point
for t = 1:length(H2O_obs.tim)
    mean_tropical = mean(arrayfun(@(i) ch4c13_obs.obs.(stations{i})(t), tropical_idx), 'omitnan');
    mean_northern_polar = mean(arrayfun(@(i) ch4c13_obs.obs.(stations{i})(t), northern_polar_idx), 'omitnan');
    mean_southern_polar = mean(arrayfun(@(i) ch4c13_obs.obs.(stations{i})(t), southern_polar_idx), 'omitnan');
    
    gradient_northern_time(t) = mean_northern_polar - mean_tropical;
    gradient_southern_time(t) = mean_southern_polar - mean_tropical;
end

% Plot the gradients over time
figure;
plot(H2O_obs.tim, gradient_northern_time, '-o', 'DisplayName', 'Northern Gradient');
hold on;
plot(H2O_obs.tim, gradient_southern_time, '-x', 'DisplayName', 'Southern Gradient');
xlabel('Time');
ylabel('Isotopic Gradient (^{13}CH_4)');
title('Temporal Evolution of Isotopic Gradient');
legend;
grid on;
hold off;
