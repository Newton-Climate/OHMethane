H2O_obs = read_monthly_H2O('global_strat_obs_monthly.nc')

% Assuming ch4c13_obs and H2O_obs have the same time inte`rvals
stations = fieldnames(ch4c13_obs.Tim);

% Initialize empty matrices for average concentration, latitude, and temporal gradient
avg_concentration = zeros(length(stations), 1);
latitudes = zeros(length(stations), 1);
temporal_gradient = zeros(length(stations), length(ch4c13_obs.Tim.(stations{1})));

% Loop through each station to calculate average concentration and gradient
for i = 1:length(stations)
    % Averaging concentration across all times for each station
    avg_concentration(i) = mean(ch4c13_obs.observations.(stations{i}));
    
    % Getting latitude for each station (assuming it's constant over time)
    latitudes(i) = ch4c13_obs.lat.(stations{i})(1);
    
    % Calculating temporal gradient (change in concentration over time)
    temporal_gradient(i, :) = [0, diff(ch4c13_obs.observations.(stations{i}))];
end

% Calculating average ^13CH₄ concentration for time bins matching H2O_obs
avg_ch4c13 = zeros(size(H2O_obs.observations));

for t = 1:length(H2O_obs.Tim)
    % Averaging across all stations for this time point
    concentrations_at_t = arrayfun(@(station) ch4c13_obs.observations.(stations{station})(t), 1:length(stations));
    avg_ch4c13(t) = mean(concentrations_at_t, 'omitnan'); 
end

% Decompose the data to highlight seasonality
window_size = 12; % For monthly data to get an annual moving average
ch4c13_trend = movmean(avg_ch4c13, window_size);
H2O_trend = movmean(H2O_obs.observations, window_size);

ch4c13_detrended = avg_ch4c13 - ch4c13_trend;
H2O_detrended = H2O_obs.observations - H2O_trend;

% Visualize seasonal variation
figure;
subplot(2,1,1);
plot(H2O_obs.Tim, H2O_detrended);
title('Seasonal Variation in H2O concentrations');
subplot(2,1,2);
plot(H2O_obs.Tim, ch4c13_detrended);
title('Seasonal Variation in ^{13}CH_4 concentrations');

% Correlate the seasonal variability
months = month(H2O_obs.Tim); % Extracting months
monthly_correlation = zeros(12,1);

for m = 1:12
    monthly_correlation(m) = corr(H2O_detrended(months == m), ch4c13_detrended(months == m), 'Type', 'Spearman');
end

% Display monthly correlations
figure;
bar(monthly_correlation);
xlabel('Month');
ylabel('Correlation between H2O and ^{13}CH_4');
title('Monthly Correlation between H2O concentrations and ^{13}CH_4');

% Display latitudinal gradient
figure;
plot(latitudes, avg_concentration, 'o-');
xlabel('Latitude');
ylabel('Average ^{13}CH_4 Concentration');
title('Latitudinal Gradient of ^{13}CH_4 Concentration');
