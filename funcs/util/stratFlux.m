% Load H2O Data
H2O_obs = read_monthly_H2O('global_strat_obs_monthly.nc');

% Extract station names from ch4c13_obs
stations = fieldnames(ch4c13_obs.tim);

% Initialize empty matrices for average concentration and latitude
avg_concentration = zeros(length(stations), 1);
latitudes = zeros(length(stations), 1);

% Create a cell array to store temporal_gradient for each station due to varying time series lengths
temporal_gradient = cell(length(stations), 1);

% Loop through each station to calculate average concentration and gradient
for i = 1:length(stations)
    % Averaging concentration across all times for each station
    avg_concentration(i) = mean(ch4c13_obs.obs.(stations{i}));
    
    % Getting latitude for each station (assuming it's constant over time)
    latitudes(i) = ch4c13_obs.lat.(stations{i})(1);
    
    % Calculating temporal gradient (change in concentration over time)
    temporal_gradient{i} = [0; diff(ch4c13_obs.obs.(stations{i}))]';
end


% Calculating average ^13CH₄ concentration for time bins matching H2O_obs
avg_ch4c13 = NaN(size(H2O_obs.obs));

for t = 1:length(H2O_obs.tim)
    concentrations_at_t = cell(size(stations));
    
    for s = 1:length(stations)
        station = stations{s};
        idx = find(ch4c13_obs.tim.(station) == H2O_obs.tim(t));
        
        if isempty(idx)
            concentrations_at_t{s} = NaN;
        else
            concentrations_at_t{s} = ch4c13_obs.obs.(station)(idx);
        end
    end
    
    avg_ch4c13(t) = mean(cell2mat(concentrations_at_t), 'omitnan'); 
end



% Decompose the data to highlight seasonality
window_size = 12; % For monthly data to get an annual moving average
ch4c13_trend = movmean(avg_ch4c13, window_size);
H2O_trend = movmean(H2O_obs.obs, window_size);

ch4c13_detrended = avg_ch4c13 - ch4c13_trend;
H2O_detrended = H2O_obs.obs - H2O_trend;

% Visualize seasonal variation
f1 = figure;
subplot(2,1,1);
plot(H2O_obs.tim, H2O_detrended);
title('Seasonal Variation in H2O concentrations');
subplot(2,1,2);
plot(H2O_obs.tim, ch4c13_detrended);
title('Seasonal Variation in ^{13}CH_4 concentrations');

% Correlate the seasonal variability
months = month(epicTimeToDatetime(H2O_obs.tim)); % Extracting months
monthly_correlation = zeros(12,1);

for m = 1:12
    monthly_correlation(m) = corr(H2O_detrended(months == m), ch4c13_detrended(months == m), 'Type', 'Spearman', 'Rows', 'complete');

end

% Display monthly correlations
f2 = figure;
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

% Extracting data based on latitude bands
tropical_idx = find(latitudes >= -15 & latitudes <= 15);
northern_polar_idx = find(latitudes >= 60 & latitudes <= 90);
southern_polar_idx = find(latitudes >= -90 & latitudes <= -60);


% Initialize matrices for storing gradients
gradient_northern_time = zeros(size(H2O_obs.tim));
gradient_southern_time = zeros(size(H2O_obs.tim));

% Loop through each time point
for t = 1:length(H2O_obs.tim)
    
    % Fetch the observation for a station at time t, if it exists
    mean_tropical = mean(arrayfun(@(i) fetchObservation(i, t, H2O_obs, ch4c13_obs, stations), tropical_idx), 'omitnan');


    mean_northern_polar = mean(arrayfun(@(i) fetchObservation(i, t, H2O_obs, ch4c13_obs, stations), northern_polar_idx), 'omitnan');
    mean_southern_polar = mean(arrayfun(@(i) fetchObservation(i, t, H2O_obs, ch4c13_obs, stations), southern_polar_idx), 'omitnan');

    
    gradient_northern_time(t) = mean_northern_polar - mean_tropical;
    gradient_southern_time(t) = mean_southern_polar - mean_tropical;
end



% Correlation of Northern Gradient with H2O concentration
correlation_northern = corr(gradient_northern_time, H2O_obs.obs, 'Rows', 'complete');

% Correlation of Southern Gradient with H2O concentration
correlation_southern = corr(gradient_southern_time, H2O_obs.obs, 'Rows', 'complete');

disp(['Correlation of Northern Gradient with H2O: ', num2str(correlation_northern)]);
disp(['Correlation of Southern Gradient with H2O: ', num2str(correlation_southern)]);


% Plot the gradients over time
f3 = figure;
plot(H2O_obs.tim, gradient_northern_time, '-o', 'DisplayName', 'Northern Gradient');
hold on;
plot(H2O_obs.tim, gradient_southern_time, '-x', 'DisplayName', 'Southern Gradient');
xlabel('Time');
ylabel('Isotopic Gradient (^{13}CH_4)');
title('Temporal Evolution of Isotopic Gradient');
legend;
grid on;
hold off;


% Assuming ch4c13_obs and H2O_obs have been loaded and gradients calculated
% Convert epic time to datetime if not already done
if ~isa(H2O_obs.tim, 'datetime')
    H2O_obs.tim = epicTimeToDatetime(H2O_obs.tim);
end

% 1. Decomposition

% Extract trend for H2O concentrations
H2O_trend = movmean(H2O_obs.obs, 12); % 12-month window for annual trend
H2O_detrended = H2O_obs.obs - H2O_trend;

% Calculate ch4c13 gradient
gradient_time = gradient_northern_time;
gradient_trend = movmean(gradient_time, 12);
gradient_detrended = gradient_time - gradient_trend;

% 2. Correlation Calculation

% Seasonal Correlation
seasonal_correlation = corr(H2O_detrended, gradient_detrended, 'Type', 'Spearman', 'Rows', 'complete');

% Trend Correlation
trend_correlation = corr(H2O_trend, gradient_trend, 'Type', 'Spearman', 'Rows', 'complete');

% Display Results
disp(['Seasonal Correlation between H2O and ch4c13 gradient: ', num2str(seasonal_correlation)]);
disp(['Trend Correlation between H2O and ch4c13 gradient: ', num2str(trend_correlation)]);

% Optional: Visualization
figure;
subplot(2, 1, 1);
plot(H2O_obs.tim, H2O_detrended, 'b');
hold on;
plot(H2O_obs.tim, gradient_detrended, 'r');
title('Seasonal Variation');
legend('H2O', 'ch4c13 Gradient');

subplot(2, 1, 2);
plot(H2O_obs.tim, H2O_trend, 'b');
hold on;
plot(H2O_obs.tim, gradient_trend, 'r');
title('Trend Over Time');
legend('H2O', 'ch4c13 Gradient');


% Interpolate
ind =169;
data = H2O_obs.obs(ind:end);
t = 1:length(data);
nonNanIndices = ~isnan(data);
dataInterp = interp1(t(nonNanIndices), data(nonNanIndices), t, 'linear', 'extrap');
Y = fft(dataInterp);

% Assuming you already computed Y as the FFT of your data
Pyy = Y .* conj(Y) / length(Y);

% Frequency vector (depends on the length of your data and its sampling rate)
f = linspace(0, 1, length(Y));

% Identify power of seasonality 
% (assuming monthly data, so looking at 1/12 frequency)
[~, idx] = min(abs(f - 1/12)); 
power_of_seasonality = Pyy(idx);
[indices, magnitudes] = manual_peak_finder(Pyy, 0.01);

% If N is the length of your time series:
N = length(Pyy);
frequencies = (0:(N-1))/N;


peak_frequencies = frequencies(indices);

peak_periods = 1 ./ peak_frequencies;

% Assuming magnitudes is the power of the peaks and peak_periods is the period of the peaks
[sorted_magnitudes, sorted_indices] = sort(magnitudes, 'descend');
sorted_periods = peak_periods(sorted_indices);

% Print the results
fprintf('Period (in units of time) - Power\n');
fprintf('------------------------------\n');
for i = 1:length(sorted_magnitudes)
    fprintf('%f - %f\n', sorted_periods(i), sorted_magnitudes(i));
end



function dt = epicTimeToDatetime(epic_time)
    % Convert epic time (days since 01/0/0000) to MATLAB's datetime

    % Reference date
    reference_date = datetime('0000-01-01', 'Format', 'yyyy-MM-dd');

    % Convert
    dt = reference_date + caldays(epic_time);
end


function observation = fetchObservation(station_idx, t, H2O_obs, ch4c13_obs, stations)
    % Check if the time t in H2O_obs exists in the current station's times
    if ismember(H2O_obs.tim(t), ch4c13_obs.tim.(stations{station_idx}))
        % If it does, retrieve the observation corresponding to that time
        idx = find(ch4c13_obs.tim.(stations{station_idx}) == H2O_obs.tim(t), 1, 'first');
        observation = ch4c13_obs.obs.(stations{station_idx})(idx);
    else
        % If it doesn't, return NaN
        observation = NaN;
    end
end


function [peak_indices, peak_magnitudes] = manual_peak_finder(Pyy, threshold)
    % Initialize empty arrays to hold results
    peak_indices = [];
    peak_magnitudes = [];

    % Iterate through the data to find peaks
    for i = 2:length(Pyy)-1
        if Pyy(i) > threshold && Pyy(i) > Pyy(i-1) && Pyy(i) > Pyy(i+1)
            peak_indices = [peak_indices; i];
            peak_magnitudes = [peak_magnitudes; Pyy(i)];
        end
    end
end
