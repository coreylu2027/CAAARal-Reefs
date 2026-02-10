clc;
clear;
close all;
% parameters
data = readtable("../../data/processed/NOAA/Surface_Data_Temperature_Cleaned.csv");
data = data(2:end, :);
data.Properties.VariableNames = {'Index', 'Time', 'Latitude', 'Longitude', 'Temperature', 'Error', 'Anomaly'};

% Converting string to number
if iscell(data.Temperature)
    data.Temperature = str2double(data.Temperature);
    data.Latitude = str2double(data.Latitude);
    data.Longitude = str2double(data.Longitude);
    data.Error = str2double(data.Error);
    data.Anomaly = str2double(data.Anomaly);
end


fprintf('Loaded %d observations\n', height(data));

fprintf('\n--- data summary ---\n');
fprintf('Temperature range: %.2f to %.2f °C\n', min(data.Temperature), max(data.Temperature));
fprintf('Mean temperature: %.2f °C\n', mean(data.Temperature));
fprintf('Std deviation: %.2f °C\n', std(data.Temperature));
fprintf('Latitude range: %.2f to %.2f\n', min(data.Latitude), max(data.Latitude));
fprintf('Longitude range: %.2f to %.2f\n', min(data.Longitude), max(data.Longitude));

%% Figure 1: Temperature Distribution
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
histogram(data.Temperature, 50, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Temperature (°C)');
ylabel('Frequency');
title('Temperature Distribution');
grid on;

subplot(1,3,2);
histogram(data.Anomaly, 50, 'FaceColor', [0.8 0.4 0.2]);
xlabel('Temperature Anomaly (°C)');
ylabel('Frequency');
title('Anomaly Distribution');
grid on;

subplot(1,3,3);
histogram(data.Error, 30, 'FaceColor', [0.4 0.8 0.4]);
xlabel('Error (°C)');
ylabel('Frequency');
title('Measurement Error Distribution');
grid on;

sgtitle('Statistical Distributions');

%% Figure 2: Spatial Distribution (2D map)
figure('Position', [100, 100, 1000, 600]);
scatter(data.Longitude, data.Latitude, 20, data.Temperature, 'filled');
colorbar;
colormap(jet);
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Sea Surface Temperature Map (2024-11-27)');
grid on;
axis equal;


%% Figures 3: temp vs lat. and long.
figure('Position', [100, 100, 1200, 400]);

subplot(1,2,1);
scatter(data.Latitude, data.Temperature, 5, data.Temperature, 'filled');
colorbar;
colormap(jet);
xlabel('Latitude (°)');
ylabel('Temperature (°C)');
title('Temperature vs Latitude');
grid on;

subplot(1,2,2);
scatter(data.Longitude, data.Temperature, 5, data.Temperature, 'filled');
colorbar;
colormap(jet);
xlabel('Longitude (°)');
ylabel('Temperature (°C)');
title('Temperature vs Longitude');
grid on;

%% 4. spatial analysis by latitude bins

% Bin data by latitude
lat_bins = -40:5:40;
lat_centers = lat_bins(1:end-1) + 2.5;
temp_by_lat = nan(size(lat_centers));

for i = 1:length(lat_centers)
    idx = data.Latitude >= lat_bins(i) & data.Latitude < lat_bins(i+1);
    if sum(idx) > 0
        temp_by_lat(i) = mean(data.Temperature(idx), 'omitnan');
    end
end

% Remove NaN values for plotting
valid_idx = ~isnan(temp_by_lat);

figure('Position', [100, 100, 800, 500]);
plot(lat_centers(valid_idx), temp_by_lat(valid_idx), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Latitude (°)');
ylabel('Mean Temperature (°C)');
title('Mean Temperature by Latitude Band');
grid on;

%% 5. ARIMA modeling
% For spatial data, we'll create a "transect" along a specific latitude or longitude
fprintf('\nARIMA MODELING\n');

% Example: Create time series along a specific latitude (around -39.875°)
target_lat = -39.875;
transect_data = data(abs(data.Latitude - target_lat) < 0.01, :);
transect_data = sortrows(transect_data, 'Longitude');

% Extract temperature series
y = transect_data.Temperature;
n = length(y);

fprintf('Created transect with %d points at latitude %.2f°\n', n, target_lat);

% Plot the transect
figure('Position', [100, 100, 1200, 800]);
subplot(3,1,1);
plot(transect_data.Longitude, y, 'o-', 'LineWidth', 1.5);
xlabel('Longitude (°)');
ylabel('Temperature (°C)');
title(sprintf('Temperature Transect at Latitude %.2f°', target_lat));
grid on;

% Check ACF and PACF
subplot(3,1,2);
autocorr(y);
title('Autocorrelation Function (ACF)');

subplot(3,1,3);
parcorr(y);
title('Partial Autocorrelation Function (PACF)');