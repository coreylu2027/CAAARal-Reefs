%% Combined SST and Coral Bleaching Analysis

clc;
clear;
close all;

%% 1. LOAD BOTH DATASETS

% Load SST data
sst_file = '../../data/raw/NOAA/sst.mnmean.v4.nc';
lon = ncread(sst_file, 'lon');
lat = ncread(sst_file, 'lat');
time = ncread(sst_file, 'time');
sst = ncread(sst_file, 'sst');

% Convert time and clean
base_date = datetime(1800,1,1);
sst_dates = base_date + days(time);
sst(sst < -100) = NaN;

% Extract tropical SST
tropical_idx = lat >= -23.5 & lat <= 23.5;
sst_tropical = sst(:, tropical_idx, :);
tropical_sst_ts = double(squeeze(mean(sst_tropical, [1 2], 'omitnan')));

% Load bleaching data
bleach_data = readtable('../../data/raw/NOAA/global_bleaching_environmental.csv');

% Filter and create dates
bleach_data.Date_Full = datetime(bleach_data.Date_Year, bleach_data.Date_Month, bleach_data.Date_Day);
valid_idx = ~isnan(bleach_data.Percent_Bleaching) & ~isnat(bleach_data.Date_Full);
bleach_data = bleach_data(valid_idx, :);

fprintf('SST data: %d months from %s to %s\n', length(sst_dates), ...
    datestr(sst_dates(1)), datestr(sst_dates(end)));
fprintf('Bleaching data: %d observations from %s to %s\n', height(bleach_data), ...
    datestr(min(bleach_data.Date_Full)), datestr(max(bleach_data.Date_Full)));

%% 2. ALIGN SST TO BLEACHING OBSERVATION DATES

% For each bleaching observation, find the SST in the same month
bleach_data.SST_Match = NaN(height(bleach_data), 1);

fprintf('\nMatching bleaching observations to SST...\n');

for i = 1:height(bleach_data)
    obs_date = bleach_data.Date_Full(i);
    
    % Find closest SST month
    [~, idx] = min(abs(sst_dates - obs_date));
    
    % Only match if within same month
    if year(sst_dates(idx)) == year(obs_date) && month(sst_dates(idx)) == month(obs_date)
        bleach_data.SST_Match(i) = tropical_sst_ts(idx);
    end
end

% Remove rows without SST match
bleach_with_sst = bleach_data(~isnan(bleach_data.SST_Match), :);
fprintf('Matched %d bleaching observations to SST\n', height(bleach_with_sst));

%% 3. VISUALIZE SST-BLEACHING RELATIONSHIP

figure('Position', [100, 100, 1400, 500]);

subplot(1,2,1);
scatter(bleach_with_sst.SST_Match, bleach_with_sst.Percent_Bleaching, 20, ...
    bleach_with_sst.Date_Year, 'filled', 'MarkerFaceAlpha', 0.6);
colorbar;
colormap(jet);
xlabel('Tropical SST (°C)');
ylabel('Bleaching (%)');
title('SST vs Bleaching (colored by year)');
grid on;

% Add regression line
p = polyfit(bleach_with_sst.SST_Match, bleach_with_sst.Percent_Bleaching, 1);
hold on;
sst_range = linspace(min(bleach_with_sst.SST_Match), max(bleach_with_sst.SST_Match), 100);
plot(sst_range, polyval(p, sst_range), 'r--', 'LineWidth', 2);
text(min(sst_range)+1, max(bleach_with_sst.Percent_Bleaching)-10, ...
    sprintf('Slope: %.2f%%/°C', p(1)), 'FontSize', 10);

subplot(1,2,2);
% Bin SST and show mean bleaching
sst_bins = 24:0.5:30;
bleach_means = zeros(length(sst_bins)-1, 1);
bleach_stds = zeros(length(sst_bins)-1, 1);
bin_centers = zeros(length(sst_bins)-1, 1);

for i = 1:length(sst_bins)-1
    idx = bleach_with_sst.SST_Match >= sst_bins(i) & ...
          bleach_with_sst.SST_Match < sst_bins(i+1);
    if sum(idx) > 0
        bleach_means(i) = mean(bleach_with_sst.Percent_Bleaching(idx));
        bleach_stds(i) = std(bleach_with_sst.Percent_Bleaching(idx));
        bin_centers(i) = (sst_bins(i) + sst_bins(i+1)) / 2;
    else
        bleach_means(i) = NaN;
        bleach_stds(i) = NaN;
        bin_centers(i) = (sst_bins(i) + sst_bins(i+1)) / 2;
    end
end

errorbar(bin_centers, bleach_means, bleach_stds, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Tropical SST (°C)');
ylabel('Mean Bleaching (%)');
title('Average Bleaching by SST Bin');
grid on;

%% 4. CLASSIFY BY REGION

fprintf('\nClassifying observations by region...\n');

major_regions = {'Caribbean/Gulf', 'North Pacific', 'South Pacific', 'Indian Ocean'};

lat_b = bleach_with_sst.Latitude_Degrees;
lon_b = bleach_with_sst.Longitude_Degrees;
region = strings(height(bleach_with_sst), 1);

for i = 1:height(bleach_with_sst)
    if lat_b(i) >= 5 && lat_b(i) <= 30 && lon_b(i) >= -100 && lon_b(i) <= -60
        region(i) = 'Caribbean/Gulf';
    elseif lat_b(i) > 0 && (lon_b(i) > 120 || lon_b(i) < -120)
        region(i) = 'North Pacific';
    elseif lat_b(i) < 0 && (lon_b(i) > 120 || lon_b(i) < -120)
        region(i) = 'South Pacific';
    elseif lon_b(i) >= 40 && lon_b(i) <= 120
        region(i) = 'Indian Ocean';
    else
        region(i) = 'Other';
    end
end
bleach_with_sst.Region = region;

% Count by region
region_counts = groupsummary(bleach_with_sst, 'Region');
for i = 1:height(region_counts)
    fprintf('  %s: %d observations\n', region_counts.Region{i}, region_counts.GroupCount(i));
end

%% 5. CREATE QUARTERLY TIME SERIES

min_date = min(bleach_with_sst.Date_Full);
max_date = max(bleach_with_sst.Date_Full);
quarterly_dates = (dateshift(min_date, 'start', 'quarter'):calmonths(3):dateshift(max_date, 'end', 'quarter'))';

fprintf('\nQuarterly aggregation: %d quarters\n', length(quarterly_dates));

%% 6. REGRESSION-BASED ARIMA MODELING

forecast_quarters = 8;
results = struct();

figure('Position', [100, 100, 1600, 1000]);

plot_idx = 0;

for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    fprintf('\n=== %s ===\n', region_name);
    
    % Filter region data
    region_data = bleach_with_sst(bleach_with_sst.Region == region_name, :);
    
    if height(region_data) < 20
        fprintf('  Insufficient data (n=%d), skipping\n', height(region_data));
        continue;
    end
    
    % Aggregate to quarterly
    bleach_q = NaN(length(quarterly_dates), 1);
    sst_q = NaN(length(quarterly_dates), 1);
    
    for q = 1:length(quarterly_dates)
        qstart = quarterly_dates(q);
        qend = qstart + calmonths(3) - days(1);
        
        idx = region_data.Date_Full >= qstart & region_data.Date_Full <= qend;
        
        if sum(idx) > 0
            bleach_q(q) = mean(region_data.Percent_Bleaching(idx));
            sst_q(q) = mean(region_data.SST_Match(idx));
        end
    end
    
    % Remove NaN values
    valid = ~isnan(bleach_q) & ~isnan(sst_q);
    y = double(bleach_q(valid));
    X = double(sst_q(valid));
    dates_valid = quarterly_dates(valid);
    
    if length(y) < 12
        fprintf('  Insufficient quarterly data (n=%d), skipping\n', length(y));
        continue;
    end
    
    fprintf('  Quarterly observations: %d\n', length(y));
    fprintf('  Mean bleaching: %.2f%%, SST: %.2f°C\n', mean(y), mean(X));
    
    plot_idx = plot_idx + 1;
    
    % Regression with ARIMA errors using regARIMA
    try
        % Model 1: Regression + AR(1)
        mdl1 = regARIMA('ARLags', 1);
        fit1 = estimate(mdl1, y, 'X', X, 'Display', 'off');
        
        % Model 2: Regression + MA(1)
        mdl2 = regARIMA('MALags', 1);
        fit2 = estimate(mdl2, y, 'X', X, 'Display', 'off');
        
        % Model 3: Regression + ARMA(1,1)
        mdl3 = regARIMA('ARLags', 1, 'MALags', 1);
        fit3 = estimate(mdl3, y, 'X', X, 'Display', 'off');
        
        % Compare AIC
        aic1 = summarize(fit1).AIC;
        aic2 = summarize(fit2).AIC;
        aic3 = summarize(fit3).AIC;
        
        [~, best_idx] = min([aic1, aic2, aic3]);
        models = {fit1, fit2, fit3};
        model_names = {'Reg+AR(1)', 'Reg+MA(1)', 'Reg+ARMA(1,1)'};
        best_fit = models{best_idx};
        
        fprintf('  Best model: %s (AIC=%.2f)\n', model_names{best_idx}, min([aic1, aic2, aic3]));
        fprintf('  SST coefficient: %.2f\n', best_fit.Beta);
        
        % Forecast future SST (linear trend)
        t_sst = (1:length(X))';
        p_sst = polyfit(t_sst, X, 1);
        future_t = (length(X)+1):(length(X)+forecast_quarters);
        X_forecast = polyval(p_sst, future_t');
        
        % Forecast bleaching
        yF = forecast(best_fit, forecast_quarters, 'Y0', y, 'X0', X, 'XF', X_forecast);
        
        % Estimate confidence intervals from residuals
        res = infer(best_fit, y, 'X', X);
        res_std = std(double(res));
        yF_upper = yF + 1.96*res_std;
        yF_lower = yF - 1.96*res_std;
        
        % Bound to [0, 100]
        yF(yF < 0) = 0; yF(yF > 100) = 100;
        yF_lower(yF_lower < 0) = 0;
        yF_upper(yF_upper > 100) = 100;
        
        % Future dates
        last_date = dates_valid(end);
        future_dates = last_date + calmonths(3*(1:forecast_quarters))';
        
        % Store results
        results.(field_name).model = best_fit;
        results.(field_name).model_name = model_names{best_idx};
        results.(field_name).forecast_bleach = yF;
        results.(field_name).forecast_sst = X_forecast;
        results.(field_name).forecast_dates = future_dates;
        results.(field_name).forecast_upper = yF_upper;
        results.(field_name).forecast_lower = yF_lower;
        results.(field_name).historical_bleach = y;
        results.(field_name).historical_sst = X;
        results.(field_name).historical_dates = dates_valid;
        results.(field_name).sst_coef = best_fit.Beta;
        results.(field_name).sst_trend = p_sst(1);
        
        % Plot
        subplot(2, 2, plot_idx);
        yyaxis left
        hold on;
        plot(dates_valid, y, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Hist Bleaching');
        plot(future_dates, yF, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Forecast Bleach');
        fill([future_dates; flipud(future_dates)], [yF_upper; flipud(yF_lower)], ...
            'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% CI');
        ylabel('Bleaching (%)');
        ylim([0 100]);
        
        yyaxis right
        plot(dates_valid, X, 'k--', 'LineWidth', 1, 'DisplayName', 'Hist SST');
        plot(future_dates, X_forecast, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Forecast SST');
        ylabel('SST (°C)');
        
        title(sprintf('%s - %s', region_name, model_names{best_idx}));
        xlabel('Year');
        legend('Location', 'best', 'FontSize', 7);
        grid on;
        datetick('x', 'yyyy', 'keeplimits');
        
    catch ME
        fprintf('  Error: %s\n', ME.message);
        continue;
    end
end

sgtitle('Regression-ARIMA Forecasts: SST vs Bleaching (2 years ahead)');

%% 7. SUMMARY TABLE

fprintf('\n--- FORECAST SUMMARY ---\n');
fprintf('%-20s %-15s %10s %12s %12s %14s\n', ...
    'Region', 'Model', 'SST Coef', 'SST Trend', 'Future SST', 'Forecast Bleach');
fprintf('%s\n', repmat('-', 1, 95));

for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    if isfield(results, field_name)
        fprintf('%-20s %-15s %9.2f %11.3f°C/Q %10.2f°C %13.2f%%\n', ...
            region_name, ...
            results.(field_name).model_name, ...
            results.(field_name).sst_coef, ...
            results.(field_name).sst_trend, ...
            mean(results.(field_name).forecast_sst), ...
            mean(results.(field_name).forecast_bleach));
    end
end

%% 8. EXPORT FORECASTS

fprintf('\n--- EXPORTING FORECASTS ---\n');

for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    if isfield(results, field_name)
        forecast_table = table(results.(field_name).forecast_dates, ...
                               results.(field_name).forecast_sst, ...
                               results.(field_name).forecast_bleach, ...
                               results.(field_name).forecast_lower, ...
                               results.(field_name).forecast_upper, ...
                               'VariableNames', {'Date', 'Forecast_SST_C', ...
                               'Forecast_Bleaching_Pct', 'Lower_95CI', 'Upper_95CI'});
        
        filename = sprintf('combined_forecast_%s.csv', strrep(region_name, '/', '_'));
        writetable(forecast_table, filename);
        fprintf('  Saved: %s\n', filename);
    end
end

%% 9. RESIDUAL DIAGNOSTICS

figure('Position', [100, 100, 1600, 800]);

subplot_idx = 0;
for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    if ~isfield(results, field_name)
        continue;
    end
    
    subplot_idx = subplot_idx + 1;
    
    % Get residuals
    y = results.(field_name).historical_bleach;
    X = results.(field_name).historical_sst;
    res = double(infer(results.(field_name).model, y, 'X', X));
    
    subplot(2, 2, subplot_idx);
    autocorr(res);
    title(sprintf('%s - Residual ACF', region_name));
end

sgtitle('Residual Diagnostics');

fprintf('\nAnalysis complete!\n');