clc;
clear;
close all;

%% 1. load data
data = readtable('../../data/raw/NOAA/global_bleaching_environmental.csv');

fprintf('Total observations: %d\n', height(data));
fprintf('Columns: %d\n', width(data));

%% 2. Define regions

% Define 6 major ocean regions based on lat/lon

lat = data.Latitude_Degrees;
lon = data.Longitude_Degrees;

% Initialize region labels
region = strings(height(data), 1);

% Classify regions
for i = 1:height(data)
    if lat(i) >= 5 && lat(i) <= 30 && lon(i) >= -100 && lon(i) <= -60
        region(i) = 'Caribbean/Gulf';
    elseif lat(i) > 0 && lon(i) >= -80 && lon(i) <= 0
        region(i) = 'North Atlantic';
    elseif lat(i) < 0 && lon(i) >= -60 && lon(i) <= 20
        region(i) = 'South Atlantic';
    elseif lat(i) > 0 && (lon(i) > 120 || lon(i) < -120)
        region(i) = 'North Pacific';
    elseif lat(i) < 0 && (lon(i) > 120 || lon(i) < -120)
        region(i) = 'South Pacific';
    elseif lon(i) >= 40 && lon(i) <= 120
        region(i) = 'Indian Ocean';
    else
        region(i) = 'Other';
    end
end

data.Region = region;

% Count observations per region
fprintf('\nObservations by region:\n');
region_counts = groupsummary(data, 'Region');
for i = 1:height(region_counts)
    fprintf('  %s: %d (%.1f%%)\n', region_counts.Region{i}, ...
        region_counts.GroupCount(i), ...
        100*region_counts.GroupCount(i)/height(data));
end

%% 3. Create dates and filter

% Create date column
data.Date_Full = datetime(data.Date_Year, data.Date_Month, data.Date_Day);

% Remove rows with missing Percent_Bleaching or invalid dates
valid_idx = ~isnan(data.Percent_Bleaching) & ~isnat(data.Date_Full);
data_clean = data(valid_idx, :);

fprintf('\nAfter filtering for valid bleaching percentages and dates:\n');
fprintf('  Total observations: %d\n', height(data_clean));

%% 4. aggregate by 3 months

% Get major regions (exclude 'Other')
major_regions = {'Caribbean/Gulf', 'North Atlantic', 'South Atlantic', ...
                 'North Pacific', 'South Pacific', 'Indian Ocean'};

% Filter to major regions only
data_filtered = data_clean(ismember(data_clean.Region, major_regions), :);

fprintf('\nMajor regions only: %d observations\n', height(data_filtered));

% Create QUARTERLY time series for better data density
fprintf('\nAggregating QUARTERLY bleaching data...\n');

% Get date range
min_date = min(data_filtered.Date_Full);
max_date = max(data_filtered.Date_Full);

% Create quarterly bins (better coverage than monthly)
quarterly_dates = (dateshift(min_date, 'start', 'quarter'):calmonths(3):dateshift(max_date, 'end', 'quarter'))';

fprintf('Date range: %s to %s (%d quarters)\n', ...
    datestr(min_date), datestr(max_date), length(quarterly_dates));

%% 5. Create time series

regional_ts = struct();

for r = 1:length(major_regions)
    region_name = major_regions{r};
    region_data = data_filtered(data_filtered.Region == region_name, :);
    
    % Aggregate by quarter
    quarterly_bleach = zeros(length(quarterly_dates), 1);
    quarterly_count = zeros(length(quarterly_dates), 1);
    
    for q = 1:length(quarterly_dates)
        quarter_start = quarterly_dates(q);
        quarter_end = quarter_start + calmonths(3) - days(1);
        
        idx = region_data.Date_Full >= quarter_start & region_data.Date_Full <= quarter_end;
        
        if sum(idx) > 0
            quarterly_bleach(q) = mean(region_data.Percent_Bleaching(idx), 'omitnan');
            quarterly_count(q) = sum(idx);
        else
            quarterly_bleach(q) = NaN;
        end
    end
    
    % Store in structure
    regional_ts.(matlab.lang.makeValidName(region_name)).dates = quarterly_dates;
    regional_ts.(matlab.lang.makeValidName(region_name)).bleaching = quarterly_bleach;
    regional_ts.(matlab.lang.makeValidName(region_name)).count = quarterly_count;
    regional_ts.(matlab.lang.makeValidName(region_name)).name = region_name;
    
    fprintf('  %s: %d quarters with data (%.1f%% coverage)\n', ...
        region_name, sum(~isnan(quarterly_bleach)), ...
        100*sum(~isnan(quarterly_bleach))/length(quarterly_dates));
end

%% 6. Visualize the regions

figure('Position', [100, 100, 1600, 1000]);

for r = 1:length(major_regions)
    subplot(3, 2, r);
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    dates = regional_ts.(field_name).dates;
    bleach = regional_ts.(field_name).bleaching;
    
    % Plot with valid data only
    valid_idx = ~isnan(bleach);
    plot(dates(valid_idx), bleach(valid_idx), 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
    
    xlabel('Year');
    ylabel('Mean Bleaching (%)');
    title(region_name);
    grid on;
    ylim([0 100]);
    datetick('x', 'yyyy', 'keeplimits');
end

sgtitle('Coral Bleaching Time Series by Ocean Region (Quarterly)');

%% 7. arima modeling with trend added,

fprintf('\nARIMA MODELING BY REGION\n');
fprintf('NOTE: Bleaching is episodic (triggered by SST anomalies), unlike continuous SST data.\n');
fprintf('Sparse quarterly data may show artificial negative trends.\n');
fprintf('Using floor of zero on negative trends to avoid unrealistic forecasts.\n\n');

forecast_periods = 8;  % 8 quarters = 2 years
results = struct();

figure('Position', [100, 100, 1600, 1000]);

for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    fprintf('\n=== %s ===\n', region_name);
    
    dates = regional_ts.(field_name).dates;
    bleach = regional_ts.(field_name).bleaching;
    
    % Remove NaN values and get continuous segments
    valid_idx = ~isnan(bleach);
    y_raw = bleach(valid_idx);
    dates_valid = dates(valid_idx);
    
    % Convert to double
    y = double(y_raw);
    
    % Check if we have enough data
    if length(y) < 12
        fprintf('  Insufficient data (n=%d), skipping ARIMA\n', length(y));
        continue;
    end
    
    fprintf('  Data points: %d\n', length(y));
    fprintf('  Mean bleaching: %.2f%%\n', mean(y));
    fprintf('  Std: %.2f%%\n', std(y));
    fprintf('  Range: %.2f%% to %.2f%%\n', min(y), max(y));
    
    % Estimate and remove linear trend
    t_numeric = (1:length(y))';
    p_trend_raw = polyfit(t_numeric, y, 1);
    
    % Floor negative trends at zero (use constant mean instead)
    % Bleaching is driven by ocean warming (see SST analysis), so shouldn't decrease long-term
    if p_trend_raw(1) < 0
        fprintf('  Negative trend detected (%.3f%% per quarter) - using constant mean instead\n', p_trend_raw(1));
        fprintf('  (Bleaching driven by rising SST, negative trends are artifacts of sparse data)\n');
        p_trend = [0, mean(y)];  % Zero slope, mean intercept
    else
        p_trend = p_trend_raw;
        fprintf('  Linear trend: %.3f%% per quarter (%.2f%% per year)\n', ...
            p_trend(1), p_trend(1)*4);
    end
    
    trend_component = polyval(p_trend, t_numeric);
    y_detrended = y - trend_component;
    
    % Try ARIMA models on detrended data
    try
        % ARIMA(1,0,0) - AR(1)
        model1 = arima(1, 0, 0);
        fit1 = estimate(model1, y_detrended, 'Display', 'off');
        aic1 = summarize(fit1).AIC;
        
        % ARIMA(1,0,1)
        model2 = arima(1, 0, 1);
        fit2 = estimate(model2, y_detrended, 'Display', 'off');
        aic2 = summarize(fit2).AIC;
        
        % ARIMA(2,0,1)
        model3 = arima(2, 0, 1);
        fit3 = estimate(model3, y_detrended, 'Display', 'off');
        aic3 = summarize(fit3).AIC;
        
        % ARIMA(0,0,1) - MA(1)
        model4 = arima(0, 0, 1);
        fit4 = estimate(model4, y_detrended, 'Display', 'off');
        aic4 = summarize(fit4).AIC;
        
        % Select best model
        [~, best_idx] = min([aic1, aic2, aic3, aic4]);
        models = {fit1, fit2, fit3, fit4};
        model_names = {'AR(1)', 'ARIMA(1,0,1)', 'ARIMA(2,0,1)', 'MA(1)'};
        best_fit = models{best_idx};
        
        fprintf('  Best model: %s (AIC=%.2f)\n', model_names{best_idx}, ...
            min([aic1, aic2, aic3, aic4]));
        
        % Forecast on detrended data
        [yF_detrended, yMSE] = forecast(best_fit, forecast_periods, 'Y0', y_detrended);
        
        % Add trend back
        last_t = length(y);
        future_t = (last_t + 1):(last_t + forecast_periods);
        trend_future = polyval(p_trend, future_t');
        yF = yF_detrended + trend_future;
        
        % Confidence intervals
        yF_upper = yF + 1.96*sqrt(yMSE);
        yF_lower = yF - 1.96*sqrt(yMSE);
        
        % Bound forecasts to [0, 100]
        yF(yF < 0) = 0;
        yF(yF > 100) = 100;
        yF_lower(yF_lower < 0) = 0;
        yF_upper(yF_upper > 100) = 100;
        
        % Create future dates
        last_date = dates_valid(end);
        future_dates = last_date + calmonths(3*(1:forecast_periods))';
        
        % Store results
        results.(field_name).model = best_fit;
        results.(field_name).model_name = model_names{best_idx};
        results.(field_name).forecast = yF;
        results.(field_name).forecast_dates = future_dates;
        results.(field_name).forecast_upper = yF_upper;
        results.(field_name).forecast_lower = yF_lower;
        results.(field_name).trend_slope = p_trend(1);
        results.(field_name).historical_data = y;
        results.(field_name).historical_dates = dates_valid;
        
        % Plot
        subplot(3, 2, r);
        hold on;
        
        % Historical data
        plot(dates_valid, y, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, ...
            'DisplayName', 'Historical');
        
        % Trend line
        trend_all = polyval(p_trend, [t_numeric; future_t']);
        dates_all = [dates_valid; future_dates];
        if p_trend(1) == 0
            plot(dates_all, trend_all, 'k--', 'LineWidth', 1, 'DisplayName', 'Mean');
        else
            plot(dates_all, trend_all, 'k--', 'LineWidth', 1, 'DisplayName', 'Trend');
        end
        
        % Forecast
        plot(future_dates, yF, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, ...
            'DisplayName', 'Forecast');
        
        % Confidence interval
        fill([future_dates; flipud(future_dates)], ...
             [yF_upper; flipud(yF_lower)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% CI');
        
        xlabel('Year');
        ylabel('Bleaching (%)');
        if p_trend(1) == 0
            title(sprintf('%s - %s (const)', region_name, model_names{best_idx}));
        else
            title(sprintf('%s - %s + Trend', region_name, model_names{best_idx}));
        end
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        ylim([0 min(100, max([y; yF_upper])*1.1)]);
        datetick('x', 'yyyy', 'keeplimits');
        
    catch ME
        fprintf('  Error fitting ARIMA: %s\n', ME.message);
    end
end

sgtitle('ARIMA Forecasts by Ocean Region (2 years ahead, quarterly)');

%% 8. Residual diagnostics

figure('Position', [100, 100, 1600, 1000]);

plot_count = 0;
for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    if ~isfield(results, field_name)
        continue;
    end
    
    plot_count = plot_count + 1;
    
    y = results.(field_name).historical_data;
    
    % Detrend
    t_numeric = (1:length(y))';
    p_trend_raw = polyfit(t_numeric, y, 1);
    if p_trend_raw(1) < 0
        p_trend = [0, mean(y)];
    else
        p_trend = p_trend_raw;
    end
    trend_component = polyval(p_trend, t_numeric);
    y_detrended = y - trend_component;
    
    % Get residuals
    res = double(infer(results.(field_name).model, y_detrended));
    
    % Plot ACF
    subplot(3, 2, plot_count);
    autocorr(res);
    title(sprintf('%s - Residual ACF', region_name));
end

sgtitle('Residual Diagnostics by Region (Detrended Data)');

%% 9. summary table

fprintf('\n FORECAST SUMMARY \n');
fprintf('%-20s %-15s %10s %12s %12s %15s\n', ...
    'Region', 'Model', 'Trend/yr', 'Current', 'Forecast', 'Change');
fprintf('%s\n', repmat('-', 1, 90));

for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    if isfield(results, field_name)
        y = results.(field_name).historical_data;
        yF = results.(field_name).forecast;
        trend = results.(field_name).trend_slope * 4;  % Convert to per year
        
        current_avg = mean(y(end-min(3, length(y)-1):end));
        forecast_avg = mean(yF);
        change = forecast_avg - current_avg;
        
        fprintf('%-20s %-15s %9.2f%% %10.2f%% %10.2f%% %+13.2f%%\n', ...
            region_name, ...
            results.(field_name).model_name, ...
            trend, ...
            current_avg, ...
            forecast_avg, ...
            change);
    else
        fprintf('%-20s %-15s %10s %12s %12s %15s\n', ...
            region_name, 'N/A', '-', '-', '-', '-');
    end
end

%% 10. regional comparison

figure('Position', [100, 100, 1400, 600]);

% Current vs Forecast comparison
subplot(1,2,1);
region_current = [];
region_forecast = [];
region_labels = {};

for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    if isfield(results, field_name)
        y = results.(field_name).historical_data;
        yF = results.(field_name).forecast;
        
        region_current(end+1) = mean(y(end-min(3, length(y)-1):end));
        region_forecast(end+1) = mean(yF);
        region_labels{end+1} = region_name;
    end
end

x = 1:length(region_labels);
width = 0.35;
bar(x - width/2, region_current, width, 'DisplayName', 'Current (last 4 qtrs)');
hold on;
bar(x + width/2, region_forecast, width, 'DisplayName', 'Forecast (next 8 qtrs)');
set(gca, 'XTick', x, 'XTickLabel', region_labels);
xtickangle(45);
ylabel('Mean Bleaching (%)');
title('Current vs Forecast Bleaching by Region');
legend('Location', 'best');
grid on;

% Trend comparison
subplot(1,2,2);
region_trends = [];
trend_labels = {};

for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    if isfield(results, field_name)
        region_trends(end+1) = results.(field_name).trend_slope * 4;  % Per year
        trend_labels{end+1} = region_name;
    end
end

bar(region_trends);
set(gca, 'XTickLabel', trend_labels);
xtickangle(45);
ylabel('Trend (% per year)');
title('Bleaching Trend by Region');
grid on;
yline(0, 'r--', 'LineWidth', 1);

%% 11. export results


for r = 1:length(major_regions)
    region_name = major_regions{r};
    field_name = matlab.lang.makeValidName(region_name);
    
    if isfield(results, field_name)
        % Create table
        forecast_table = table(results.(field_name).forecast_dates, ...
                               results.(field_name).forecast, ...
                               results.(field_name).forecast_lower, ...
                               results.(field_name).forecast_upper, ...
                               'VariableNames', {'Date', 'Forecast_Bleaching', ...
                               'Lower_95CI', 'Upper_95CI'});
        
        % Save to CSV
        filename = sprintf('bleaching_forecast_%s.csv', ...
            strrep(region_name, '/', '_'));
        writetable(forecast_table, filename);
        fprintf('  Saved: %s\n', filename);
    end
end