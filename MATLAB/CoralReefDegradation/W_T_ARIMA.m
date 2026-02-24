clc;
clear;
close all;

%% ── Publication figure defaults ─────────────────────────────────────────
set(groot, 'defaultAxesColor',        'w');
set(groot, 'defaultFigureColor',      'w');
set(groot, 'defaultAxesBox',          'off');
set(groot, 'defaultAxesXGrid',        'off');
set(groot, 'defaultAxesYGrid',        'off');
set(groot, 'defaultAxesFontName',     'Helvetica');
set(groot, 'defaultAxesFontSize',     11);
set(groot, 'defaultAxesFontColor',    'k');
set(groot, 'defaultAxesXColor',       'k');
set(groot, 'defaultAxesYColor',       'k');
set(groot, 'defaultAxesZColor',       'k');
set(groot, 'defaultAxesLineWidth',    0.8);
set(groot, 'defaultTextFontName',     'Helvetica');
set(groot, 'defaultTextColor',        'k');
set(groot, 'defaultLegendBox',        'off');
set(groot, 'defaultLegendTextColor',  'k');
set(groot, 'defaultLegendFontSize',   10);

BLUE  = [0.13 0.39 0.68];
RED   = [0.80 0.17 0.17];
GREEN = [0.13 0.60 0.32];
LGRAY = [0.85 0.85 0.85];

%% ── Load data ────────────────────────────────────────────────────────────
file = '../../data/raw/NOAA/sst.mnmean.v4.nc';

lon  = ncread(file, 'lon');
lat  = ncread(file, 'lat');
time = ncread(file, 'time');
sst  = ncread(file, 'sst');

fprintf('Data dimensions:\n');
fprintf('  Longitude: %d points (%.1f to %.1f deg)\n', length(lon), min(lon), max(lon));
fprintf('  Latitude:  %d points (%.1f to %.1f deg)\n', length(lat), min(lat), max(lat));
fprintf('  Time:      %d months\n', length(time));
fprintf('  SST array: %dx%dx%d\n', size(sst,1), size(sst,2), size(sst,3));

base_date = datetime(1800,1,1);
dates     = base_date + days(time);
fprintf('Time range: %s to %s\n', datestr(dates(1)), datestr(dates(end)));

%% ── Clean missing data ───────────────────────────────────────────────────
sst(sst < -100) = NaN;

%% ── Extract tropical band ────────────────────────────────────────────────
tropical_idx  = lat >= -23.5 & lat <= 23.5;
lat_tropical  = lat(tropical_idx);
sst_tropical  = sst(:, tropical_idx, :);

fprintf('Tropical latitudes: %.1f S to %.1f N\n', min(lat_tropical), max(lat_tropical));

%% ── Tropical mean time series ────────────────────────────────────────────
tropical_mean_sst = squeeze(mean(sst_tropical, [1 2], 'omitnan'));

fprintf('Tropical SST range: %.2f to %.2f deg C\n', min(tropical_mean_sst), max(tropical_mean_sst));
fprintf('Mean: %.2f deg C,  Std: %.2f deg C\n', mean(tropical_mean_sst), std(tropical_mean_sst));

%% ── Helper: scrub every axes in a figure to be publication-clean ─────────
%   Fixes box, grid, tick direction, colors, and font on all axes and text
%   objects, including those created internally by autocorr/parcorr/qqplot.
function cleanAxes(fig)
    axList = findall(fig, 'Type', 'axes');
    for k = 1:numel(axList)
        ax = axList(k);
        set(ax, 'Box','off', 'XGrid','off', 'YGrid','off', ...
            'TickDir','out', 'XColor','k', 'YColor','k', 'ZColor','k', ...
            'FontName','Helvetica', 'FontSize',11, 'Color','w');
    end
    txtList = findall(fig, 'Type', 'text');
    set(txtList, 'Color','k', 'FontName','Helvetica');
    legList = findall(fig, 'Type', 'legend');
    for k = 1:numel(legList)
        set(legList(k), 'Box','off', 'TextColor','k');
    end
    set(fig, 'Color', 'w');
end

%% ── Figure 1 · Time series overview ─────────────────────────────────────
fig1 = figure('Position', [100 100 1400 900]);

subplot(4,1,1);
plot(dates, tropical_mean_sst, '-', 'Color', BLUE, 'LineWidth', 1.2);
xlabel('Year', 'FontWeight','bold');
ylabel('SST (deg C)', 'FontWeight','bold');
title('Tropical Band Mean SST (23.5 S to 23.5 N) - Full Record', 'FontWeight','bold');
datetick('x','yyyy','keeplimits');

subplot(4,1,2);
plot(dates(year(dates)>=1950), tropical_mean_sst(year(dates)>=1950), ...
    '-', 'Color', RED, 'LineWidth', 1.2);
xlabel('Year', 'FontWeight','bold');
ylabel('SST (deg C)', 'FontWeight','bold');
title('Tropical SST (1950-Present)', 'FontWeight','bold');
datetick('x','yyyy','keeplimits');

subplot(4,1,3);
plot(dates(year(dates)>=2000), tropical_mean_sst(year(dates)>=2000), ...
    '-', 'Color', GREEN, 'LineWidth', 1.2);
xlabel('Year', 'FontWeight','bold');
ylabel('SST (deg C)', 'FontWeight','bold');
title('Tropical SST (2000-Present)', 'FontWeight','bold');
datetick('x','yyyy','keeplimits');

t_numeric = 1:length(tropical_mean_sst);
p = polyfit(t_numeric', tropical_mean_sst, 1);
trend = polyval(p, t_numeric);

subplot(4,1,4);
hold on;
plot(dates, tropical_mean_sst, '-', 'Color', BLUE, 'LineWidth', 1.0, 'DisplayName', 'Observed');
plot(dates, trend, '--', 'Color', RED, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Trend: %.4f deg C/month (%.2f deg C/century)', p(1), p(1)*12*100));
xlabel('Year', 'FontWeight','bold');
ylabel('SST (deg C)', 'FontWeight','bold');
title('Tropical SST with Linear Trend', 'FontWeight','bold');
legend('Location','northwest');
datetick('x','yyyy','keeplimits');

cleanAxes(fig1);
fprintf('\nLinear trend: %.4f deg C/month (%.2f deg C/century)\n', p(1), p(1)*12*100);

%% ── Figure 2 · ACF / PACF ───────────────────────────────────────────────
y = double(tropical_mean_sst);

fig2 = figure('Position', [100 100 1200 800]);

subplot(3,1,1);
plot(y, '-', 'Color', BLUE, 'LineWidth', 1.0);
xlabel('Month Index', 'FontWeight','bold');
ylabel('SST (deg C)', 'FontWeight','bold');
title('Tropical SST Time Series', 'FontWeight','bold');

subplot(3,1,2);
autocorr(y);
title('Autocorrelation Function (ACF)', 'FontWeight','bold');

subplot(3,1,3);
parcorr(y);
title('Partial Autocorrelation Function (PACF)', 'FontWeight','bold');

cleanAxes(fig2);

%% ── ARIMA modelling ──────────────────────────────────────────────────────
fprintf('\n-- ARIMA MODELLING WITH TREND --\n');

t_numeric   = (1:length(y))';
p_trend     = polyfit(t_numeric, y, 1);
y_detrended = y - polyval(p_trend, t_numeric);

fprintf('Trend removed: %.6f deg C/month  (%.4f/yr, %.2f/century)\n', ...
    p_trend(1), p_trend(1)*12, p_trend(1)*12*100);

model1 = arima(1,0,1); fit1 = estimate(model1, y_detrended, 'Display','off');
model2 = arima(2,0,1); fit2 = estimate(model2, y_detrended, 'Display','off');
model3 = arima(1,0,2); fit3 = estimate(model3, y_detrended, 'Display','off');
model4 = arima(2,0,2); fit4 = estimate(model4, y_detrended, 'Display','off');
model5 = arima('ARLags',1,'MALags',1,'SARLags',12,'SMALags',12);
fit5   = estimate(model5, y_detrended, 'Display','off');

aics = [summarize(fit1).AIC, summarize(fit2).AIC, summarize(fit3).AIC, ...
        summarize(fit4).AIC, summarize(fit5).AIC];
model_names = {'ARIMA(1,0,1)','ARIMA(2,0,1)','ARIMA(1,0,2)', ...
               'ARIMA(2,0,2)','SARIMA(1,0,1)x(1,0,1,12)'};

fprintf('\nModel AIC comparison:\n');
for m = 1:5, fprintf('  %-35s %.2f\n', model_names{m}, aics(m)); end

[~, best_idx] = min(aics);
fits = {fit1,fit2,fit3,fit4,fit5};
best_fit = fits{best_idx};
fprintf('\nBest: %s  (AIC = %.2f)\n', model_names{best_idx}, aics(best_idx));

%% ── Figure 3 · Forecast ──────────────────────────────────────────────────
numPeriods     = 240;
[yF_det, yMSE] = forecast(best_fit, numPeriods, 'Y0', y_detrended);

future_t_idx = (length(y)+1 : length(y)+numPeriods)';
yF           = yF_det + polyval(p_trend, future_t_idx);
yF_upper     = yF + 1.96*sqrt(yMSE);
yF_lower     = yF - 1.96*sqrt(yMSE);
future_dates = dates(end) + calmonths(1:numPeriods)';

fig3 = figure('Position', [100 100 1400 560]);
hold on;

fill([future_dates; flipud(future_dates)], [yF_upper; flipud(yF_lower)], ...
    RED, 'FaceAlpha',0.15, 'EdgeColor','none', 'DisplayName','95% CI');

plot_start = find(year(dates) >= year(dates(end))-50, 1);
dates_all  = [dates(plot_start:end); future_dates];
trend_all  = polyval(p_trend, [plot_start:length(y), future_t_idx']');
plot(dates_all, trend_all, '--', 'Color','k', 'LineWidth',1.4, 'DisplayName','Linear trend');

plot(dates(plot_start:end), y(plot_start:end), '-', 'Color',BLUE, ...
    'LineWidth',1.6, 'DisplayName','Historical');
plot(future_dates, yF, '-', 'Color',RED, 'LineWidth',2.0, 'DisplayName','Forecast');

xlabel('Date', 'FontWeight','bold');
ylabel('Tropical Mean SST (deg C)', 'FontWeight','bold');
title(sprintf('Tropical SST Forecast - %s + Linear Trend', model_names{best_idx}), ...
    'FontWeight','bold');
legend('Location','northwest');
datetick('x','yyyy','keeplimits');

cleanAxes(fig3);

fprintf('\nForecast summary:\n');
for mo = [12 60 120 240]
    fprintf('  %s: %.2f +/- %.2f deg C\n', datestr(future_dates(mo),'mmm-yyyy'), ...
        yF(mo), 1.96*sqrt(yMSE(mo)));
end
fprintf('Projected increase over 20 years: %.2f deg C\n', yF(240)-y(end));

%% ── Figure 4 · Residual diagnostics ─────────────────────────────────────
res = double(infer(best_fit, y_detrended));
[h_lb, pValue] = lbqtest(res, 'Lags', 20);

fig4 = figure('Position', [100 100 1400 800]);

subplot(2,3,1);
plot(dates, res, '-', 'Color','k', 'LineWidth',0.6);
yline(0, '--', 'Color',RED, 'LineWidth',1.5);
xlabel('Date', 'FontWeight','bold');
ylabel('Residuals (deg C)', 'FontWeight','bold');
title('Residuals Over Time', 'FontWeight','bold');
datetick('x','yyyy','keeplimits');

subplot(2,3,2);
histogram(res, 50, 'Normalization','pdf', 'FaceColor',LGRAY, 'EdgeColor','none');
hold on;
xr = linspace(min(res), max(res), 100);
plot(xr, normpdf(xr, mean(res), std(res)), '-', 'Color',RED, 'LineWidth',2);
xlabel('Residuals (deg C)', 'FontWeight','bold');
ylabel('Density', 'FontWeight','bold');
title('Residual Distribution', 'FontWeight','bold');
legend({'Residuals','Normal'}, 'Location','best');

subplot(2,3,3);
qqplot(res);
title('Q-Q Plot of Residuals', 'FontWeight','bold');

subplot(2,3,4);
autocorr(res);
title('Residual ACF', 'FontWeight','bold');

subplot(2,3,5);
parcorr(res);
title('Residual PACF', 'FontWeight','bold');

subplot(2,3,6);
axis off;
if h_lb == 0, wn_str = 'White noise (pass)'; else, wn_str = 'Autocorrelation detected'; end
stats_lines = { ...
    'Residual Statistics', ...
    sprintf('Mean:    %.4f deg C', mean(res)), ...
    sprintf('Std Dev: %.4f deg C', std(res)), ...
    sprintf('Min:     %.4f deg C', min(res)), ...
    sprintf('Max:     %.4f deg C', max(res)), ...
    '', ...
    'Ljung-Box (20 lags)', ...
    sprintf('p = %.4f', pValue), ...
    wn_str};
for k = 1:numel(stats_lines)
    fw = 'normal';
    if k == 1 || k == 7, fw = 'bold'; end
    text(0.05, 1.0-(k-1)*0.10, stats_lines{k}, ...
        'Units','normalized', 'FontSize',10, 'FontWeight',fw, ...
        'Color','k', 'VerticalAlignment','top');
end

cleanAxes(fig4);

fprintf('\nResidual diagnostics:\n');
fprintf('  Mean: %.4f,  Std: %.4f deg C\n', mean(res), std(res));
fprintf('  Ljung-Box p-value: %.4f\n', pValue);

%% ── Export forecast ──────────────────────────────────────────────────────
forecast_table = table(future_dates, yF, yF_lower, yF_upper, sqrt(yMSE), ...
    'VariableNames', {'Date','Forecast_SST','Lower_95CI','Upper_95CI','Std_Error'});
writetable(forecast_table, 'tropical_sst_forecast.csv');
fprintf('\nForecast saved to tropical_sst_forecast.csv\n');