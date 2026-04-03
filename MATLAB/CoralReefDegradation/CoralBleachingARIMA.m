clc;
clear;
close all;


pwd
folderPath = '../../data/raw/CORIS_DATA/';
files = dir(fullfile(folderPath, '**', '*Benthic*'));
if isempty(files)
    error('No benthic files found in %s - check folderPath.', folderPath);
end

fprintf('Found %d file(s), loading...\n', length(files));

raw_tables   = {};
all_col_sets = {};

for f = 1:length(files)
    fpath = fullfile(files(f).folder, files(f).name);
    fprintf('  %s\n', files(f).name);
    try
        tbl = readtable(fpath, 'VariableNamingRule', 'preserve', 'TextType', 'string');


        if height(tbl) > 0 && ismember('YEAR', tbl.Properties.VariableNames)
            if isnan(str2double(tbl.YEAR(1)))
                tbl(1,:) = [];
            end
        end

        raw_tables{end+1}   = tbl;
        all_col_sets{end+1} = tbl.Properties.VariableNames;
        fprintf('    %d rows\n', height(tbl));
    catch ME
        fprintf('  skipped %s - %s\n', files(f).name, ME.message);
    end
end

if isempty(raw_tables)
    error('No files loaded successfully.');
end


all_cols = unique([all_col_sets{:}], 'stable');

for f = 1:length(raw_tables)
    missing_cols = setdiff(all_cols, raw_tables{f}.Properties.VariableNames);
    for mc = 1:length(missing_cols)
        raw_tables{f}.(missing_cols{mc}) = repmat("", height(raw_tables{f}), 1);
    end
    raw_tables{f} = raw_tables{f}(:, all_cols);
end

data_all = vertcat(raw_tables{:});
fprintf('Total rows: %d   Columns: %d\n', height(data_all), width(data_all));

col_names = data_all.Properties.VariableNames;
fprintf('\nColumns:\n');
for c = 1:length(col_names)
    fprintf('  %d. %s\n', c, col_names{c});
end


num_cols = {'YEAR','HARDBOTTOM_P','SOFTBOTTOM_P','RUBBLE_P','latitude', ...
            'longitude','MIN_DEPTH','MAX_DEPTH','PROT','WTD_RUG','STATION_NR'};
for nc = 1:length(num_cols)
    col = num_cols{nc};
    if ismember(col, col_names)
        data_all.(col) = str2double(data_all.(col));
    end
end


hard_coral_codes = [
    "ACR CERV","ACR PALM","AGA AGAR","AGA FRAG","AGA HUMI","AGA LAMA","AGA SPE.", ...
    "COL NATA","COL SPE.","DEN CYLI","DIC STOK","DIP LABY","DIP SPE.","EUS FAST", ...
    "FAV FRAG","HEL CUCU","ISO SINU","MAD AURE","MAD CARM","MAD DECA","MAD FORM", ...
    "MAD SPE.","MAN AREO","MEA MEAN","MEA SPE.","MIL SPE.","MON CAVE","MUS ANGU", ...
    "MYC ALIC","MYC FERO","MYC LAMA","OCU DIFF","OCU ROBU","ORB ANNU","ORB FAVE", ...
    "ORB FRAN","OTH CORA","POR ASTR","POR DIVA","POR FURC","POR PORI","POR SPE.", ...
    "PSE CLIV","PSE STRI","SCO CUBE","SCO SPE.","SID RADI","SID SIDE","SID SPE.", ...
    "SOL BOUR","SOL SPE.","STE INTE"
];

algae_codes = [
    "DIC SPE.","HAL SPE.","LOB SPE.","MAC FLES","MAC CALC","RHO CRUS", ...
    "CYA SPE.","TUR SEDI","TUR FREE","MAG SPE.","PEY SPE.","RAM SPE."
];

sponge_codes = ["SPO OTHE","CLI SPE."];


n_rows = height(data_all);

if ismember('PRIMARY_SAMPLE_UNIT', col_names)
    psu_col = data_all.PRIMARY_SAMPLE_UNIT;
else
    psu_col = repmat("0", n_rows, 1);
end


if ismember('STATION_NR', col_names)
    stn_col = string(data_all.STATION_NR);
else
    stn_col = repmat("0", n_rows, 1);
end

year_col_num = data_all.YEAR;
year_col_str = string(year_col_num);

transect_id      = year_col_str + "_" + psu_col + "_" + stn_col;
unique_transects = unique(transect_id);
n_transects      = length(unique_transects);
fprintf('\nUnique transects: %d\n', n_transects);


T = struct();
T.year        = zeros(n_transects, 1);
T.psu         = strings(n_transects, 1);
T.station     = strings(n_transects, 1);
T.lat         = NaN(n_transects, 1);
T.lon         = NaN(n_transects, 1);
T.hard_coral  = zeros(n_transects, 1);
T.algae       = zeros(n_transects, 1);
T.sponge      = zeros(n_transects, 1);
T.bare        = zeros(n_transects, 1);
T.n_coral_spp = zeros(n_transects, 1);
T.depth_m     = NaN(n_transects, 1);
T.rugosity    = NaN(n_transects, 1);
T.protection  = NaN(n_transects, 1);
T.sub_region  = strings(n_transects, 1);
T.habitat     = strings(n_transects, 1);
T.zone        = strings(n_transects, 1);

hb_col  = data_all.HARDBOTTOM_P;
cat_col = data_all.COVER_CAT_CD;

for i = 1:n_transects
    rows = find(transect_id == unique_transects(i));

    T.year(i)    = year_col_num(rows(1));
    T.psu(i)     = psu_col(rows(1));
    T.station(i) = stn_col(rows(1));

    if ismember('latitude',  col_names), T.lat(i) = data_all.latitude(rows(1));  end
    if ismember('longitude', col_names), T.lon(i) = data_all.longitude(rows(1)); end

    if ismember('MIN_DEPTH', col_names) && ismember('MAX_DEPTH', col_names)
        T.depth_m(i) = (data_all.MIN_DEPTH(rows(1)) + data_all.MAX_DEPTH(rows(1))) / 2;
    end
    if ismember('WTD_RUG',        col_names), T.rugosity(i)   = mean(data_all.WTD_RUG(rows),'omitnan');       end
    if ismember('PROT',           col_names), T.protection(i)  = data_all.PROT(rows(1));                      end
    if ismember('SUB_REGION_NAME',col_names), T.sub_region(i)  = data_all.SUB_REGION_NAME(rows(1));           end
    if ismember('HABITAT_CD',     col_names), T.habitat(i)     = data_all.HABITAT_CD(rows(1));                end
    if ismember('ZONE_NAME',      col_names), T.zone(i)        = data_all.ZONE_NAME(rows(1));                 end

    cats = cat_col(rows);
    hb   = hb_col(rows);

    hc_mask   = ismember(cats, hard_coral_codes);
    alg_mask  = ismember(cats, algae_codes);
    spo_mask  = ismember(cats, sponge_codes);
    bare_mask = (cats == "BAR SUB.");

    T.hard_coral(i)  = sum(hb(hc_mask),   'omitnan');
    T.algae(i)       = sum(hb(alg_mask),  'omitnan');
    T.sponge(i)      = sum(hb(spo_mask),  'omitnan');
    T.bare(i)        = sum(hb(bare_mask), 'omitnan');
    T.n_coral_spp(i) = sum(hc_mask);
end

transect_tbl = struct2table(T);


transect_tbl = transect_tbl(~isnan(transect_tbl.depth_m) & transect_tbl.depth_m > 0, :);
fprintf('Transects with valid depth: %d\n', height(transect_tbl));

years_all = unique(transect_tbl.year);
years_all = years_all(~isnan(years_all));

fprintf('\n%-6s %-8s %-8s %-8s %-10s %-10s\n', ...
    'Year','n','HardCor%','Algae%','Sponge%','CoralSpp');
fprintf('%s\n', repmat('-',1,58));
for y = 1:length(years_all)
    idx = transect_tbl.year == years_all(y);
    fprintf('%-6d %-8d %-8.2f %-8.2f %-10.2f %-10.2f\n', ...
        years_all(y), sum(idx), ...
        mean(transect_tbl.hard_coral(idx),'omitnan'), ...
        mean(transect_tbl.algae(idx),     'omitnan'), ...
        mean(transect_tbl.sponge(idx),    'omitnan'), ...
        mean(transect_tbl.n_coral_spp(idx),'omitnan'));
end

annual_mean_hc  = arrayfun(@(y) mean(transect_tbl.hard_coral(transect_tbl.year==y),'omitnan'), years_all);
annual_se_hc    = arrayfun(@(y) std(transect_tbl.hard_coral(transect_tbl.year==y),'omitnan') / ...
                   sqrt(sum(transect_tbl.year==y)), years_all);
annual_mean_alg = arrayfun(@(y) mean(transect_tbl.algae(transect_tbl.year==y),'omitnan'), years_all);
annual_mean_spo = arrayfun(@(y) mean(transect_tbl.sponge(transect_tbl.year==y),'omitnan'), years_all);


figure('Position', [50 50 1400 800], 'Name', 'Raw Benthic Time Series');

subplot(2,2,1);
errorbar(years_all, annual_mean_hc, annual_se_hc, 'o-', 'LineWidth', 2, ...
    'MarkerSize', 9, 'MarkerFaceColor', [0.1 0.4 0.8], 'Color', [0.1 0.4 0.8]);
xlabel('Year'); ylabel('Hard Coral Cover (%)');
title('Hard Coral Cover');
grid on;
text(years_all, annual_mean_hc + annual_se_hc + 0.3, ...
    arrayfun(@(n) sprintf('n=%d',n), ...
    arrayfun(@(y) sum(transect_tbl.year==y), years_all), 'UniformOutput', false), ...
    'HorizontalAlignment','center','FontSize',8);

subplot(2,2,2);
plot(years_all, annual_mean_alg, 's-', 'LineWidth', 2, 'MarkerSize', 9, ...
    'MarkerFaceColor', [0.2 0.7 0.2], 'Color', [0.2 0.7 0.2]);
hold on;
plot(years_all, annual_mean_spo, 'd-', 'LineWidth', 2, 'MarkerSize', 9, ...
    'MarkerFaceColor', [0.9 0.5 0.1], 'Color', [0.9 0.5 0.1]);
xlabel('Year'); ylabel('Cover (%)');
title('Algae and Sponge Cover');
legend({'Algae','Sponge'},'Location','best');
grid on;

subplot(2,2,3);
scatter(transect_tbl.algae, transect_tbl.hard_coral, 30, transect_tbl.year, ...
    'filled', 'MarkerFaceAlpha', 0.5);
colormap(subplot(2,2,3), cool);
cb = colorbar; cb.Label.String = 'Year';
p_corr = polyfit(transect_tbl.algae, transect_tbl.hard_coral, 1);
hold on;
xr = linspace(0, max(transect_tbl.algae), 100);
plot(xr, polyval(p_corr, xr), 'r--', 'LineWidth', 2);
r = corr(transect_tbl.algae, transect_tbl.hard_coral, 'rows', 'complete');
text(5, max(transect_tbl.hard_coral)*0.9, sprintf('r = %.3f', r), 'FontSize', 10);
xlabel('Algae Cover (%)'); ylabel('Hard Coral Cover (%)');
title('Coral vs Algae');
grid on;

subplot(2,2,4);
scatter(transect_tbl.depth_m, transect_tbl.hard_coral, 25, transect_tbl.year, ...
    'filled', 'MarkerFaceAlpha', 0.4);
colormap(subplot(2,2,4), cool);
colorbar;
xlabel('Depth (m)'); ylabel('Hard Coral Cover (%)');
title('Coral Cover vs Depth');
grid on;

sgtitle('Benthic Cover - Exploratory');


fprintf('\nPanel regression...\n');

year_base = min(years_all);
transect_tbl.year_norm = transect_tbl.year - year_base;

depth_mean = mean(transect_tbl.depth_m,'omitnan');
depth_std  = std(transect_tbl.depth_m, 'omitnan');
alg_mean   = mean(transect_tbl.algae,  'omitnan');
alg_std    = std(transect_tbl.algae,   'omitnan');

transect_tbl.depth_z = (transect_tbl.depth_m - depth_mean) ./ depth_std;
transect_tbl.algae_z = (transect_tbl.algae   - alg_mean)   ./ alg_std;
transect_tbl.prot_z  = transect_tbl.protection;

use_idx = ~isnan(transect_tbl.depth_z) & ~isnan(transect_tbl.algae_z) & ...
          ~isnan(transect_tbl.hard_coral) & ~isnan(transect_tbl.prot_z);

y_panel = transect_tbl.hard_coral(use_idx);
X_panel = [ones(sum(use_idx),1), ...
           transect_tbl.year_norm(use_idx), ...
           transect_tbl.depth_z(use_idx), ...
           transect_tbl.algae_z(use_idx), ...
           transect_tbl.prot_z(use_idx)];

[b_ols, ~, ~, ~, stats_ols] = regress(y_panel, X_panel);

fprintf('  N = %d   R2 = %.4f   F = %.2f   p = %.4f\n', ...
    sum(use_idx), stats_ols(1), stats_ols(2), stats_ols(3));

coef_names = {'Intercept','Year trend (per yr)','Depth z','Algae z','Protection'};
for c = 1:length(b_ols)
    fprintf('  %-22s: %8.4f\n', coef_names{c}, b_ols(c));
end
fprintf('  Trend: %.4f%% per year (covariate-adjusted)\n', b_ols(2));


fprintf('\nARIMA on annual time series...\n');

y_annual         = annual_mean_hc;
y_algae_annual   = annual_mean_alg;
forecast_horizon = 6;

t_ann   = (1:length(y_annual))';
p_trend = polyfit(t_ann, y_annual, 1);

if p_trend(1) < 0
    fprintf('  Declining trend: %.4f%%/obs - keeping (real ecological signal).\n', p_trend(1));
else
    fprintf('  Positive trend: %.4f%%/obs\n', p_trend(1));
end

trend_component = polyval(p_trend, t_ann);
y_detrended     = y_annual - trend_component;

model_candidates = {};
aic_vals         = [];
model_names_list = {};

try
    fit_ar1 = estimate(arima(1,0,0), y_detrended, 'Display','off');
    model_candidates{end+1} = fit_ar1;
    aic_vals(end+1)         = summarize(fit_ar1).AIC;
    model_names_list{end+1} = 'AR(1)';
    fprintf('  AR(1)  AIC = %.4f\n', aic_vals(end));
catch ME
    fprintf('  AR(1) failed: %s\n', ME.message);
end

try
    fit_ma1 = estimate(arima(0,0,1), y_detrended, 'Display','off');
    model_candidates{end+1} = fit_ma1;
    aic_vals(end+1)         = summarize(fit_ma1).AIC;
    model_names_list{end+1} = 'MA(1)';
    fprintf('  MA(1)  AIC = %.4f\n', aic_vals(end));
catch ME
    fprintf('  MA(1) failed: %s\n', ME.message);
end

try
    fit_wn = estimate(arima(0,0,0), y_detrended, 'Display','off');
    model_candidates{end+1} = fit_wn;
    aic_vals(end+1)         = summarize(fit_wn).AIC;
    model_names_list{end+1} = 'White Noise';
    fprintf('  WN     AIC = %.4f\n', aic_vals(end));
catch ME
    fprintf('  WN failed: %s\n', ME.message);
end

if isempty(model_candidates)
    error('No ARIMA models converged.');
end

[~, best_idx]   = min(aic_vals);
best_fit        = model_candidates{best_idx};
best_model_name = model_names_list{best_idx};
fprintf('\n  Best: %s (AIC = %.4f)\n', best_model_name, aic_vals(best_idx));

[yF_det, yMSE] = forecast(best_fit, forecast_horizon, 'Y0', y_detrended);

last_t       = length(y_annual);
future_t     = (last_t+1):(last_t+forecast_horizon);
trend_future = polyval(p_trend, future_t');

yF       = max(0, min(100, yF_det + trend_future));
yF_upper = min(100, max(0, yF_det + trend_future + 1.96*sqrt(yMSE)));
yF_lower = max(0,        yF_det + trend_future - 1.96*sqrt(yMSE));

last_year        = years_all(end);
future_years_vec = (last_year+1):(last_year+forecast_horizon);

fprintf('\n  %-6s %-10s %-10s %-10s\n','Year','Forecast','Lower95','Upper95');
fprintf('  %s\n', repmat('-',1,40));
for f = 1:forecast_horizon
    fprintf('  %-6d %-10.2f %-10.2f %-10.2f\n', ...
        future_years_vec(f), yF(f), yF_lower(f), yF_upper(f));
end


fprintf('\nARIMAX with algae as exogenous predictor...\n');

results_arimax = struct();

try
    X_alg = (y_algae_annual - mean(y_algae_annual)) / std(y_algae_annual);

    fit_ax1 = estimate(regARIMA('ARLags',1), y_annual, 'X', X_alg, 'Display','off');
    fit_ax2 = estimate(regARIMA('MALags',1), y_annual, 'X', X_alg, 'Display','off');

    aic_ax1 = summarize(fit_ax1).AIC;
    aic_ax2 = summarize(fit_ax2).AIC;
    fprintf('  RegARIMAX+AR(1): AIC = %.4f\n', aic_ax1);
    fprintf('  RegARIMAX+MA(1): AIC = %.4f\n', aic_ax2);

    if aic_ax1 <= aic_ax2
        best_arimax = fit_ax1;  arimax_name = 'RegARIMAX+AR(1)';
    else
        best_arimax = fit_ax2;  arimax_name = 'RegARIMAX+MA(1)';
    end

    fprintf('  Best: %s   algae coef: %.4f\n', arimax_name, best_arimax.Beta);

    t_alg = (1:length(y_algae_annual))';
    p_alg = polyfit(t_alg, y_algae_annual, 1);
    future_t_alg     = (length(y_algae_annual)+1):(length(y_algae_annual)+forecast_horizon);
    alg_forecast_raw = polyval(p_alg, future_t_alg');
    alg_forecast_z   = (alg_forecast_raw - mean(y_algae_annual)) / std(y_algae_annual);

    yF_arimax = forecast(best_arimax, forecast_horizon, ...
        'Y0', y_annual, 'X0', X_alg, 'XF', alg_forecast_z);

    res_ax     = double(infer(best_arimax, y_annual, 'X', X_alg));
    res_std_ax = std(res_ax);

    results_arimax.forecast       = max(0, min(100, yF_arimax));
    results_arimax.forecast_upper = min(100, max(0, yF_arimax + 1.96*res_std_ax));
    results_arimax.forecast_lower = max(0,        yF_arimax - 1.96*res_std_ax);
    results_arimax.algae_forecast = alg_forecast_raw;
    results_arimax.model_name     = arimax_name;
    results_arimax.beta_algae     = best_arimax.Beta;
    results_arimax.success        = true;

    fprintf('\n  %-6s %-10s %-10s %-12s\n','Year','Forecast','Algae%','Change');
    for f = 1:forecast_horizon
        fprintf('  %-6d %-10.2f %-10.2f %-12.2f\n', future_years_vec(f), ...
            results_arimax.forecast(f), alg_forecast_raw(f), ...
            results_arimax.forecast(f) - y_annual(end));
    end

catch ME
    fprintf('  ARIMAX failed: %s - using ARIMA only.\n', ME.message);
    results_arimax.success    = false;
    results_arimax.forecast   = yF;
    results_arimax.model_name = 'ARIMA (ARIMAX fallback)';
end


fprintf('\nRegional analysis...\n');

sub_regions   = unique(transect_tbl.sub_region);
sub_regions   = sub_regions(sub_regions ~= "");
sub_yr_count  = arrayfun(@(s) length(unique(transect_tbl.year(transect_tbl.sub_region==s))), sub_regions);
model_regions = sub_regions(sub_yr_count >= 3);
fprintf('  Regions with >= 3 years: %d\n', length(model_regions));

regional_results = struct();

figure('Position', [50 50 1600 900], 'Name', 'Regional Forecasts');
n_plot_rows = ceil(length(model_regions) / 3);
n_plot_cols = min(3, length(model_regions));

for s = 1:length(model_regions)
    rname = model_regions(s);
    field = matlab.lang.makeValidName(rname);
    fprintf('\n  %s\n', rname);

    sub_data  = transect_tbl(transect_tbl.sub_region == rname, :);
    sub_years = unique(sub_data.year)';
    y_reg  = arrayfun(@(y) mean(sub_data.hard_coral(sub_data.year==y),'omitnan'), sub_years);
    se_reg = arrayfun(@(y) std(sub_data.hard_coral(sub_data.year==y),'omitnan') / ...
             sqrt(sum(sub_data.year==y)), sub_years);
    n_reg  = arrayfun(@(y) sum(sub_data.year==y), sub_years);

    fprintf('    years: %s   mean: %.2f%%\n', num2str(sub_years), mean(y_reg));

    t_reg     = (1:length(y_reg))';
    p_reg     = polyfit(t_reg, y_reg, 1);
    y_det_reg = y_reg - polyval(p_reg, t_reg);

    try
        if length(y_reg) >= 4
            fit_reg = estimate(arima(1,0,0), y_det_reg, 'Display','off');
        else
            fit_reg = estimate(arima(0,0,0), y_det_reg, 'Display','off');
        end

        [yF_det_reg, yMSE_reg] = forecast(fit_reg, forecast_horizon, 'Y0', y_det_reg);

        ft_reg   = (length(y_reg)+1):(length(y_reg)+forecast_horizon);
        trend_ft = polyval(p_reg, ft_reg');
        yF_reg    = max(0, min(100, yF_det_reg + trend_ft));
        yF_reg_up = min(100, max(0, yF_det_reg + trend_ft + 1.96*sqrt(yMSE_reg)));
        yF_reg_lo = max(0,        yF_det_reg + trend_ft - 1.96*sqrt(yMSE_reg));

        future_yr_reg = (sub_years(end)+1):(sub_years(end)+forecast_horizon);

        regional_results.(field).sub_region     = rname;
        regional_results.(field).years          = sub_years;
        regional_results.(field).observed       = y_reg;
        regional_results.(field).se             = se_reg;
        regional_results.(field).trend_slope    = p_reg(1);
        regional_results.(field).forecast       = yF_reg;
        regional_results.(field).forecast_upper = yF_reg_up;
        regional_results.(field).forecast_lower = yF_reg_lo;
        regional_results.(field).forecast_years = future_yr_reg;

        fprintf('    trend: %.4f%%/yr\n', p_reg(1));

        subplot(n_plot_rows, n_plot_cols, s);
        hold on;
        errorbar(sub_years, y_reg, se_reg, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
            'MarkerFaceColor','b','DisplayName','Observed');
        yr_all    = [sub_years(:); future_yr_reg(:)];
        trend_all = polyval(p_reg, (1:length(yr_all))');
        plot(yr_all, trend_all, 'k--', 'LineWidth', 1, 'DisplayName','Trend');
        plot(future_yr_reg, yF_reg, 'r-o', 'LineWidth', 2, 'MarkerSize', 7, ...
            'MarkerFaceColor','r','DisplayName','Forecast');
        fill([future_yr_reg, fliplr(future_yr_reg)], ...
             [yF_reg_up', fliplr(yF_reg_lo')], ...
             'r','FaceAlpha',0.15,'EdgeColor','none','DisplayName','95% CI');
        xlabel('Year'); ylabel('Hard Coral (%)');
        title(sprintf('%s (avg n=%d)', strtrim(rname), round(mean(n_reg))));
        legend('Location','best','FontSize',7);
        grid on;
        xlim([min(sub_years)-1, max(future_yr_reg)+0.5]);
        ylim([0, max([y_reg(:); yF_reg_up(:)])*1.2 + 2]);

    catch ME
        fprintf('    error: %s\n', ME.message);
    end
end

sgtitle('Hard Coral Cover - Regional ARIMA Forecasts');


figure('Position', [50 50 1400 700], 'Name', 'Forecast');

subplot(1,2,1);
hold on;
errorbar(years_all, annual_mean_hc, annual_se_hc, 'b-o', 'LineWidth', 2.5, ...
    'MarkerSize', 10, 'MarkerFaceColor', [0.1 0.4 0.8], 'Color', [0.1 0.4 0.8], ...
    'DisplayName','Observed ± SE');
trend_x = linspace(years_all(1), future_years_vec(end), 100);
trend_t = linspace(1, length(y_annual)+forecast_horizon, 100);
plot(trend_x, polyval(p_trend, trend_t), 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Trend (%.4f%%/yr)', p_trend(1)));
plot(future_years_vec, yF, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 9, ...
    'MarkerFaceColor','r','DisplayName', sprintf('ARIMA (%s)', best_model_name));
fill([future_years_vec, fliplr(future_years_vec)], [yF_upper', fliplr(yF_lower')], ...
    'r','FaceAlpha',0.15,'EdgeColor','none','DisplayName','95% CI');
if results_arimax.success
    plot(future_years_vec, results_arimax.forecast, 'g-s', 'LineWidth', 2, 'MarkerSize', 8, ...
        'MarkerFaceColor','g','DisplayName', results_arimax.model_name);
end
xlabel('Year'); ylabel('Hard Coral Cover (%)');
title('Hard Coral Forecast');
legend('Location','best','FontSize',9);
grid on;
ylim([0, max([annual_mean_hc; yF_upper])*1.3 + 1]);
xlim([years_all(1)-1, future_years_vec(end)+0.5]);

subplot(1,2,2);
hold on;
plot(years_all, annual_mean_alg, 'g-o', 'LineWidth', 2.5, 'MarkerSize', 10, ...
    'MarkerFaceColor',[0.1 0.7 0.1],'DisplayName','Observed Algae');
if results_arimax.success
    plot(future_years_vec, results_arimax.algae_forecast, 'g--s', 'LineWidth', 2, ...
        'MarkerSize',8,'DisplayName','Projected Algae');
end
yyaxis right
plot(years_all, annual_mean_hc, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor','b','DisplayName','Observed Coral');
plot(future_years_vec, yF, 'r--', 'LineWidth', 2, 'DisplayName','Forecast Coral');
ylabel('Hard Coral Cover (%)');
yyaxis left;
ylabel('Algae Cover (%)');
xlabel('Year');
title('Algae-Coral Dynamics');
legend('Location','best','FontSize',9);
grid on;

sgtitle('Benthic Cover - ARIMA + ARIMAX Forecasts');


top_species = ["ORB FAVE","SID SIDE","MON CAVE","POR ASTR","PSE STRI", ...
               "COL NATA","STE INTE","ACR CERV","MIL SPE.","SID RADI"];
top_species_names = {'O. faveolata','S. siderea','M. cavernosa','P. astreoides', ...
                     'P. strigosa','C. natans','S. intersepta','A. cervicornis', ...
                     'Millepora spp.','S. radians'};

spp_annual = zeros(length(top_species), length(years_all));
for sp = 1:length(top_species)
    sp_rows = (data_all.COVER_CAT_CD == top_species(sp));
    for y = 1:length(years_all)
        mask = sp_rows & (data_all.YEAR == years_all(y));
        if any(mask)
            spp_annual(sp,y) = mean(data_all.HARDBOTTOM_P(mask),'omitnan');
        end
    end
end

figure('Position', [50 50 1400 600], 'Name', 'Species Trends');

subplot(1,2,1);
area(years_all, spp_annual');
legend(top_species_names,'Location','eastoutside','FontSize',8);
xlabel('Year'); ylabel('Mean Cover (pts)');
title('Top Coral Species Composition');
xlim([years_all(1)-0.5, years_all(end)+0.5]);
grid on;

subplot(1,2,2);
colors_spp = lines(length(top_species));
hold on;
for sp = 1:length(top_species)
    if any(spp_annual(sp,:) > 0)
        plot(years_all, spp_annual(sp,:), 'o-', 'LineWidth', 2, 'MarkerSize', 7, ...
            'Color', colors_spp(sp,:), 'DisplayName', top_species_names{sp});
    end
end
xlabel('Year'); ylabel('Mean Cover (pts)');
title('Per-Species Trends');
legend('Location','eastoutside','FontSize',8);
grid on; xlim([years_all(1)-0.5, years_all(end)+0.5]);

sgtitle('Hard Coral Species Composition');


y_hat     = X_panel * b_ols;
res_panel = y_panel - y_hat;

figure('Position', [50 50 1400 700], 'Name', 'Residual Diagnostics');

subplot(2,3,1);
scatter(y_hat, res_panel, 20, 'filled', 'MarkerFaceAlpha', 0.5);
yline(0, 'r--', 'LineWidth', 2);
xlabel('Fitted (%)'); ylabel('Residuals (%)');
title('Panel: Residuals vs Fitted');
grid on;

subplot(2,3,2);
histogram(res_panel, 40, 'Normalization','pdf','FaceColor',[0.6 0.6 0.9]);
hold on;
xr = linspace(min(res_panel), max(res_panel), 100);
plot(xr, normpdf(xr, mean(res_panel), std(res_panel)), 'r-', 'LineWidth', 2);
xlabel('Residuals (%)'); ylabel('Density');
title('Residual Distribution');
legend('Data','Normal','Location','best');
grid on;

subplot(2,3,3);
qqplot(res_panel);
title('Q-Q Plot');

try
    res_arima = double(infer(best_fit, y_detrended));
    subplot(2,3,4);
    stem(years_all, res_arima, 'LineWidth', 2);
    yline(0,'r--');
    xlabel('Year'); ylabel('Residuals (%)');
    title(sprintf('ARIMA Residuals (%s)', best_model_name));
    grid on;

    subplot(2,3,5);
    if length(res_arima) >= 4
        autocorr(res_arima);
        title('Residual ACF');
    else
        text(0.5, 0.5, 'Too few points for ACF', ...
            'HorizontalAlignment','center','Units','normalized');
        axis off; title('ACF (insufficient data)');
    end
catch
    subplot(2,3,4); axis off;
    subplot(2,3,5); axis off;
end

subplot(2,3,6);
axis off;
text(0.05, 0.95, sprintf([ ...
    'Diagnostics\n\n' ...
    'Panel (N=%d)\n' ...
    '  mean res: %.4f%%\n' ...
    '  std res:  %.4f%%\n' ...
    '  R2:       %.4f\n\n' ...
    'ARIMA: %s\n' ...
    '  trend:    %.4f%%/yr\n' ...
    '  AIC:      %.4f'], ...
    sum(use_idx), mean(res_panel), std(res_panel), stats_ols(1), ...
    best_model_name, p_trend(1), aic_vals(best_idx)), ...
    'Units','normalized','VerticalAlignment','top','FontSize',9,'FontName','Courier');

sgtitle('Residual Diagnostics');


fprintf('\n%s\n', repmat('-',1,72));
fprintf('Forecast summary\n');
fprintf('%s\n', repmat('-',1,72));
fprintf('%-8s %-14s %-14s %-14s %-14s\n','Year','ARIMA Fcst','Lower95','Upper95','ARIMAX Fcst');
fprintf('%s\n', repmat('-',1,72));
for f = 1:forecast_horizon
    arimax_val = NaN;
    if results_arimax.success, arimax_val = results_arimax.forecast(f); end
    fprintf('%-8d %-14.2f %-14.2f %-14.2f %-14.2f\n', ...
        future_years_vec(f), yF(f), yF_lower(f), yF_upper(f), arimax_val);
end
fprintf('%s\n', repmat('-',1,72));
fprintf('Trend: %.4f%% per year\n', p_trend(1));
if results_arimax.success
    fprintf('ARIMAX algae coef: %.4f%% per SD\n', results_arimax.beta_algae);
end

fprintf('\n%-35s %-8s %-10s %-12s\n','Sub-Region','N yrs','Trend/yr','Fcst mean%');
fprintf('%s\n', repmat('-',1,72));
for s = 1:length(model_regions)
    field = matlab.lang.makeValidName(model_regions(s));
    if isfield(regional_results, field)
        fprintf('%-35s %-8d %-10.4f %-12.2f\n', ...
            strtrim(model_regions(s)), ...
            length(regional_results.(field).years), ...
            regional_results.(field).trend_slope, ...
            mean(regional_results.(field).forecast));
    end
end


fprintf('\nExporting...\n');

global_tbl = table(future_years_vec', yF, yF_lower, yF_upper, ...
    'VariableNames',{'Year','ARIMA_Forecast_HardCoral_Pct','Lower_95CI','Upper_95CI'});
if results_arimax.success
    global_tbl.ARIMAX_Forecast     = results_arimax.forecast;
    global_tbl.ARIMAX_Lower        = results_arimax.forecast_lower;
    global_tbl.ARIMAX_Upper        = results_arimax.forecast_upper;
    global_tbl.Projected_Algae_Pct = results_arimax.algae_forecast;
end
writetable(global_tbl, 'benthic_cover_forecast.csv');
fprintf('  benthic_cover_forecast.csv\n');

for s = 1:length(model_regions)
    field = matlab.lang.makeValidName(model_regions(s));
    if isfield(regional_results, field)
        reg_tbl = table( ...
            regional_results.(field).forecast_years', ...
            regional_results.(field).forecast, ...
            regional_results.(field).forecast_lower, ...
            regional_results.(field).forecast_upper, ...
            'VariableNames',{'Year','Forecast_HardCoral_Pct','Lower_95CI','Upper_95CI'});
        fname = sprintf('benthic_forecast_%s.csv', ...
            strrep(strtrim(model_regions(s)), ' ', '_'));
        writetable(reg_tbl, fname);
        fprintf('  %s\n', fname);
    end
end

writetable(table(coef_names', b_ols, 'VariableNames',{'Predictor','Coefficient'}), ...
    'benthic_panel_regression_coefficients.csv');
fprintf('  benthic_panel_regression_coefficients.csv\n');

fprintf('\nDone.\n');
if p_trend(1) < 0
    fprintf('Declining trend (%.4f%%/yr) - monitor closely.\n', p_trend(1));
else
    fprintf('Stable or slight increase (%.4f%%/yr).\n', p_trend(1));
end