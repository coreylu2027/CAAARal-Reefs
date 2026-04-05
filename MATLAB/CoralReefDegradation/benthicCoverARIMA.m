clc;
clear;
close all;


set(groot, 'defaultAxesColor',       'w');
set(groot, 'defaultFigureColor',     'w');
set(groot, 'defaultAxesBox',         'off');
set(groot, 'defaultAxesXGrid',       'off');
set(groot, 'defaultAxesYGrid',       'off');
set(groot, 'defaultAxesFontName',    'Helvetica');
set(groot, 'defaultAxesFontSize',    11);
set(groot, 'defaultAxesXColor',      'k');
set(groot, 'defaultAxesYColor',      'k');
set(groot, 'defaultAxesZColor',      'k');
set(groot, 'defaultAxesLineWidth',   0.9);
set(groot, 'defaultTextFontName',    'Helvetica');
set(groot, 'defaultTextColor',       'k');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendBox',       'off');
set(groot, 'defaultLegendTextColor', 'k');
set(groot, 'defaultLegendFontSize',  10);

CLR_CORAL   = [0.13 0.39 0.68];
CLR_ALGAE   = [0.13 0.60 0.32];
CLR_SPONGE  = [0.90 0.62 0.00];
CLR_FCST    = [0.80 0.17 0.17];
CLR_SARIMAX = [0.49 0.18 0.56];
CLR_TREND   = [0.30 0.30 0.30];
CLR_LGRAY   = [0.88 0.88 0.88];
SHOW_ONLY_CORE_FIGURES = true;
SHOW_SPECIES_FIGURE = false;


function cleanFig(fig)
    for ax = findall(fig, 'Type', 'axes')'
        set(ax, 'Box','off', 'XGrid','off', 'YGrid','off', ...
            'TickDir','out', 'LineWidth',0.9, ...
            'XColor','k', 'YColor','k', 'ZColor','k', ...
            'FontName','Helvetica', 'FontSize',11, 'Color','w');
        if isprop(ax,'Title')
            ax.Title.Color      = 'k';
            ax.Title.FontWeight = 'normal';
            ax.Title.FontSize   = 12;
        end
        if isprop(ax,'XLabel'), ax.XLabel.Color = 'k'; end
        if isprop(ax,'YLabel'), ax.YLabel.Color = 'k'; end
    end
    set(findall(fig,'Type','text'),  'Color','k', 'FontName','Helvetica');
    for lg = findall(fig,'Type','legend')'
        set(lg, 'Box','off', 'TextColor','k');
    end
    set(fig, 'Color','w');
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function y = apply_mean_reversion(y, anchor, max_strength, tau_years)
    y = y(:);
    for k = 1:length(y)
        w = max_strength * (1 - exp(-k / tau_years));
        y(k) = y(k) + w * (anchor - y(k));
    end
end

function [y, y_lower, y_upper] = enforce_ecological_bounds(y, y_lower, y_upper, last_obs, collapse_threshold, resilience_band, decline_damping, min_cover)
    % Keep forecasts ecologically plausible: no rebound after collapse and softer decline near low cover.
    y       = y(:);
    y_lower = y_lower(:);
    y_upper = y_upper(:);
    prev = last_obs;

    for k = 1:length(y)
        cand = y(k);


        if prev <= collapse_threshold || cand <= collapse_threshold
            y(k:end)       = 0;
            y_lower(k:end) = 0;
            y_upper(k:end) = 0;
            break;
        end


        if cand < prev && prev < resilience_band
            cand = prev + decline_damping * (cand - prev);
        end

        cand = max(min_cover, min(100, cand));
        y(k) = cand;

        lo_k = max(0, min(100, y_lower(k)));
        hi_k = max(0, min(100, y_upper(k)));
        y_lower(k) = min(lo_k, y(k));
        y_upper(k) = max(hi_k, y(k));

        prev = y(k);
    end
end

scriptPath = mfilename('fullpath');
if isempty(scriptPath) && usejava('desktop')
    scriptPath = matlab.desktop.editor.getActiveFilename;
end
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end

repoRoot = fileparts(fileparts(scriptDir));
if ~isfolder(fullfile(repoRoot, 'MATLAB'))
    repoRoot = fileparts(scriptDir);
end

% Prefer raw CORIS inputs and fall back to processed if raw is unavailable.
corisRoots = {
    fullfile(repoRoot, 'data', 'raw', 'CORIS_DATA'), ...
    fullfile(repoRoot, 'data', 'processed', 'CORIS_DATA')
};

sourceRoot = "";
for r = 1:numel(corisRoots)
    if isfolder(corisRoots{r})
        sourceRoot = string(corisRoots{r});
        break;
    end
end

if strlength(sourceRoot) == 0
    searched = strjoin(corisRoots, ', ');
    error('No CORIS_DATA root found. Searched: %s', searched);
end

files = dir(fullfile(char(sourceRoot), '**', '*Benthic*Cover*.csv'));

if isempty(files)
    error('No benthic cover files found under %s', char(sourceRoot));
end

fprintf('\nUsing CORIS_DATA source: %s\n', char(sourceRoot));
fprintf('%d benthic cover file(s) found\n', length(files));

raw_tables   = {};
all_col_sets = {};

for f = 1:length(files)
    fpath = fullfile(files(f).folder, files(f).name);
    try
        tbl = readtable(fpath, 'VariableNamingRule','preserve', 'TextType','string');
        if height(tbl) > 0 && ismember('YEAR', tbl.Properties.VariableNames)
            if isnan(str2double(tbl.YEAR(1)))
                tbl(1,:) = [];
            end
        end
        raw_tables{end+1}   = tbl;
        all_col_sets{end+1} = tbl.Properties.VariableNames;
        fprintf('  %-50s  %d rows\n', files(f).name, height(tbl));
    catch ME
        fprintf('  %-50s  skipped (%s)\n', files(f).name, ME.message);
    end
end

if isempty(raw_tables)
    error('No files loaded.');
end


all_cols = unique([all_col_sets{:}], 'stable');
for f = 1:length(raw_tables)
    for col = setdiff(all_cols, raw_tables{f}.Properties.VariableNames)
        raw_tables{f}.(col{1}) = repmat("", height(raw_tables{f}), 1);
    end
    raw_tables{f} = raw_tables{f}(:, all_cols);
end

% Align columns before concatenation so every source file contributes consistently.
data_all  = vertcat(raw_tables{:});
rows_before_dedup = height(data_all);
data_all  = unique(data_all, 'rows', 'stable');
rows_dropped = rows_before_dedup - height(data_all);
col_names = data_all.Properties.VariableNames;
fprintf('\n%d rows total,  %d columns\n', height(data_all), width(data_all));
if rows_dropped > 0
    fprintf('Removed %d exact duplicate rows\n', rows_dropped);
end

fprintf('\nColumns:\n');
for c = 1:length(col_names)
    fprintf('  %2d.  %s\n', c, col_names{c});
end


num_cols = {'YEAR','HARDBOTTOM_P','SOFTBOTTOM_P','RUBBLE_P','latitude', ...
            'longitude','MIN_DEPTH','MAX_DEPTH','PROT','WTD_RUG','STATION_NR'};
n_data = height(data_all);
for nc = 1:length(num_cols)
    col = num_cols{nc};
    if ~ismember(col, col_names), continue; end
    v = data_all.(col);
    if isnumeric(v), continue; end
    v = v(:);
    if numel(v) ~= n_data, v = repmat("", n_data, 1); end
    v(ismissing(v)) = "";
    data_all.(col) = str2double(v);
end


hard_coral_codes = [ ...
    "ACR CERV","ACR PALM","AGA AGAR","AGA FRAG","AGA HUMI","AGA LAMA","AGA SPE.", ...
    "COL NATA","COL SPE.","DEN CYLI","DIC STOK","DIP LABY","DIP SPE.","EUS FAST", ...
    "FAV FRAG","HEL CUCU","ISO SINU","MAD AURE","MAD CARM","MAD DECA","MAD FORM", ...
    "MAD SPE.","MAN AREO","MEA MEAN","MEA SPE.","MIL SPE.","MON CAVE","MUS ANGU", ...
    "MYC ALIC","MYC FERO","MYC LAMA","OCU DIFF","OCU ROBU","ORB ANNU","ORB FAVE", ...
    "ORB FRAN","OTH CORA","POR ASTR","POR DIVA","POR FURC","POR PORI","POR SPE.", ...
    "PSE CLIV","PSE STRI","SCO CUBE","SCO SPE.","SID RADI","SID SIDE","SID SPE.", ...
    "SOL BOUR","SOL SPE.","STE INTE"];

algae_codes = [ ...
    "DIC SPE.","HAL SPE.","LOB SPE.","MAC FLES","MAC CALC","RHO CRUS", ...
    "CYA SPE.","TUR SEDI","TUR FREE","MAG SPE.","PEY SPE.","RAM SPE."];

sponge_codes = ["SPO OTHE","CLI SPE."];


data_all  = data_all(~isnan(data_all.YEAR), :);
col_names = data_all.Properties.VariableNames;
n_rows    = height(data_all);

if ismember('PRIMARY_SAMPLE_UNIT', col_names)
    psu_col = data_all.PRIMARY_SAMPLE_UNIT;
    psu_col(ismissing(psu_col)) = "UNKNOWN";
else
    psu_col = repmat("0", n_rows, 1);
end

if ismember('STATION_NR', col_names)
    stn_raw = data_all.STATION_NR;
    stn_raw(isnan(stn_raw)) = -1;
    stn_col = string(stn_raw);
else
    stn_col = repmat("0", n_rows, 1);
end

year_col_num = data_all.YEAR;
transect_id  = string(year_col_num) + "_" + psu_col + "_" + stn_col;

bad = contains(transect_id,"NaN") | ismissing(transect_id);
if any(bad)
    fprintf('Dropping %d rows with malformed transect IDs\n', sum(bad));
    data_all     = data_all(~bad,:);
    transect_id  = transect_id(~bad);
    psu_col      = psu_col(~bad);
    stn_col      = stn_col(~bad);
    year_col_num = year_col_num(~bad);
end

unique_transects = unique(transect_id);
n_transects      = length(unique_transects);
fprintf('\n%d unique transects\n', n_transects);

T = struct( ...
    'year',        zeros(n_transects,1), ...
    'psu',         strings(n_transects,1), ...
    'station',     strings(n_transects,1), ...
    'lat',         NaN(n_transects,1), ...
    'lon',         NaN(n_transects,1), ...
    'hard_coral',  zeros(n_transects,1), ...
    'algae',       zeros(n_transects,1), ...
    'sponge',      zeros(n_transects,1), ...
    'bare',        zeros(n_transects,1), ...
    'n_coral_spp', zeros(n_transects,1), ...
    'depth_m',     NaN(n_transects,1), ...
    'rugosity',    NaN(n_transects,1), ...
    'protection',  NaN(n_transects,1), ...
    'sub_region',  strings(n_transects,1), ...
    'habitat',     strings(n_transects,1), ...
    'zone',        strings(n_transects,1));


if ismember('HARDBOTTOM_P', data_all.Properties.VariableNames)
    if ~isnumeric(data_all.HARDBOTTOM_P)
        data_all.HARDBOTTOM_P = str2double(data_all.HARDBOTTOM_P);
    end
end


if ismember('COVER_CAT_CD', data_all.Properties.VariableNames)
    data_all.COVER_CAT_CD = strtrim(data_all.COVER_CAT_CD);
end

hb_col  = data_all.HARDBOTTOM_P;
cat_col = data_all.COVER_CAT_CD;


unique_cats = unique(cat_col);
fprintf('\nUnique COVER_CAT_CD values (%d total, first 50 shown):\n', length(unique_cats));
disp(unique_cats(1:min(50,end)));


n_coral_rows = sum(ismember(cat_col, hard_coral_codes));
fprintf('Rows matching hard_coral_codes: %d / %d\n', n_coral_rows, length(cat_col));

for i = 1:n_transects
    rows = find(transect_id == unique_transects(i));
    if isempty(rows), continue; end

    T.year(i)    = year_col_num(rows(1));
    T.psu(i)     = psu_col(rows(1));
    T.station(i) = stn_col(rows(1));

    if ismember('latitude',        col_names), T.lat(i)        = data_all.latitude(rows(1));  end
    if ismember('longitude',       col_names), T.lon(i)        = data_all.longitude(rows(1)); end
    if ismember('MIN_DEPTH',       col_names) && ismember('MAX_DEPTH',col_names)
        T.depth_m(i) = (data_all.MIN_DEPTH(rows(1)) + data_all.MAX_DEPTH(rows(1))) / 2;
    end
    if ismember('WTD_RUG',         col_names), T.rugosity(i)   = mean(data_all.WTD_RUG(rows),'omitnan');  end
    if ismember('PROT',            col_names), T.protection(i) = data_all.PROT(rows(1));                  end
    if ismember('SUB_REGION_NAME', col_names), T.sub_region(i) = data_all.SUB_REGION_NAME(rows(1));       end
    if ismember('HABITAT_CD',      col_names), T.habitat(i)    = data_all.HABITAT_CD(rows(1));            end
    if ismember('ZONE_NAME',       col_names), T.zone(i)       = data_all.ZONE_NAME(rows(1));             end

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
fprintf('%d transects with valid depth\n', height(transect_tbl));

years_all = unique(transect_tbl.year);
years_all = years_all(~isnan(years_all));


annual_mean_hc  = arrayfun(@(y) mean(transect_tbl.hard_coral(transect_tbl.year==y),'omitnan'), years_all);
annual_se_hc    = arrayfun(@(y) std( transect_tbl.hard_coral(transect_tbl.year==y),'omitnan') / ...
                                sqrt(sum(transect_tbl.year==y)), years_all);
annual_mean_alg = arrayfun(@(y) mean(transect_tbl.algae( transect_tbl.year==y),'omitnan'), years_all);
annual_mean_spo = arrayfun(@(y) mean(transect_tbl.sponge(transect_tbl.year==y),'omitnan'), years_all);
annual_n        = arrayfun(@(y) sum( transect_tbl.year==y), years_all);

fprintf('\n%-6s  %-6s  %-10s  %-10s  %-10s  %-10s\n', ...
    'Year','n','HardCoral%','Algae%','Sponge%','CoralSpp');
fprintf('%s\n', repmat('-',1,60));
for y = 1:length(years_all)
    idx = transect_tbl.year == years_all(y);
    fprintf('%-6d  %-6d  %-10.2f  %-10.2f  %-10.2f  %-10.2f\n', ...
        years_all(y), sum(idx), ...
        mean(transect_tbl.hard_coral(idx),'omitnan'), ...
        mean(transect_tbl.algae(idx),'omitnan'), ...
        mean(transect_tbl.sponge(idx),'omitnan'), ...
        mean(transect_tbl.n_coral_spp(idx),'omitnan'));
end

if ~SHOW_ONLY_CORE_FIGURES


    fig1 = figure('Position',[50 50 1400 800], 'Name','Exploratory');

    ax1 = subplot(2,2,1);
    errorbar(years_all, annual_mean_hc, annual_se_hc, 'o-', ...
        'LineWidth',1.8, 'MarkerSize',7, 'CapSize',4, ...
        'Color',CLR_CORAL, 'MarkerFaceColor',CLR_CORAL);
    xlabel('Year'); ylabel('Hard Coral Cover (\%)');
    title('Hard Coral Cover');
    text(years_all, annual_mean_hc + annual_se_hc + 0.3, ...
        arrayfun(@(n) sprintf('n=%d',n), annual_n, 'UniformOutput',false), ...
        'HorizontalAlignment','center', 'FontSize',8, 'Color','k');

    ax2 = subplot(2,2,2);
    plot(years_all, annual_mean_alg, 's-', 'LineWidth',1.8, 'MarkerSize',7, ...
        'Color',CLR_ALGAE, 'MarkerFaceColor',CLR_ALGAE, 'DisplayName','Algae');
    hold on;
    plot(years_all, annual_mean_spo, 'd-', 'LineWidth',1.8, 'MarkerSize',7, ...
        'Color',CLR_SPONGE, 'MarkerFaceColor',CLR_SPONGE, 'DisplayName','Sponge');
    xlabel('Year'); ylabel('Cover (\%)');
    title('Algae and Sponge Cover');
    legend('Location','best');

    ax3 = subplot(2,2,3);
    scatter(transect_tbl.algae, transect_tbl.hard_coral, 18, transect_tbl.year, ...
        'filled', 'MarkerFaceAlpha',0.4);
    colormap(ax3, parula);
    cb = colorbar; cb.Label.String = 'Year'; cb.Color = 'k';
    hold on;

    if std(transect_tbl.algae) > 0 && std(transect_tbl.hard_coral) > 0
        xr = linspace(0, max(transect_tbl.algae), 100);
        plot(xr, polyval(polyfit(transect_tbl.algae, transect_tbl.hard_coral,1), xr), ...
            '--', 'Color',CLR_FCST, 'LineWidth',1.5);
        r_corr = corr(transect_tbl.algae, transect_tbl.hard_coral, 'rows','complete');
        text(0.05, 0.92, sprintf('r = %.3f', r_corr), 'Units','normalized', ...
            'FontSize',10, 'Color','k');
    else
        text(0.05, 0.92, 'r = N/A (zero variance)', 'Units','normalized', ...
            'FontSize',10, 'Color','k');
    end
    xlabel('Algae Cover (\%)'); ylabel('Hard Coral Cover (\%)');
    title('Coral vs. Algae');

    ax4 = subplot(2,2,4);
    scatter(transect_tbl.depth_m, transect_tbl.hard_coral, 18, transect_tbl.year, ...
        'filled', 'MarkerFaceAlpha',0.4);
    colormap(ax4, parula); colorbar;
    xlabel('Depth (m)'); ylabel('Hard Coral Cover (\%)');
    title('Coral Cover vs. Depth');

    sg1 = sgtitle('Benthic Cover - Exploratory Overview');
    sg1.FontSize = 13; sg1.Color = 'k'; sg1.FontWeight = 'normal';
    cleanFig(fig1);
end


fprintf('\nPanel regression\n%s\n', repmat('-',1,50));

transect_tbl.year_norm = transect_tbl.year - min(years_all);
transect_tbl.depth_z   = (transect_tbl.depth_m - mean(transect_tbl.depth_m,'omitnan')) ./ ...
                           std(transect_tbl.depth_m,'omitnan');
transect_tbl.algae_z   = (transect_tbl.algae - mean(transect_tbl.algae,'omitnan')) ./ ...
                           std(transect_tbl.algae,'omitnan');
transect_tbl.prot_z    = transect_tbl.protection;

use_idx = ~isnan(transect_tbl.depth_z)   & ~isnan(transect_tbl.algae_z) & ...
          ~isnan(transect_tbl.hard_coral) & ~isnan(transect_tbl.prot_z);

y_panel = transect_tbl.hard_coral(use_idx);
X_panel = [ones(sum(use_idx),1), ...
           transect_tbl.year_norm(use_idx), ...
           transect_tbl.depth_z(use_idx), ...
           transect_tbl.algae_z(use_idx), ...
           transect_tbl.prot_z(use_idx)];

[b_ols, ~, ~, ~, stats_ols] = regress(y_panel, X_panel);

fprintf('N = %d   R2 = %.4f   F = %.2f   p = %.4f\n', ...
    sum(use_idx), stats_ols(1), stats_ols(2), stats_ols(3));

coef_names = {'Intercept','Year (per yr)','Depth z','Algae z','Protection'};
for c = 1:length(b_ols)
    fprintf('  %-20s  %8.4f\n', coef_names{c}, b_ols(c));
end
fprintf('Adjusted trend: %.4f%%/yr\n', b_ols(2));


fprintf('\nACF / PACF\n%s\n', repmat('-',1,50));

y_annual       = annual_mean_hc;
y_algae_annual = annual_mean_alg;
t_ann          = (1:length(y_annual))';

p_trend     = polyfit(t_ann, y_annual, 1);
trend_cmp   = polyval(p_trend, t_ann);
y_detrended = y_annual - trend_cmp;

fprintf('Linear trend: %.4f%%/yr\n', p_trend(1));

fig2 = figure('Position',[100 100 1200 700], 'Name','ACF / PACF');

    subplot(3,1,1);
    plot(years_all, y_annual, 'o-', 'Color',CLR_CORAL, 'LineWidth',1.8, ...
        'MarkerSize',7, 'MarkerFaceColor',CLR_CORAL, 'DisplayName','Observed');
    hold on;
    plot(years_all, trend_cmp, '--', 'Color',CLR_TREND, 'LineWidth',1.4, ...
        'DisplayName', sprintf('Trend (%.4f%%/yr)', p_trend(1)));
    xlabel('Year'); ylabel('Hard Coral Cover (\%)');
    title('Annual Mean - Hard Coral with Trend');
    legend('Location','best');

    subplot(3,1,2);
    if length(y_detrended) >= 4
        autocorr(y_detrended);
        title('ACF - Detrended Annual Series');
    else
        axis off; title('ACF (< 4 observations)');
    end

    subplot(3,1,3);
    if length(y_detrended) >= 4
        parcorr(y_detrended);
        title('PACF - Detrended Annual Series');
    else
        axis off; title('PACF (< 4 observations)');
    end

sg2 = sgtitle('ACF / PACF - Annual Hard Coral (Detrended)');
sg2.FontSize = 13; sg2.Color = 'k'; sg2.FontWeight = 'normal';
cleanFig(fig2);


fprintf('\nARIMA / SARIMA model selection\n%s\n', repmat('-',1,50));

target_forecast_year = 2100;
forecast_horizon     = target_forecast_year - years_all(end);
if forecast_horizon < 1
    error('Latest observed year (%d) is already >= target year (%d).', years_all(end), target_forecast_year);
end

collapse_threshold    = 0.05;
resilience_band       = 8.00;
low_cover_damping     = 0.15;
min_cover_precollapse = 1.00;

trend_decay_tau       = 18;
reversion_tau         = 10;
reversion_strength    = 0.22;
anchor_window         = min(7, length(y_annual));
coral_anchor          = mean(y_annual(end-anchor_window+1:end), 'omitnan');

model_defs = { ...
    arima('ARLags',1,     'D',0),              'AR(1)';
    arima('MALags',1,     'D',0),              'MA(1)';
    arima('ARLags',1,'MALags',1,'D',0),        'ARMA(1,1)';
    arima('ARLags',[1 2], 'D',0),              'AR(2)';
    arima('MALags',[1 2], 'D',0),              'MA(2)';
    arima('ARLags',[1 2],'MALags',1,'D',0),    'ARMA(2,1)';
    arima('ARLags',1,'MALags',[1 2],'D',0),    'ARMA(1,2)';
    arima('ARLags',[1 2],'MALags',[1 2],'D',0),'ARMA(2,2)';
    arima('D',0),                              'White Noise';
    arima('ARLags',1,'MALags',1,'SARLags',4,'SMALags',4,'D',0), 'SARMA(1,1)x(1,1,4)';
    arima('ARLags',1,'MALags',1,'SARLags',5,'SMALags',5,'D',0), 'SARMA(1,1)x(1,1,5)';
    arima('ARLags',1,'SARLags',4,'D',0),       'SAR(1)x(1,4)';
    arima('ARLags',1,'SARLags',5,'D',0),       'SAR(1)x(1,5)'};

n_models         = size(model_defs,1);
model_candidates = cell(n_models,1);
aic_vals         = NaN(n_models,1);
bic_vals         = NaN(n_models,1);
model_names_list = model_defs(:,2);

fprintf('  %-32s  %10s  %10s\n', 'Model', 'AIC', 'BIC');
fprintf('  %s\n', repmat('-',1,56));

for m = 1:n_models
    try
        fit_m = estimate(model_defs{m,1}, y_detrended, 'Display','off');
        s_m   = summarize(fit_m);
        model_candidates{m} = fit_m;
        aic_vals(m)         = s_m.AIC;
        bic_vals(m)         = s_m.BIC;
        fprintf('  %-32s  %10.3f  %10.3f\n', model_names_list{m}, aic_vals(m), bic_vals(m));
    catch
        fprintf('  %-32s  did not converge\n', model_names_list{m});
    end
end

valid_mask = ~isnan(aic_vals) & ~cellfun(@isempty, model_candidates);
if ~any(valid_mask)
    error('No models converged.');
end

valid_idx       = find(valid_mask);
[~, rel_best]   = min(aic_vals(valid_mask));
best_model_idx  = valid_idx(rel_best);
best_fit        = model_candidates{best_model_idx};
best_model_name = model_names_list{best_model_idx};

fprintf('\n  Selected: %s   AIC = %.3f   BIC = %.3f\n', ...
    best_model_name, aic_vals(best_model_idx), bic_vals(best_model_idx));


[yF_det, yMSE] = forecast(best_fit, forecast_horizon, 'Y0', y_detrended);

h = (1:forecast_horizon)';
trend_future = trend_cmp(end) + p_trend(1) * trend_decay_tau * (1 - exp(-h / trend_decay_tau));

yF       = max(0, min(100, yF_det + trend_future));
yF_upper = min(100, max(0, yF_det + trend_future + 1.96*sqrt(yMSE)));
yF_lower = max(0,         yF_det + trend_future - 1.96*sqrt(yMSE));
yF = apply_mean_reversion(yF, coral_anchor, reversion_strength, reversion_tau);
yF_upper = apply_mean_reversion(yF_upper, coral_anchor, reversion_strength, reversion_tau);
yF_lower = apply_mean_reversion(yF_lower, coral_anchor, reversion_strength, reversion_tau);
[yF, yF_lower, yF_upper] = enforce_ecological_bounds( ...
    yF, yF_lower, yF_upper, y_annual(end), ...
    collapse_threshold, resilience_band, low_cover_damping, min_cover_precollapse);

last_year        = years_all(end);
future_years_vec = (last_year+1 : last_year+forecast_horizon);

fprintf('\n  %-6s  %-10s  %-10s  %-10s\n', 'Year','Forecast','Lower 95','Upper 95');
fprintf('  %s\n', repmat('-',1,42));
for f = 1:forecast_horizon
    fprintf('  %-6d  %-10.2f  %-10.2f  %-10.2f\n', ...
        future_years_vec(f), yF(f), yF_lower(f), yF_upper(f));
end


try
    res_arima    = double(infer(best_fit, y_detrended));
    max_lag      = min(10, floor(length(res_arima)/3));
    [h_lb, p_lb] = lbqtest(res_arima, 'Lags', max_lag);
    fprintf('\n  Ljung-Box (%d lags): p = %.4f  [%s]\n', max_lag, p_lb, ...
        ternary(h_lb==0, 'pass', 'residual autocorrelation remains'));
catch ME
    fprintf('  Ljung-Box failed: %s\n', ME.message);
    res_arima = [];
    h_lb = NaN; p_lb = NaN;
end


fprintf('\nSARIMAX (algae exogenous)\n%s\n', repmat('-',1,50));

results_sarimax = struct('success', false);

try
    X_alg = (y_algae_annual - mean(y_algae_annual)) / std(y_algae_annual);

    sarimax_defs = { ...
        regARIMA('ARLags',1),  'ARIMAX + AR(1)';
        regARIMA('MALags',1),  'ARIMAX + MA(1)'};

    if length(y_annual) >= 8
        sarimax_defs(end+1,:) = {regARIMA('ARLags',1,'SARLags',4), 'SARIMAX + SAR(1,4)'};
        sarimax_defs(end+1,:) = {regARIMA('ARLags',1,'SARLags',5), 'SARIMAX + SAR(1,5)'};
    end

    best_aic_sx  = Inf;
    best_sarimax = [];
    best_sx_name = '';

    fprintf('  %-30s  %10s\n', 'Model', 'AIC');
    fprintf('  %s\n', repmat('-',1,44));

    for sx = 1:size(sarimax_defs,1)
        try
            fit_sx = estimate(sarimax_defs{sx,1}, y_annual, 'X',X_alg, 'Display','off');
            aic_sx = summarize(fit_sx).AIC;
            fprintf('  %-30s  %10.3f\n', sarimax_defs{sx,2}, aic_sx);
            if aic_sx < best_aic_sx
                best_aic_sx  = aic_sx;
                best_sarimax = fit_sx;
                best_sx_name = sarimax_defs{sx,2};
            end
        catch
            fprintf('  %-30s  did not converge\n', sarimax_defs{sx,2});
        end
    end

    if isempty(best_sarimax), error('No SARIMAX model converged.'); end
    fprintf('\n  Selected: %s   AIC = %.3f\n', best_sx_name, best_aic_sx);
    fprintf('  Algae coefficient: %.4f\n', best_sarimax.Beta);

    p_alg            = polyfit((1:length(y_algae_annual))', y_algae_annual, 1);
    future_t_alg     = (length(y_algae_annual)+1 : length(y_algae_annual)+forecast_horizon)';
    alg_forecast_raw = polyval(p_alg, future_t_alg);
    alg_forecast_z   = (alg_forecast_raw - mean(y_algae_annual)) / std(y_algae_annual);

    yF_sx      = forecast(best_sarimax, forecast_horizon, ...
        'Y0',y_annual, 'X0',X_alg, 'XF',alg_forecast_z);
    res_sx     = double(infer(best_sarimax, y_annual, 'X',X_alg));
    res_std_sx = std(res_sx);

    results_sarimax.forecast       = max(0, min(100, yF_sx));
    results_sarimax.forecast_upper = min(100, max(0, yF_sx + 1.96*res_std_sx));
    results_sarimax.forecast_lower = max(0,         yF_sx - 1.96*res_std_sx);
    results_sarimax.forecast = apply_mean_reversion( ...
        results_sarimax.forecast, coral_anchor, reversion_strength, reversion_tau);
    results_sarimax.forecast_upper = apply_mean_reversion( ...
        results_sarimax.forecast_upper, coral_anchor, reversion_strength, reversion_tau);
    results_sarimax.forecast_lower = apply_mean_reversion( ...
        results_sarimax.forecast_lower, coral_anchor, reversion_strength, reversion_tau);
    [results_sarimax.forecast, results_sarimax.forecast_lower, results_sarimax.forecast_upper] = ...
        enforce_ecological_bounds(results_sarimax.forecast, ...
        results_sarimax.forecast_lower, results_sarimax.forecast_upper, y_annual(end), ...
        collapse_threshold, resilience_band, low_cover_damping, min_cover_precollapse);
    results_sarimax.algae_forecast = alg_forecast_raw;
    results_sarimax.model_name     = best_sx_name;
    results_sarimax.beta_algae     = best_sarimax.Beta;
    results_sarimax.success        = true;

    fprintf('\n  %-6s  %-10s  %-10s  %-10s  %-10s\n', ...
        'Year','Forecast','Lower 95','Upper 95','Change');
    for f = 1:forecast_horizon
        fprintf('  %-6d  %-10.2f  %-10.2f  %-10.2f  %+.2f\n', ...
            future_years_vec(f), ...
            results_sarimax.forecast(f), ...
            results_sarimax.forecast_lower(f), ...
            results_sarimax.forecast_upper(f), ...
            results_sarimax.forecast(f) - y_annual(end));
    end

catch ME
    fprintf('  SARIMAX failed: %s\n  Falling back to ARIMA forecast.\n', ME.message);
    results_sarimax.forecast   = yF;
    results_sarimax.model_name = sprintf('ARIMA fallback (%s)', best_model_name);
end


fig3 = figure('Position',[50 50 1400 600], 'Name','Forecast');

axf1 = subplot(1,2,1);
hold on;
fill([future_years_vec, fliplr(future_years_vec)], [yF_upper', fliplr(yF_lower')], ...
    CLR_FCST, 'FaceAlpha',0.12, 'EdgeColor','none', 'DisplayName','95% CI');
trend_future_curve = trend_cmp(end) + p_trend(1) * trend_decay_tau * ...
    (1 - exp(-(1:forecast_horizon)' / trend_decay_tau));
plot([years_all(:); future_years_vec(:)], [trend_cmp(:); trend_future_curve(:)], ...
    '--', 'Color',CLR_TREND, 'LineWidth',1.2, ...
    'DisplayName', sprintf('Damped trend (%.4f\\%%/yr initial)', p_trend(1)));
plot(future_years_vec, yF, 'o-', 'LineWidth',1.8, 'MarkerSize',7, ...
    'Color',CLR_FCST, 'MarkerFaceColor',CLR_FCST, ...
    'DisplayName', sprintf('ARIMA (%s)', best_model_name));
if results_sarimax.success
    plot(future_years_vec, results_sarimax.forecast, 's-', 'LineWidth',1.8, 'MarkerSize',7, ...
        'Color',CLR_SARIMAX, 'MarkerFaceColor',CLR_SARIMAX, ...
        'DisplayName', results_sarimax.model_name);
end
errorbar(years_all, annual_mean_hc, annual_se_hc, 'o', 'LineWidth',1.8, ...
    'MarkerSize',7, 'CapSize',4, ...
    'Color',CLR_CORAL, 'MarkerFaceColor',CLR_CORAL, 'DisplayName','Observed +/- SE');
xlabel('Year'); ylabel('Hard Coral Cover (\%)');
title('Hard Coral Cover - Forecast');
legend('Location','best','FontSize',9);
ylim([0, max([annual_mean_hc(:); yF_upper(:)])*1.35 + 1]);
xlim([years_all(1)-1, future_years_vec(end)+0.5]);

axf2 = subplot(1,2,2);
hold on;
yyaxis left;
plot(years_all, annual_mean_alg, 'o-', 'LineWidth',1.8, 'MarkerSize',7, ...
    'Color',CLR_ALGAE, 'MarkerFaceColor',CLR_ALGAE, 'DisplayName','Algae (obs)');
if results_sarimax.success
    plot(future_years_vec, results_sarimax.algae_forecast, '--', ...
        'LineWidth',1.4, 'Color',CLR_ALGAE, 'DisplayName','Algae (proj)');
end
set(axf2, 'YColor','k');
ylabel('Algae Cover (\%)');

yyaxis right;
plot(years_all, annual_mean_hc, 'o-', 'LineWidth',1.8, 'MarkerSize',7, ...
    'Color',CLR_CORAL, 'MarkerFaceColor',CLR_CORAL, 'DisplayName','Coral (obs)');
plot(future_years_vec, yF, '--', 'LineWidth',1.4, 'Color',CLR_FCST, ...
    'DisplayName','Coral (fcst)');
set(axf2, 'YColor','k');
ylabel('Hard Coral Cover (\%)');
xlabel('Year');
title('Algae-Coral Dynamics');
legend('Location','best','FontSize',9);

sg3 = sgtitle('Benthic Cover - ARIMA + SARIMAX Forecasts');
sg3.FontSize = 13; sg3.Color = 'k'; sg3.FontWeight = 'normal';
cleanFig(fig3);


fprintf('\nRegional analysis\n%s\n', repmat('-',1,50));

sub_regions   = unique(transect_tbl.sub_region);
sub_regions   = sub_regions(sub_regions ~= "");
sub_yr_count  = arrayfun(@(s) length(unique(transect_tbl.year(transect_tbl.sub_region==s))), sub_regions);
model_regions = sub_regions(sub_yr_count >= 3);
fprintf('%d regions with >= 3 survey years\n', length(model_regions));

regional_results = struct();

if ~SHOW_ONLY_CORE_FIGURES
    fig4 = figure('Position',[50 50 1600 900], 'Name','Regional Forecasts');
end
n_plot_rows = ceil(length(model_regions)/3);
n_plot_cols = min(3, length(model_regions));

for s = 1:length(model_regions)
    rname = model_regions(s);
    field = matlab.lang.makeValidName(rname);

    sub_data  = transect_tbl(transect_tbl.sub_region == rname, :);
    sub_years = unique(sub_data.year)';
    y_reg  = arrayfun(@(y) mean(sub_data.hard_coral(sub_data.year==y),'omitnan'), sub_years);
    se_reg = arrayfun(@(y) std( sub_data.hard_coral(sub_data.year==y),'omitnan') / ...
             sqrt(sum(sub_data.year==y)), sub_years);
    n_reg  = arrayfun(@(y) sum(sub_data.year==y), sub_years);

    t_reg     = (1:length(y_reg))';
    p_reg     = polyfit(t_reg, y_reg, 1);
    y_det_reg = y_reg - polyval(p_reg, t_reg);

    fprintf('  %-35s  %d yrs  mean %.2f%%  trend %.4f%%/yr\n', ...
        strtrim(rname), length(sub_years), mean(y_reg), p_reg(1));

    try
        if length(y_reg) >= 8
            reg_defs = { ...
                arima('ARLags',1,'MALags',1,'D',0),                            'ARMA(1,1)';
                arima('ARLags',1,'D',0),                                     'AR(1)';
                arima('MALags',1,'D',0),                                     'MA(1)';
                arima('ARLags',1,'MALags',1,'SARLags',4,'SMALags',4,'D',0), 'SARMA x(4)'};
        elseif length(y_reg) >= 5
            reg_defs = { ...
                arima('ARLags',1,'MALags',1,'D',0),  'ARMA(1,1)';
                arima('ARLags',1,'D',0), 'AR(1)';
                arima('MALags',1,'D',0), 'MA(1)'};
        elseif length(y_reg) >= 4
            reg_defs = { ...
                arima('ARLags',1,'D',0), 'AR(1)';
                arima('MALags',1,'D',0), 'MA(1)'};
        else
            reg_defs = {arima('D',0), 'White Noise'};
        end

        best_aic_reg = Inf;
        best_fit_reg = [];
        for rm = 1:size(reg_defs,1)
            try
                fm = estimate(reg_defs{rm,1}, y_det_reg, 'Display','off');
                am = summarize(fm).AIC;
                if am < best_aic_reg
                    best_aic_reg = am;
                    best_fit_reg = fm;
                end
            catch; end
        end
        if isempty(best_fit_reg), error('no model converged'); end

        [yF_det_reg, yMSE_reg] = forecast(best_fit_reg, forecast_horizon, 'Y0',y_det_reg);

        ft_reg        = (length(y_reg)+1 : length(y_reg)+forecast_horizon)';
        trend_ft      = polyval(p_reg, ft_reg);
        yF_reg        = max(0, min(100, yF_det_reg + trend_ft));
        yF_reg_up     = min(100, max(0, yF_det_reg + trend_ft + 1.96*sqrt(yMSE_reg)));
        yF_reg_lo     = max(0,         yF_det_reg + trend_ft - 1.96*sqrt(yMSE_reg));
        [yF_reg, yF_reg_lo, yF_reg_up] = enforce_ecological_bounds( ...
            yF_reg, yF_reg_lo, yF_reg_up, y_reg(end), ...
            collapse_threshold, resilience_band, low_cover_damping, min_cover_precollapse);
        future_yr_reg = sub_years(end)+1 : sub_years(end)+forecast_horizon;

        regional_results.(field).sub_region     = rname;
        regional_results.(field).years          = sub_years;
        regional_results.(field).observed       = y_reg;
        regional_results.(field).se             = se_reg;
        regional_results.(field).trend_slope    = p_reg(1);
        regional_results.(field).forecast       = yF_reg;
        regional_results.(field).forecast_upper = yF_reg_up;
        regional_results.(field).forecast_lower = yF_reg_lo;
        regional_results.(field).forecast_years = future_yr_reg;

        if ~SHOW_ONLY_CORE_FIGURES
            axr = subplot(n_plot_rows, n_plot_cols, s);
            hold on;
            fill([future_yr_reg, fliplr(future_yr_reg)], [yF_reg_up', fliplr(yF_reg_lo')], ...
                CLR_FCST, 'FaceAlpha',0.12, 'EdgeColor','none', 'DisplayName','95% CI');
            yr_all    = [sub_years(:); future_yr_reg(:)];
            trend_all = polyval(p_reg, (1:length(yr_all))');
            plot(yr_all, trend_all, '--', 'Color',CLR_TREND, 'LineWidth',1, 'DisplayName','Trend');
            errorbar(sub_years, y_reg, se_reg, 'o-', 'LineWidth',1.6, 'MarkerSize',6, 'CapSize',3, ...
                'Color',CLR_CORAL, 'MarkerFaceColor',CLR_CORAL, 'DisplayName','Observed');
            plot(future_yr_reg, yF_reg, 'o-', 'LineWidth',1.6, 'MarkerSize',6, ...
                'Color',CLR_FCST, 'MarkerFaceColor',CLR_FCST, 'DisplayName','Forecast');
            xlabel('Year'); ylabel('Hard Coral (\%)');
            title(sprintf('%s  (n=%d)', strtrim(rname), round(mean(n_reg))));
            lg_r = legend('Location','best','FontSize',7);
            lg_r.Box = 'off'; lg_r.TextColor = 'k';
            xlim([min(sub_years)-1, max(future_yr_reg)+0.5]);
            ylim([0, max([y_reg(:); yF_reg_up(:)])*1.25 + 2]);
        end

    catch ME
        fprintf('    %s: %s\n', strtrim(rname), ME.message);
    end
end

if ~SHOW_ONLY_CORE_FIGURES
    sg4 = sgtitle('Hard Coral Cover - Regional Forecasts');
    sg4.FontSize = 13; sg4.Color = 'k'; sg4.FontWeight = 'normal';
    cleanFig(fig4);
end

if SHOW_SPECIES_FIGURE


    top_species = ["ORB FAVE","SID SIDE","MON CAVE","POR ASTR","PSE STRI", ...
                   "COL NATA","STE INTE","ACR CERV","MIL SPE.","SID RADI"];
    top_species_names = {'O. faveolata','S. siderea','M. cavernosa','P. astreoides', ...
                         'P. strigosa','C. natans','S. intersepta','A. cervicornis', ...
                         'Millepora spp.','S. radians'};

    spp_annual = zeros(length(top_species), length(years_all));


    if ~isnumeric(data_all.HARDBOTTOM_P)
        data_all.HARDBOTTOM_P = str2double(data_all.HARDBOTTOM_P);
    end
    data_all.COVER_CAT_CD = strtrim(data_all.COVER_CAT_CD);

    for sp = 1:length(top_species)
        sp_rows = (data_all.COVER_CAT_CD == top_species(sp));
        for y = 1:length(years_all)
            mask = sp_rows & (data_all.YEAR == years_all(y));
            if any(mask)
                spp_annual(sp,y) = mean(data_all.HARDBOTTOM_P(mask),'omitnan');
            end
        end
    end

    colors_spp = [0.00 0.45 0.70; 0.90 0.62 0.00; 0.00 0.62 0.45; 0.80 0.14 0.15; ...
                  0.49 0.18 0.56; 0.34 0.71 0.91; 0.94 0.89 0.26; 0.10 0.10 0.10; ...
                  0.65 0.34 0.16; 0.56 0.56 0.56];

    fig5 = figure('Position',[50 50 1400 580], 'Name','Species Composition');

    subplot(1,2,1);
    area(years_all, spp_annual');
    legend(top_species_names, 'Location','eastoutside', 'FontSize',8);
    xlabel('Year'); ylabel('Mean Cover (pts)');
    title('Composition - Stacked');
    xlim([years_all(1)-0.5, years_all(end)+0.5]);

    subplot(1,2,2);
    hold on;
    for sp = 1:length(top_species)
        if any(spp_annual(sp,:) > 0)
            plot(years_all, spp_annual(sp,:), 'o-', 'LineWidth',1.6, 'MarkerSize',6, ...
                'Color',colors_spp(sp,:), 'MarkerFaceColor',colors_spp(sp,:), ...
                'DisplayName',top_species_names{sp});
        end
    end
    xlabel('Year'); ylabel('Mean Cover (pts)');
    title('Per-Species Trends');
    legend('Location','eastoutside','FontSize',8);
    xlim([years_all(1)-0.5, years_all(end)+0.5]);

    sg5 = sgtitle('Hard Coral Species Composition');
    sg5.FontSize = 13; sg5.Color = 'k'; sg5.FontWeight = 'normal';
    cleanFig(fig5);
end

if ~SHOW_ONLY_CORE_FIGURES


    y_hat     = X_panel * b_ols;
    res_panel = y_panel - y_hat;

    fig6 = figure('Position',[50 50 1400 720], 'Name','Residual Diagnostics');

    axd1 = subplot(2,3,1);
    scatter(y_hat, res_panel, 14, CLR_TREND, 'filled', 'MarkerFaceAlpha',0.3);
    yline(0, '--', 'Color',CLR_FCST, 'LineWidth',1.4);
    xlabel('Fitted (\%)'); ylabel('Residuals (\%)');
    title('Panel OLS - Residuals vs. Fitted');

    axd2 = subplot(2,3,2);
    histogram(res_panel, 40, 'Normalization','pdf', 'FaceColor',CLR_LGRAY, 'EdgeColor','none');
    hold on;
    xr = linspace(min(res_panel), max(res_panel), 100);
    plot(xr, normpdf(xr, mean(res_panel), std(res_panel)), '-', 'Color',CLR_FCST, 'LineWidth',1.8);
    xlabel('Residuals (\%)'); ylabel('Density');
    title('Panel OLS - Residual Distribution');
    legend('Data','Normal','Location','best');

    axd3 = subplot(2,3,3);
    qqplot(res_panel);
    title('Panel OLS - Q-Q Plot');
    try
        axd3.Children(1).Color = CLR_CORAL;
        axd3.Children(2).Color = CLR_TREND;
    catch; end

    if ~isempty(res_arima)
        axd4 = subplot(2,3,4);
        stem(years_all, res_arima, 'Color',CLR_CORAL, 'LineWidth',1.4, ...
            'MarkerFaceColor',CLR_CORAL, 'MarkerSize',5);
        yline(0, '--', 'Color',CLR_TREND, 'LineWidth',1);
        xlabel('Year'); ylabel('Residuals (\%)');
        title(sprintf('ARIMA Residuals - %s', best_model_name));

        axd5 = subplot(2,3,5);
        if length(res_arima) >= 4
            autocorr(res_arima);
            title('ARIMA Residual ACF');
        else
            axis off; title('ACF (< 4 obs)');
        end
    else
        subplot(2,3,4); axis off;
        subplot(2,3,5); axis off;
    end

    subplot(2,3,6);
    axis off;
    lb_str = 'n/a';
    if ~isnan(h_lb)
        lb_str = sprintf('p = %.4f  [%s]', p_lb, ternary(h_lb==0,'pass','fail'));
    end
    text(0.05, 0.97, sprintf([ ...
        'Summary\n\n' ...
        'Panel OLS  (N = %d)\n' ...
        '  R2            %.4f\n' ...
        '  mean resid    %.4f%%\n' ...
        '  std  resid    %.4f%%\n\n' ...
        'ARIMA  %s\n' ...
        '  trend         %.4f%%/yr\n' ...
        '  AIC           %.3f\n' ...
        '  Ljung-Box     %s'], ...
        sum(use_idx), stats_ols(1), mean(res_panel), std(res_panel), ...
        best_model_name, p_trend(1), aic_vals(best_model_idx), lb_str), ...
        'Units','normalized', 'VerticalAlignment','top', ...
        'FontSize',9, 'FontName','Courier', 'Color','k');

    sg6 = sgtitle('Residual Diagnostics');
    sg6.FontSize = 13; sg6.Color = 'k'; sg6.FontWeight = 'normal';
    cleanFig(fig6);
end


fprintf('\n%s\n', repmat('=',1,72));
fprintf('Forecast summary\n');
fprintf('%s\n', repmat('=',1,72));
fprintf('  %-6s  %-12s  %-10s  %-10s  %-14s\n', ...
    'Year','ARIMA','Lower 95','Upper 95','SARIMAX');
fprintf('  %s\n', repmat('-',1,60));
for f = 1:forecast_horizon
    sx_val = NaN;
    if results_sarimax.success, sx_val = results_sarimax.forecast(f); end
    fprintf('  %-6d  %-12.2f  %-10.2f  %-10.2f  %-14.2f\n', ...
        future_years_vec(f), yF(f), yF_lower(f), yF_upper(f), sx_val);
end
fprintf('%s\n', repmat('=',1,72));
fprintf('Trend: %.4f%%/yr\n', p_trend(1));
if results_sarimax.success
    fprintf('SARIMAX model: %s   beta(algae) = %.4f\n', ...
        results_sarimax.model_name, results_sarimax.beta_algae);
end

fprintf('\n  %-35s  %-6s  %-10s  %-12s\n', 'Region','N yrs','Trend/yr','Fcst mean%');
fprintf('  %s\n', repmat('-',1,68));
for s = 1:length(model_regions)
    field = matlab.lang.makeValidName(model_regions(s));
    if isfield(regional_results, field)
        fprintf('  %-35s  %-6d  %-10.4f  %-12.2f\n', ...
            strtrim(model_regions(s)), ...
            length(regional_results.(field).years), ...
            regional_results.(field).trend_slope, ...
            mean(regional_results.(field).forecast));
    end
end


fprintf('\nWriting CSVs\n%s\n', repmat('-',1,50));

global_tbl = table(future_years_vec', yF, yF_lower, yF_upper, ...
    'VariableNames',{'Year','ARIMA_Forecast_HardCoral_Pct','Lower_95CI','Upper_95CI'});
if results_sarimax.success
    global_tbl.SARIMAX_Forecast    = results_sarimax.forecast;
    global_tbl.SARIMAX_Lower       = results_sarimax.forecast_lower;
    global_tbl.SARIMAX_Upper       = results_sarimax.forecast_upper;
    global_tbl.Projected_Algae_Pct = results_sarimax.algae_forecast;
end
forecastPath = fullfile(scriptDir, 'benthic_cover_forecast.csv');
writetable(global_tbl, forecastPath);
fprintf('  %s\n', forecastPath);

obs_tbl = table(years_all(:), annual_mean_hc(:), annual_mean_alg(:), annual_mean_spo(:), ...
    'VariableNames', {'Year','Observed_HardCoral_Pct','Observed_Algae_Pct','Observed_Sponge_Pct'});
observedPath = fullfile(scriptDir, 'benthic_cover_observed_annual.csv');
writetable(obs_tbl, observedPath);
fprintf('  %s\n', observedPath);

for s = 1:length(model_regions)
    field = matlab.lang.makeValidName(model_regions(s));
    if isfield(regional_results, field)
        reg_tbl = table(regional_results.(field).forecast_years', ...
            regional_results.(field).forecast, ...
            regional_results.(field).forecast_lower, ...
            regional_results.(field).forecast_upper, ...
            'VariableNames',{'Year','Forecast_HardCoral_Pct','Lower_95CI','Upper_95CI'});
        fname = sprintf('benthic_forecast_%s.csv', ...
            strrep(strtrim(model_regions(s)),' ','_'));
        fpath = fullfile(scriptDir, fname);
        writetable(reg_tbl, fpath);
        fprintf('  %s\n', fpath);
    end
end

writetable(table(coef_names', b_ols, 'VariableNames',{'Predictor','Coefficient'}), ...
    fullfile(scriptDir, 'benthic_panel_regression_coefficients.csv'));
fprintf('  %s\n', fullfile(scriptDir, 'benthic_panel_regression_coefficients.csv'));

figuresDir = fullfile(repoRoot, 'MATLAB', 'Figures');
if ~isfolder(figuresDir)
    mkdir(figuresDir);
end

fprintf('\nVector PDFs written:\n');
pdf_path = fullfile(figuresDir, 'benthic_cover_forecast_main.pdf');
if ~isempty(fig3) && isgraphics(fig3)
    exportgraphics(fig3, pdf_path, 'ContentType','vector');
    fprintf('  %s\n', pdf_path);
else
    warning('Forecast figure (fig3) is not available; no PDF exported.');
end

if p_trend(1) < 0
    fprintf('\nDeclining trend (%.4f%%/yr)\n', p_trend(1));
else
    fprintf('\nTrend: %.4f%%/yr\n', p_trend(1));
end