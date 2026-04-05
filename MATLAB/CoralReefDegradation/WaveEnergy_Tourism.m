clc;
clear;
close all;


C_CORAL  = [0.13 0.39 0.68];
C_WAVE   = [0.20 0.60 0.86];
C_TOUR   = [0.90 0.45 0.00];
C_TOTAL  = [0.80 0.17 0.17];
C_REF    = [0.40 0.40 0.40];
C_FCST   = [0.65 0.65 0.65];
C_CI     = [0.85 0.85 0.85];

FONT     = 'Helvetica';
FS_AX    = 10;
FS_TTL   = 12;
FS_LBL   = 10;
FS_TXT   =  8;
LW_MAIN  = 1.8;
LW_REF   = 1.2;
MS       = 6;


set(groot,'defaultAxesColor',       'w');
set(groot,'defaultFigureColor',     'w');
set(groot,'defaultAxesBox',         'off');
set(groot,'defaultAxesXGrid',       'off');
set(groot,'defaultAxesYGrid',       'off');
set(groot,'defaultAxesFontName',    FONT);
set(groot,'defaultAxesFontSize',    FS_AX);
set(groot,'defaultAxesXColor',      'k');
set(groot,'defaultAxesYColor',      'k');
set(groot,'defaultAxesZColor',      'k');
set(groot,'defaultAxesLineWidth',   0.9);
set(groot,'defaultAxesTickDir',     'out');
set(groot,'defaultTextFontName',    FONT);
set(groot,'defaultTextColor',       'k');
set(groot,'defaultTextInterpreter', 'latex');
set(groot,'defaultLegendInterpreter', 'latex');
set(groot,'defaultAxesTickLabelInterpreter', 'latex');
set(groot,'defaultLegendBox',       'off');
set(groot,'defaultLegendTextColor', 'k');
set(groot,'defaultLegendFontSize',  FS_TXT+1);


function cleanFig(fig)
    set(fig, 'Color','w');

    all_ax = findall(fig, 'Type','axes');
    for ax = all_ax'
        set(ax, ...
            'Color',          'w', ...
            'XColor',         'k', ...
            'YColor',         'k', ...
            'ZColor',         'k', ...
            'Box',            'off', ...
            'XGrid',          'off', ...
            'YGrid',          'off', ...
            'ZGrid',          'off', ...
            'GridColor',      'k', ...
            'MinorGridColor', 'k', ...
            'TickDir',        'out', ...
            'LineWidth',      0.9,  ...
            'FontName',       'Helvetica', ...
            'FontSize',       10);
        ax.Title.Color      = 'k';
        ax.Title.FontName   = 'Helvetica';
        ax.Title.FontSize   = 11;
        ax.Title.FontWeight = 'normal';
        ax.XLabel.Color     = 'k';
        ax.XLabel.FontName  = 'Helvetica';
        ax.YLabel.Color     = 'k';
        ax.YLabel.FontName  = 'Helvetica';
    end

    all_txt = findall(fig, 'Type','text');
    set(all_txt, 'Color','k', 'FontName','Helvetica');

    all_lg = findall(fig, 'Type','legend');
    for lg = all_lg'
        set(lg, 'Box','off', 'TextColor','k', ...
            'Color','none', 'EdgeColor','none');
    end

    all_cb = findall(fig, 'Type','colorbar');
    for cb = all_cb'
        cb.Color     = 'k';
        cb.FontName  = 'Helvetica';
        cb.FontSize  = 9;
        if ~isempty(cb.Label)
            cb.Label.Color   = 'k';
            cb.Label.FontName = 'Helvetica';
        end
    end

    all_fig_txt = findall(fig, 'Type','text');
    set(all_fig_txt, 'Color','k', 'FontName','Helvetica');
end


function styleAx(ax, xl, yl, ttl)
    ax.Box      = 'off';
    ax.XGrid    = 'off';
    ax.YGrid    = 'off';
    ax.Color    = 'w';
    ax.XColor   = 'k';
    ax.YColor   = 'k';
    ax.ZColor   = 'k';
    ax.TickDir  = 'out';
    ax.FontName = 'Helvetica';
    ax.FontSize = 10;
    if nargin >= 2 && ~isempty(xl)
        xlabel(ax, xl, 'Color','k', 'FontName','Helvetica', 'FontSize',10, 'Interpreter','latex');
    end
    if nargin >= 3 && ~isempty(yl)
        ylabel(ax, yl, 'Color','k', 'FontName','Helvetica', 'FontSize',10, 'Interpreter','latex');
    end
    if nargin >= 4 && ~isempty(ttl)
        title(ax, ttl, 'Color','k', 'FontName','Helvetica', ...
            'FontSize',11, 'FontWeight','normal', 'Interpreter','latex');
    end
end

function sg = makeSgtitle(txt)
    sg = sgtitle(txt, 'Interpreter','latex');
    sg.FontName   = 'Helvetica';
    sg.FontSize   = 12;
    sg.FontWeight = 'normal';
    sg.Color      = 'k';
end

function idx = marker_idx_from_years(years, step_years)
    idx = 1:step_years:numel(years);
    if isempty(idx) || idx(end) ~= numel(years)
        idx = [idx numel(years)];
    end
end


scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptPath = which('WaveEnergy_Tourism');
end
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

forecast_candidates = {
    fullfile(scriptDir, 'benthic_cover_forecast.csv'), ...
    fullfile(repoRoot, 'data', 'analyzed', 'benthic_cover_forecast.csv'), ...
    'benthic_cover_forecast.csv'
};
observed_candidates = {
    fullfile(scriptDir, 'benthic_cover_observed_annual.csv'), ...
    fullfile(repoRoot, 'data', 'analyzed', 'benthic_cover_observed_annual.csv'), ...
    'benthic_cover_observed_annual.csv'
};

forecast_tbl = [];
for i = 1:numel(forecast_candidates)
    if isfile(forecast_candidates{i})
        forecast_tbl = readtable(forecast_candidates{i});
        break;
    end
end

observed_tbl = [];
for i = 1:numel(observed_candidates)
    if isfile(observed_candidates{i})
        observed_tbl = readtable(observed_candidates{i});
        break;
    end
end

if ~isempty(forecast_tbl) && ~isempty(observed_tbl)
    obs_years = observed_tbl.Year(:)';
    obs_cover = observed_tbl.Observed_HardCoral_Pct(:)';

    fcst_years = forecast_tbl.Year(:)';
    if ismember('SARIMAX_Forecast', forecast_tbl.Properties.VariableNames)
        fcst_cover = forecast_tbl.SARIMAX_Forecast(:)';
    else
        fcst_cover = forecast_tbl.ARIMA_Forecast_HardCoral_Pct(:)';
    end
    fcst_lo = forecast_tbl.Lower_95CI(:)';
    fcst_hi = forecast_tbl.Upper_95CI(:)';
    source_tag = 'Benthic forecast CSV';

elseif exist('years_all','var')        && exist('annual_mean_hc','var') && ...
       exist('future_years_vec','var') && exist('yF','var')

    obs_years  = years_all(:)';
    obs_cover  = annual_mean_hc(:)';
    fcst_years = future_years_vec(:)';
    fcst_cover = yF(:)';
    fcst_lo    = max(0,   yF_lower(:)');
    fcst_hi    = min(100, yF_upper(:)');
    source_tag = 'CORIS / ARIMA pipeline';

else
    fprintf('[INFO] Using synthetic FL Keys series (pipeline workspace not found).\n\n');
    rng(42);
    obs_years = 2000:2023;
    n_s       = length(obs_years);
    obs_cover = max(1, 28 - 0.55*(0:n_s-1) + 0.9*randn(1,n_s) + ...
                    1.3*sin(2*pi*(0:n_s-1)/4.8));

    fcst_years = 2024:2100;
    n_f        = length(fcst_years);

    fcst_cover = max(1, obs_cover(end) - 0.55*(1:n_f));
    fcst_lo    = max(0,   fcst_cover - 3.5);
    fcst_hi    = min(100, fcst_cover + 3.5);
    source_tag = 'Synthetic (standalone)';
end

target_forecast_year = 2100;
if fcst_years(end) < target_forecast_year
    n_add = target_forecast_year - fcst_years(end);
    yr_add = (fcst_years(end)+1):target_forecast_year;

    if numel(fcst_cover) >= 5
        tail_slope = mean(diff(fcst_cover(end-4:end)));
    elseif numel(fcst_cover) >= 2
        tail_slope = fcst_cover(end) - fcst_cover(end-1);
    else
        tail_slope = 0;
    end

    cov_add = zeros(1, n_add);
    lo_add  = zeros(1, n_add);
    hi_add  = zeros(1, n_add);
    prev_c  = fcst_cover(end);
    prev_lo = fcst_lo(end);
    prev_hi = fcst_hi(end);
    for ii = 1:n_add
        cand = prev_c + tail_slope;
        if prev_c <= 0 || cand <= 0
            cand = 0;
        end
        cov_add(ii) = max(0, min(100, cand));
        lo_add(ii)  = max(0, prev_lo + tail_slope);
        hi_add(ii)  = max(cov_add(ii), min(100, prev_hi + tail_slope));
        prev_c  = cov_add(ii);
        prev_lo = lo_add(ii);
        prev_hi = hi_add(ii);
    end

    fcst_years = [fcst_years, yr_add];
    fcst_cover = [fcst_cover, cov_add];
    fcst_lo    = [fcst_lo, lo_add];
    fcst_hi    = [fcst_hi, hi_add];
elseif fcst_years(end) > target_forecast_year
    keep_idx = fcst_years <= target_forecast_year;
    fcst_years = fcst_years(keep_idx);
    fcst_cover = fcst_cover(keep_idx);
    fcst_lo    = fcst_lo(keep_idx);
    fcst_hi    = fcst_hi(keep_idx);
end

all_years = [obs_years,  fcst_years];
all_cover = [obs_cover,  fcst_cover];
n_obs     = length(obs_years);
n_fcst    = length(fcst_years);


A0         = obs_cover(1);
T0         = 2.1e6;
M0         = 1800;
M_base     = T0 * M0;

E0         = 1.0e6;
k_atten    = 0.004;
A_reef_km2 = 9.0;
A_reef_m2  = A_reef_km2 * 1e6;
C_ref      = A0;
D_max      = 5.0e9;
lambda_dmg = 2.5;

r_disc     = 0.015;
if abs(r_disc - 0.015) > 1e-12
    error('Discount rate must be 1.5%% (r_disc = 0.015).');
end


E_shore   = @(C) E0 .* exp(-k_atten .* C .* (A_reef_m2/1e6));
E_ref_val = E_shore(C_ref);


E_t      = E_shore(all_cover);
dE_frac  = max(0, (E_t - E_ref_val) / E_ref_val);
D_wave   = D_max .* (1 - exp(-lambda_dmg .* dE_frac));


T_t    = T0 .* (all_cover / A0);
M_t    = T_t .* M0;
dM_t   = max(0, M_base - M_t);


L_total = D_wave + dM_t;


D_wave_obs  = D_wave(1:n_obs);    dM_obs  = dM_t(1:n_obs);
D_wave_fcst = D_wave(n_obs+1:end);dM_fcst = dM_t(n_obs+1:end);
T_obs       = T_t(1:n_obs);       T_fcst  = T_t(n_obs+1:end);
L_obs       = L_total(1:n_obs);   L_fcst  = L_total(n_obs+1:end);
dT_obs      = max(0, T0-T_obs);   dT_fcst = max(0, T0-T_fcst);


Dw_from_lo = D_max .* (1 - exp(-lambda_dmg .* ...
              max(0,(E_shore(fcst_lo)-E_ref_val)/E_ref_val)));
Dw_from_hi = D_max .* (1 - exp(-lambda_dmg .* ...
              max(0,(E_shore(fcst_hi)-E_ref_val)/E_ref_val)));


dM_from_lo = max(0, M_base - T0.*(fcst_lo/A0).*M0);
dM_from_hi = max(0, M_base - T0.*(fcst_hi/A0).*M0);


L_lo = min(Dw_from_lo + dM_from_lo, Dw_from_hi + dM_from_hi);
L_hi = max(Dw_from_lo + dM_from_lo, Dw_from_hi + dM_from_hi);


D_wave_lo = min(Dw_from_lo, Dw_from_hi);
D_wave_hi = max(Dw_from_lo, Dw_from_hi);
dM_lo     = min(dM_from_lo, dM_from_hi);
dM_hi     = max(dM_from_lo, dM_from_hi);


t_disc    = (1:n_fcst)';
NPV_wave  = sum(D_wave_fcst(:) ./ (1+r_disc).^t_disc);
NPV_tour  = sum(dM_fcst(:)     ./ (1+r_disc).^t_disc);
NPV_total = NPV_wave + NPV_tour;
NPV_lo    = sum((D_wave_lo(:)+dM_lo(:)) ./ (1+r_disc).^t_disc);
NPV_hi    = sum((D_wave_hi(:)+dM_hi(:)) ./ (1+r_disc).^t_disc);


cover_sq = fcst_cover;

sq_decline_rate = max(0, (cover_sq(1) - cover_sq(end)) / max(1, n_fcst-1));
hold_decline_fraction = 0.35;
hold_decline_rate = hold_decline_fraction * sq_decline_rate;
cover_hold = max(0, obs_cover(end) - hold_decline_rate * (0:n_fcst-1));

restoration_lag_years = 6;
rest_year_idx = 0:n_fcst-1;
rest_progress = max(0, rest_year_idx - restoration_lag_years) ./ ...
    max(1, (n_fcst-1 - restoration_lag_years));
rest_shape = rest_progress .^ 0.7;
cover_rest = obs_cover(end) + (A0 - obs_cover(end)) * rest_shape;
cover_rest = min(100, max(0, cover_rest));

scenarios      = {cover_sq, cover_hold, cover_rest};
scenario_names = {'Status Quo', 'Hold Cover (Managed)', 'Restoration'};
npv_scenarios  = zeros(1,3);
for sc = 1:3
    C_sc  = scenarios{sc};
    dE_sc = max(0, (E_shore(C_sc) - E_ref_val) / E_ref_val);
    Dw_sc = D_max .* (1 - exp(-lambda_dmg .* dE_sc));
    dM_sc = max(0, M_base - T0.*(C_sc/A0).*M0);
    npv_scenarios(sc) = sum((Dw_sc(:)+dM_sc(:))./(1+r_disc).^t_disc)/1e6;
end
npv_saved = npv_scenarios(1) - npv_scenarios;


if any(cover_sq < 0) || any(cover_hold < 0) || any(cover_rest < 0)
    error('Scenario cover values must be non-negative.');
end
if any(cover_sq > 100) || any(cover_hold > 100) || any(cover_rest > 100)
    error('Scenario cover values must stay within 0-100%%.');
end
if npv_scenarios(2) <= npv_scenarios(3)
    warning('Hold-cover scenario is not worse than restoration; verify assumptions.');
end

fprintf('Mitigation assumptions:\n');
fprintf('  Hold-cover decline fraction vs status quo: %.2f\n', hold_decline_fraction);
fprintf('  Restoration implementation lag (years): %d\n', restoration_lag_years);


sens_labels = {'$k$ (attenuation)', '$D_{\\mathrm{max}}$ (\$M/yr)', '$M_0$ (\$/tourist)'};
mults       = [0.5, 1.0, 1.5];
sens_npv    = zeros(3,3);
for lv = 1:3
    for p = 1:3
        k_s  = k_atten * (1 + (p==1)*(mults(lv)-1));
        D_s  = D_max   * (1 + (p==2)*(mults(lv)-1));
        M0_s = M0      * (1 + (p==3)*(mults(lv)-1));
        Es   = @(C) E0 .* exp(-k_s .* C .* (A_reef_m2/1e6));
        dE_s = max(0, (Es(fcst_cover) - Es(C_ref)) / Es(C_ref));
        Dw_s = D_s .* (1 - exp(-lambda_dmg .* dE_s));
        dM_s = max(0, T0*M0_s*(1 - fcst_cover/A0));
        sens_npv(lv,p) = sum((Dw_s(:)+dM_s(:))./(1+r_disc).^t_disc)/1e6;
    end
end


DIV = repmat('=',1,76);
fprintf('\n%s\n  CORAL REEF ECONOMIC LOSS MODEL  |  %s\n%s\n\n',DIV,source_tag,DIV);
fprintf('LINEAR TOURISM MODEL\n%s\n',repmat('-',1,56));
fprintf('  T(A) = T0*(A/A0)   T0=%.2fM/yr  A0=%.2f%%\n',T0/1e6,A0);
fprintf('  M(T) = T*M0        M0=$%.0f/tourist\n',M0);
fprintf('  M_base = T0*M0   = $%.3fM/yr\n\n',M_base/1e6);
fprintf('FORECAST LOSSES\n%s\n',repmat('-',1,100));
fprintf('%-6s  %-8s  %-12s  %-14s  %-14s  %-14s  %-12s  %s\n',...
    'Year','Cover%','Tourists M','TouristLoss','WaveDmg $M','Tourism $M','Total $M','95% CI');
fprintf('%s\n',repmat('-',1,100));
for i = 1:n_fcst
    fprintf('%-6d  %-8.2f  %-12.3f  %-14.0f  %-14.2f  %-14.2f  %-12.2f  [%.1f-%.1f]\n',...
        fcst_years(i),fcst_cover(i),T_fcst(i)/1e6,dT_fcst(i),...
        D_wave_fcst(i)/1e6,dM_fcst(i)/1e6,L_fcst(i)/1e6,L_lo(i)/1e6,L_hi(i)/1e6);
end
fprintf('\nNPV  (%d yr, r=%.0f%%)\n',n_fcst,r_disc*100);
fprintf('  Wave  $%.2fM  |  Tourism  $%.2fM  |  TOTAL  $%.2fM  [CI $%.1f-$%.1fM]\n',...
    NPV_wave/1e6,NPV_tour/1e6,NPV_total/1e6,NPV_lo/1e6,NPV_hi/1e6);
fprintf('  Wave  $%.2fB  |  Tourism  $%.2fB  |  TOTAL  $%.2fB  [CI $%.2f-$%.2fB]\n',...
    NPV_wave/1e9,NPV_tour/1e9,NPV_total/1e9,NPV_lo/1e9,NPV_hi/1e9);
fprintf('\nMITIGATION\n');
for sc=1:3, fprintf('  %-20s  NPV $%.1fM ($%.2fB)  (save $%.1fM / $%.2fB)\n',...
    scenario_names{sc},npv_scenarios(sc),npv_scenarios(sc)/1e3,npv_saved(sc),npv_saved(sc)/1e3); end
fprintf('%s\n\n',DIV);


xlim_all  = [obs_years(1)-0.5, fcst_years(end)+0.5];
xbreak    = fcst_years(1) - 0.5;
plot_step_years = 5;
obs_mark_idx  = marker_idx_from_years(obs_years, plot_step_years);
fcst_mark_idx = marker_idx_from_years(fcst_years, plot_step_years);
all_mark_idx  = marker_idx_from_years(all_years, plot_step_years);
xticks_vec    = obs_years(1):plot_step_years:fcst_years(end);
xticks_vec    = xticks_vec(xticks_vec ~= 2098);
if xticks_vec(end) ~= fcst_years(end)
    xticks_vec = [xticks_vec fcst_years(end)];
end


% =========================================================================
% FIGURE E1 — Linear Tourism Model
% =========================================================================
figE1 = figure('Position',[60 60 1300 540],'Color','w');

ax1 = subplot(1,2,1); hold(ax1,'on');
fill([fcst_years, fliplr(fcst_years)], ...
     [T0.*(fcst_hi/A0)/1e6, fliplr(T0.*(fcst_lo/A0)/1e6)], ...
     C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax1, all_years, T_t/1e6, 'o-', ...
    'Color',C_TOUR,'MarkerFaceColor',C_TOUR, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',all_mark_idx, ...
    'DisplayName','$T(A)=T_0\,(A/A_0)$');
yline(ax1, T0/1e6, '--', 'Color',C_REF,'LineWidth',LW_REF, ...
    'HandleVisibility','off');
xline(ax1, xbreak, ':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax1, fcst_years(1)+0.5, T0/1e6*0.97, '$\leftarrow$ forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex');
styleAx(ax1,'Year','Tourists (millions\,yr$^{-1}$)','Visitor Count from Coral Cover');
legend(ax1,'Location','southwest');
xlim(ax1, xlim_all); ylim(ax1, [0, T0/1e6*1.18]);
xticks(ax1, xticks_vec);
text(ax1, fcst_years(end)-8, T0/1e6 * 1.035, ...
    sprintf('$T_0 = %.1f\\,\\mathrm{M\\,yr}^{-1}$', T0/1e6), ...
    'Interpreter','latex', 'FontSize',FS_TXT+1, ...
    'Color',C_REF, 'HorizontalAlignment','right', 'VerticalAlignment','bottom');

ax2 = subplot(1,2,2); hold(ax2,'on');
fill([fcst_years, fliplr(fcst_years)], ...
     [T0.*(fcst_hi/A0).*M0/1e6, fliplr(T0.*(fcst_lo/A0).*M0/1e6)], ...
     C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax2, all_years, M_t/1e6, 'o-', ...
    'Color',C_TOUR,'MarkerFaceColor',C_TOUR, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',all_mark_idx, ...
    'DisplayName','$M(A)=T_0 M_0\,(A/A_0)$');
yline(ax2, M_base/1e6, '--','Color',C_REF,'LineWidth',LW_REF, ...
    'HandleVisibility','off');
xline(ax2, xbreak, ':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax2, fcst_years(1)+0.5, M_base/1e6*0.97, '$\leftarrow$ forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex');
styleAx(ax2,'Year','Economic Output (\$M\,yr$^{-1}$)','Economic Output from Visitor Count');
legend(ax2,'Location','southwest');
xlim(ax2, xlim_all); ylim(ax2, [0, M_base/1e6*1.18]);
xticks(ax2, xticks_vec);
text(ax2, fcst_years(end)-8, M_base/1e6 * 1.035, ...
    sprintf('$M_{\\mathrm{base}} = %.1f\\,\\mathrm{\\$M\\,yr}^{-1}$', M_base/1e6), ...
    'Interpreter','latex', 'FontSize',FS_TXT+1, ...
    'Color',C_REF, 'HorizontalAlignment','right', 'VerticalAlignment','bottom');

makeSgtitle('Linear Tourism Model --- Visitor Count and Economic Output');
cleanFig(figE1);


% =========================================================================
% FIGURE E2 — Annual Losses Decomposed
% =========================================================================
figE2 = figure('Position',[60 60 1300 560],'Color','w');
y_ceil = max(L_hi/1e6)*1.35 + 2;

ax3 = subplot(1,2,1); hold(ax3,'on');
fill([fcst_years, fliplr(fcst_years)], [L_hi/1e6, fliplr(L_lo/1e6)], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
ah = area(ax3, obs_years, [D_wave_obs(:)/1e6, dM_obs(:)/1e6]);
ah(1).FaceColor = C_WAVE;  ah(1).FaceAlpha = 0.70;
ah(1).EdgeColor = 'w';     ah(1).LineWidth  = 0.5;
ah(2).FaceColor = C_TOUR;  ah(2).FaceAlpha = 0.70;
ah(2).EdgeColor = 'w';     ah(2).LineWidth  = 0.5;
plot(ax3, fcst_years, (D_wave_fcst+dM_fcst)/1e6, 'o--', ...
    'Color',C_TOTAL,'MarkerFaceColor','w', ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',fcst_mark_idx, ...
    'DisplayName','Total (forecast)');
xline(ax3, xbreak,':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax3, fcst_years(1)+0.5, y_ceil*0.83, '$\leftarrow$ forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex');
styleAx(ax3,'Year','Annual Economic Loss (\$M)','Stacked Annual Losses');
legend(ax3, [ah(1),ah(2)], {'Wave Damage','Tourism Loss'}, 'Location','northwest');
xlim(ax3, xlim_all); ylim(ax3, [0, y_ceil]);
xticks(ax3, xticks_vec);

ax4 = subplot(1,2,2); hold(ax4,'on');
fill([fcst_years, fliplr(fcst_years)], [L_hi/1e6, fliplr(L_lo/1e6)], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'DisplayName','95\% CI');
plot(ax4, all_years, D_wave/1e6, 's-', ...
    'Color',C_WAVE,'MarkerFaceColor',C_WAVE, ...
    'LineWidth',1.4,'MarkerSize',MS-1,'MarkerIndices',all_mark_idx, ...
    'DisplayName','Wave damage');
plot(ax4, all_years, dM_t/1e6, 'd-', ...
    'Color',C_TOUR,'MarkerFaceColor',C_TOUR, ...
    'LineWidth',1.4,'MarkerSize',MS-1,'MarkerIndices',all_mark_idx, ...
    'DisplayName','Tourism loss');
plot(ax4, all_years, L_total/1e6, 'o-', ...
    'Color',C_TOTAL,'MarkerFaceColor',C_TOTAL, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',all_mark_idx, ...
    'DisplayName','Total loss');
xline(ax4, xbreak,':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax4, fcst_years(1)+0.5, y_ceil*0.83, '$\leftarrow$ forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex');
styleAx(ax4,'Year','Annual Economic Loss (\$M)','Loss Components');
legend(ax4,'Location','northwest');
xlim(ax4, xlim_all); ylim(ax4, [0, y_ceil]);
xticks(ax4, xticks_vec);

makeSgtitle('Coral Reef Economic Losses: Wave Damage $+$ Tourism Revenue');
cleanFig(figE2);


% =========================================================================
% FIGURE E3 — Transfer Functions
% =========================================================================
figE3 = figure('Position',[60 60 1300 520],'Color','w');
cover_range = linspace(0, A0*1.1, 500);
dot_obs  = repmat(C_CORAL, n_obs,  1);
dot_fcst = repmat(C_FCST,  n_fcst, 1);

ax5 = subplot(1,3,1); hold(ax5,'on');
plot(ax5, cover_range, T0.*(cover_range/A0)/1e6, '-', ...
    'Color',C_TOUR,'LineWidth',LW_MAIN+0.2);
xline(ax5, A0,'--','Color',C_REF,'LineWidth',LW_REF, ...
    'HandleVisibility','off');
scatter(ax5, all_cover, T_t/1e6, 22, [dot_obs; dot_fcst], ...
    'filled','MarkerFaceAlpha',0.65);
text(ax5, A0*0.08, T0/1e6*0.78, ...
    sprintf('Slope $= T_0/A_0 = %.3f$ M per \\%%', T0/A0/1e6), ...
    'FontSize',FS_TXT,'Color',C_TOUR,'FontName',FONT,'Interpreter','latex');
styleAx(ax5,'Hard Coral Cover (\%)','Tourists (millions\,yr$^{-1}$)', ...
    '$T(A) = T_0 \cdot A/A_0$');
ylim(ax5,[0, T0/1e6*1.15]);
yl5 = ylim(ax5);
text(ax5, A0 + 0.45, yl5(2)*0.98, '$A_0$ baseline', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex', ...
    'BackgroundColor','w','Margin',1,'HorizontalAlignment','left');

ax6 = subplot(1,3,2); hold(ax6,'on');
dM_range = max(0, M_base*(1 - cover_range/A0)) / 1e6;
plot(ax6, cover_range, dM_range, '-', ...
    'Color',C_TOUR,'LineWidth',LW_MAIN+0.2);
xline(ax6, A0,'--','Color',C_REF,'LineWidth',LW_REF, ...
    'HandleVisibility','off');
scatter(ax6, all_cover, dM_t/1e6, 22, [dot_obs; dot_fcst], ...
    'filled','MarkerFaceAlpha',0.65);
text(ax6, 0.56, 0.90, ...
    sprintf('Slope $= -M_{\\mathrm{base}}/A_0$\n$= -%.2f$ M per \\%%', M_base/A0/1e6), ...
    'FontSize',FS_TXT,'Color',C_TOUR,'FontName',FONT,'Interpreter','latex', ...
    'Units','normalized','BackgroundColor','w','Margin',1);
styleAx(ax6,'Hard Coral Cover (\%)','Tourism Loss (\$M\,yr$^{-1}$)', ...
    '$\Delta M = M_{\mathrm{base}}(1 - A/A_0)$');
yl6 = ylim(ax6);
text(ax6, A0 + 0.45, yl6(2)*0.98, '$A_0$ baseline', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex', ...
    'BackgroundColor','w','Margin',1,'HorizontalAlignment','left');

ax7 = subplot(1,3,3); hold(ax7,'on');
dE_range = max(0, (E_shore(cover_range) - E_ref_val) / E_ref_val);
Dw_range = D_max .* (1 - exp(-lambda_dmg .* dE_range)) / 1e6;
plot(ax7, cover_range, Dw_range, '-', ...
    'Color',C_WAVE,'LineWidth',LW_MAIN+0.2);
xline(ax7, A0,'--','Color',C_REF,'LineWidth',LW_REF, ...
    'HandleVisibility','off');
scatter(ax7, all_cover, D_wave/1e6, 22, [dot_obs; dot_fcst], ...
    'filled','MarkerFaceAlpha',0.65);
text(ax7, 0.62, 0.90, ...
    sprintf('$D_{\\mathrm{max}} = %.0f$\,M', D_max/1e6), ...
    'FontSize',FS_TXT,'Color',C_WAVE,'FontName',FONT,'Interpreter','latex', ...
    'Units','normalized','BackgroundColor','w','Margin',1);
styleAx(ax7,'Hard Coral Cover (\%)','Wave Damage (\$M\,yr$^{-1}$)', ...
    '$D_{\mathrm{wave}} = D_{\mathrm{max}}[1 - e^{-\lambda\,\Delta E}]$');
yl7 = ylim(ax7);
text(ax7, A0 + 0.45, yl7(2)*0.98, '$A_0$ baseline', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex', ...
    'BackgroundColor','w','Margin',1,'HorizontalAlignment','left');

makeSgtitle('Transfer Functions: Coral Cover to Economic Loss');
cleanFig(figE3);


% =========================================================================
% FIGURE E4 — Cumulative Discounted NPV + Mitigation
% =========================================================================
figE4 = figure('Position',[60 60 1300 560],'Color','w');

ax8 = subplot(1,2,1); hold(ax8,'on');
cum_Dw = cumsum(D_wave_fcst(:)./(1+r_disc).^t_disc)/1e6;
cum_dM = cumsum(dM_fcst(:)./(1+r_disc).^t_disc)/1e6;
cum_lo = cumsum((D_wave_lo(:)+dM_lo(:))./(1+r_disc).^t_disc)/1e6;
cum_hi = cumsum((D_wave_hi(:)+dM_hi(:))./(1+r_disc).^t_disc)/1e6;
fill([fcst_years, fliplr(fcst_years)],[cum_hi', fliplr(cum_lo')], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
ac = area(ax8, fcst_years, [cum_Dw(:), cum_dM(:)]);
ac(1).FaceColor = C_WAVE; ac(1).FaceAlpha = 0.70;
ac(1).EdgeColor = 'w';    ac(1).LineWidth  = 0.5;
ac(2).FaceColor = C_TOUR; ac(2).FaceAlpha = 0.70;
ac(2).EdgeColor = 'w';    ac(2).LineWidth  = 0.5;
plot(ax8, fcst_years, cum_Dw+cum_dM, 'o-', ...
    'Color',C_TOTAL,'MarkerFaceColor',C_TOTAL, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',fcst_mark_idx, ...
    'DisplayName','Total NPV');
label_idx = marker_idx_from_years(fcst_years, 10);
for i = label_idx
    text(ax8, fcst_years(i), cum_Dw(i)+cum_dM(i)+max(cum_hi)*0.045, ...
        sprintf('\\$%.0fM', cum_Dw(i)+cum_dM(i)), ...
        'HorizontalAlignment','center','FontSize',FS_TXT,'Color',C_TOTAL, ...
        'FontName',FONT,'Interpreter','latex');
end
styleAx(ax8,'Year','Cumulative NPV (\$M)', ...
    sprintf('Cumulative Discounted Losses  ($r = %.0f\\%%$)', r_disc*100));
legend(ax8, [ac(1),ac(2)], {'Wave NPV','Tourism NPV'},'Location','northwest');
xlim(ax8, [fcst_years(1)-0.5, fcst_years(end)+0.5]);
ylim(ax8, [0, max(cum_hi)*1.25]);
xticks(ax8, xticks_vec);

ax9 = subplot(1,2,2); hold(ax9,'on');
bar_clrs = [C_TOTAL; C_WAVE; C_CORAL];
for sc = 1:3
    bh = bar(ax9, sc, npv_scenarios(sc), 0.55);
    bh.FaceColor = bar_clrs(sc,:);
    bh.FaceAlpha = 0.80;
    bh.EdgeColor = 'none';
end
for sc = 1:3
    text(ax9, sc, npv_scenarios(sc) + max(npv_scenarios)*0.03, ...
        sprintf('\\$%.1fM', npv_scenarios(sc)), ...
        'HorizontalAlignment','center','FontSize',FS_TXT,'Color','k', ...
        'FontName',FONT,'Interpreter','latex');
    if sc > 1 && npv_saved(sc) > 0.5
        text(ax9, sc, npv_scenarios(sc)/2, ...
            sprintf('Save \\$%.1fM', npv_saved(sc)), ...
            'HorizontalAlignment','center','FontSize',FS_TXT-1, ...
            'Color','w','FontWeight','bold','FontName',FONT,'Interpreter','latex');
    end
end
yline(ax9, npv_scenarios(1),'--','Color',C_REF,'LineWidth',LW_REF);
set(ax9,'XTick',1:3,'XTickLabel',scenario_names, ...
    'XTickLabelRotation',8,'TickDir','out', ...
    'Box','off','XColor','k','YColor','k', ...
    'FontName',FONT,'FontSize',FS_AX,'Color','w', ...
    'TickLabelInterpreter','latex');
ylabel(ax9,'NPV of Total Losses (\$M)', ...
    'Color','k','FontName',FONT,'FontSize',FS_LBL,'Interpreter','latex');
title(ax9,'Mitigation Scenarios --- NPV Comparison', ...
    'Color','k','FontName',FONT,'FontSize',12,'FontWeight','normal','Interpreter','latex');
ylim(ax9, [0, max(npv_scenarios)*1.22]);

makeSgtitle('NPV of Reef-Degradation Losses $+$ Mitigation Value');
cleanFig(figE4);


% =========================================================================
% FIGURE E5 — Dual-Axis Cover vs Total Loss
% =========================================================================
figE5 = figure('Position',[60 60 1100 500],'Color','w');
ax10 = axes(figE5); hold(ax10,'on');

left_ylim_hi  = max([A0, obs_cover, fcst_hi]) * 1.20;
right_ylim_hi = max(L_hi/1e6) * 1.20;

yyaxis(ax10,'left');
hCoverCI = fill([fcst_years, fliplr(fcst_years)], [fcst_hi, fliplr(fcst_lo)], ...
    C_CI,'EdgeColor','none','FaceAlpha',0.55,'HandleVisibility','off');
hCoverObs = plot(ax10, obs_years,  obs_cover,'o-', ...
    'Color',C_CORAL,'MarkerFaceColor',C_CORAL, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',obs_mark_idx, ...
    'DisplayName','Cover (observed)');
hCoverFcst = plot(ax10, fcst_years, fcst_cover,'o--', ...
    'Color',C_FCST,'MarkerFaceColor','w', ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',fcst_mark_idx, ...
    'DisplayName','Cover (forecast)');
yline(ax10, A0,':','Color',C_REF,'LineWidth',LW_REF,'HandleVisibility','off');
text(ax10, obs_years(1)+1.20, A0 + 0.06*left_ylim_hi, ...
    sprintf('$A_0 = %.1f\\%%$', A0), ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex');
ax10.YAxis(1).Color = 'k';
ylabel(ax10,'Hard Coral Cover (\%)','Color','k','FontName',FONT,'FontSize',FS_LBL,'Interpreter','latex');
ylim(ax10, [0, left_ylim_hi]);

yyaxis(ax10,'right');
hLossCI = fill([fcst_years, fliplr(fcst_years)], [L_hi/1e6, fliplr(L_lo/1e6)], ...
    C_CI,'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');
hLossObs = plot(ax10, obs_years,  L_obs/1e6,'s-', ...
    'Color',C_TOTAL,'MarkerFaceColor',C_TOTAL, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',obs_mark_idx, ...
    'DisplayName','Loss (observed)');
hLossFcst = plot(ax10, fcst_years, L_fcst/1e6,'s--', ...
    'Color',C_TOTAL,'MarkerFaceColor','w', ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'MarkerIndices',fcst_mark_idx, ...
    'DisplayName','Loss (forecast)');
ax10.YAxis(2).Color = 'k';
ylabel(ax10,'Total Annual Loss (\$M)','Color','k','FontName',FONT,'FontSize',FS_LBL,'Interpreter','latex');
ylim(ax10, [0, right_ylim_hi]);

xline(ax10, xbreak,':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax10, fcst_years(1)+0.8, right_ylim_hi*0.80, '$\leftarrow$ forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT,'Interpreter','latex');
xlabel(ax10,'Year','Color','k','FontName',FONT,'FontSize',FS_LBL,'Interpreter','latex');
lg5 = legend(ax10, [hCoverObs, hCoverFcst, hLossObs, hLossFcst], ...
    {'Cover (observed)','Cover (forecast)', ...
     'Total loss (observed)','Total loss (forecast)'}, ...
     'Location','northwest');
lg5.Box = 'off';
lg5.TextColor = 'k';
xlim(ax10, xlim_all);
xticks(ax10, xticks_vec);

ax10.Box      = 'off';
ax10.XGrid    = 'off';
ax10.YGrid    = 'off';
ax10.Color    = 'w';
ax10.XColor   = 'k';
ax10.TickDir  = 'out';
ax10.FontName = FONT;
ax10.FontSize = FS_AX;
ax10.LineWidth = 0.9;
ax10.TickLabelInterpreter = 'latex';

makeSgtitle('Coral Cover Decline and Associated Economic Losses');
cleanFig(figE5);


% =========================================================================
% CSV outputs
% =========================================================================
csv_obs_path  = fullfile(scriptDir, 'reef_economic_loss_observed.csv');
csv_fcst_path = fullfile(scriptDir, 'reef_economic_loss_forecast.csv');
csv_npv_path  = fullfile(scriptDir, 'reef_economic_npv.csv');
csv_mit_path  = fullfile(scriptDir, 'reef_economic_mitigation.csv');

writetable(table(obs_years',obs_cover',T_obs'/1e6,dT_obs', ...
    D_wave_obs'/1e6,dM_obs'/1e6,L_obs'/1e6, ...
    'VariableNames',{'Year','Cover_Pct','Tourists_M','Tourist_Loss', ...
        'WaveDamage_$M','TourismLoss_$M','TotalLoss_$M'}), ...
    csv_obs_path);

writetable(table(fcst_years',fcst_cover',T_fcst'/1e6,dT_fcst', ...
    D_wave_fcst'/1e6,dM_fcst'/1e6,L_fcst'/1e6,L_lo'/1e6,L_hi'/1e6, ...
    'VariableNames',{'Year','Cover_Pct','Tourists_M','Tourist_Loss', ...
        'WaveDamage_$M','TourismLoss_$M','TotalLoss_$M', ...
        'TotalLoss_Lower_$M','TotalLoss_Upper_$M'}), ...
    csv_fcst_path);

writetable(table({'Wave Protection';'Tourism Revenue';'TOTAL';'CI Low';'CI High'}, ...
    [NPV_wave;NPV_tour;NPV_total;NPV_lo;NPV_hi]/1e6, ...
    'VariableNames',{'Component','NPV_$M'}), csv_npv_path);

writetable(table(scenario_names',npv_scenarios',npv_saved', ...
    'VariableNames',{'Scenario','NPV_Loss_$M','NPV_Saved_$M'}), ...
    csv_mit_path);

fprintf('CSVs written:\n');
fprintf('  %s\n', csv_obs_path);
fprintf('  %s\n', csv_fcst_path);
fprintf('  %s\n', csv_npv_path);
fprintf('  %s\n\n', csv_mit_path);

figuresDir = fullfile(repoRoot, 'MATLAB', 'Figures');
if ~isfolder(figuresDir)
    mkdir(figuresDir);
end
fprintf('Figures output directory: %s\n', figuresDir);

pdf_paths = {
    fullfile(figuresDir, 'waveenergy_tourism_E1_linear_tourism_model.pdf'), ...
    fullfile(figuresDir, 'waveenergy_tourism_E2_annual_losses_decomposed_channels.pdf'), ...
    fullfile(figuresDir, 'waveenergy_tourism_E3_linear_tourism_transfer_function.pdf'), ...
    fullfile(figuresDir, 'waveenergy_tourism_E4_cumulative_discounted_npv_wave_tourism.pdf'), ...
    fullfile(figuresDir, 'waveenergy_tourism_E5_dual_axis_cover_vs_total_loss_full_series.pdf')
};
fig_handles = [figE1, figE2, figE3, figE4, figE5];
for i = 1:numel(fig_handles)
    if ~isgraphics(fig_handles(i), 'figure')
        warning('Figure handle is invalid or deleted for export: %s', pdf_paths{i});
        continue;
    end
    try
        exportgraphics(fig_handles(i), pdf_paths{i}, 'ContentType','vector');
    catch ME
        try
            set(fig_handles(i), 'PaperPositionMode','auto');
            print(fig_handles(i), pdf_paths{i}, '-dpdf', '-painters');
        catch ME2
            warning('Failed to export %s (%s | %s)', pdf_paths{i}, ME.message, ME2.message);
            continue;
        end
    end
    if isfile(pdf_paths{i})
        fprintf('  [OK] %s\n', pdf_paths{i});
    else
        warning('PDF was not written: %s', pdf_paths{i});
    end
end

fprintf('Vector PDFs written:\n');
for i = 1:numel(pdf_paths)
    fprintf('  %s\n', pdf_paths{i});
end
fprintf('\n');

fprintf('Figures saved:  E1 Linear Tourism Model\n');
fprintf('                E2 Annual Losses by Channel\n');
fprintf('                E3 Transfer Functions\n');
fprintf('                E4 Cumulative Discounted NPV (Wave + Tourism)\n');
fprintf('                E5 Dual-Axis Cover vs Total Loss (Full Series)\n');