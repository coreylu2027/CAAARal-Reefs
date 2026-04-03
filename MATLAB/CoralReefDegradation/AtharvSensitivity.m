

clc; close all;


C_CORAL = [0.13 0.39 0.68];
C_WAVE  = [0.20 0.60 0.86];
C_TOUR  = [0.90 0.45 0.00];
C_TOTAL = [0.80 0.17 0.17];
C_REF   = [0.40 0.40 0.40];
C_CI    = [0.85 0.85 0.85];
C_OPT   = [0.13 0.60 0.32];
C_PESS  = [0.80 0.17 0.17];
FONT    = 'Helvetica';
LW      = 1.8;
MS      = 7;

set(groot,'defaultAxesColor','w','defaultFigureColor','w', ...
    'defaultAxesBox','off','defaultAxesXGrid','off','defaultAxesYGrid','off', ...
    'defaultAxesFontName',FONT,'defaultAxesFontSize',11, ...
    'defaultAxesXColor','k','defaultAxesYColor','k', ...
    'defaultAxesLineWidth',0.9,'defaultAxesTickDir','out', ...
    'defaultTextFontName',FONT,'defaultTextColor','k', ...
    'defaultTextInterpreter','latex', ...
    'defaultLegendInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendBox','off','defaultLegendTextColor','k');


scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));

% Keep sensitivity inputs synchronized with the benthic outputs generated in MATLAB/CoralReefDegradation.
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
else
    obs_years  = 2013:2023;
    obs_cover  = [6.39,5.88,6.88,7.82,7.94,6.36,6.92,1.03,6.78,5.56,6.22];

    fcst_years = 2024:2049;
    fcst_cover = [5.32,6.18,5.58,4.74,2.97,3.95,3.48,4.27,5.62,4.38, ...
                  4.53,3.44,1.79,2.71,2.26,3.03,4.36,3.14,3.28,2.20, ...
                  0.58,1.48,1.04,1.79,3.09,1.89];
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
    prev_c = fcst_cover(end);
    for ii = 1:n_add
        cand = prev_c + tail_slope;
        if prev_c <= 0 || cand <= 0
            cand = 0;
        end
        cov_add(ii) = max(0, min(100, cand));
        prev_c = cov_add(ii);
    end

    fcst_years = [fcst_years, yr_add];
    fcst_cover = [fcst_cover, cov_add];
elseif fcst_years(end) > target_forecast_year
    keep_idx = fcst_years <= target_forecast_year;
    fcst_years = fcst_years(keep_idx);
    fcst_cover = fcst_cover(keep_idx);
end


trend_slope_base = -0.6397;
trend_slope_sd   =  0.2606;


T0         = 2.1e6;
M0_param   = 1800;
M_base     = T0 * M0_param;
E0         = 1.0e6;
k_base     = 0.004;
A_reef_km  = 9.0;
A0         = obs_cover(1);
D_max_base = 5.0e9;
lam        = 2.5;
r_base     = 0.015;

n_fcst   = length(fcst_years);
t_disc_v = (1:n_fcst)';      % Column vector for discounting


E_fn     = @(c, k)     E0 .* exp(-k .* c(:) .* A_reef_km);
Dw_fn    = @(c, k, Dm) Dm .* (1 - exp(-lam .* max(0, ...
               (E_fn(c,k) - E_fn(A0,k)) ./ E_fn(A0,k))));
dM_fn    = @(c, m0)    max(0, M_base - T0 .* m0 .* (c(:) ./ A0));
calc_npv = @(c, k, Dm, m0, r) ...
               sum((Dw_fn(c,k,Dm) + dM_fn(c,m0)) ./ (1+r).^t_disc_v) / 1e6;


calc_npv_wave = @(c, k, Dm, r) ...
               sum(Dw_fn(c,k,Dm) ./ (1+r).^t_disc_v) / 1e6;
calc_npv_tour = @(c, m0, r) ...
               sum(dM_fn(c,m0) ./ (1+r).^t_disc_v) / 1e6;


NPV_base = calc_npv(fcst_cover, k_base, D_max_base, M0_param, r_base);
NPV_wave_base = calc_npv_wave(fcst_cover, k_base, D_max_base, r_base);
NPV_tour_base = calc_npv_tour(fcst_cover, M0_param, r_base);
fprintf('Baseline NPV: $%.1fM  (Wave: $%.1fM  Tourism: $%.1fM)\n', ...
    NPV_base, NPV_wave_base, NPV_tour_base);


npv_k_lo = calc_npv(fcst_cover, k_base*0.5,      D_max_base,     M0_param,     r_base);
npv_k_hi = calc_npv(fcst_cover, k_base*1.5,      D_max_base,     M0_param,     r_base);
npv_D_lo = calc_npv(fcst_cover, k_base,           D_max_base*0.5, M0_param,     r_base);
npv_D_hi = calc_npv(fcst_cover, k_base,           D_max_base*1.5, M0_param,     r_base);
npv_M_lo = calc_npv(fcst_cover, k_base,           D_max_base,     M0_param*0.5, r_base);
npv_M_hi = calc_npv(fcst_cover, k_base,           D_max_base,     M0_param*1.5, r_base);


pnames  = {'D_{max}  (coastal damage ceiling)', ...
           'k  (wave attenuation coeff.)', ...
           'M_0  (per-tourist output)'};
lo_vals = [npv_D_lo; npv_k_lo; npv_M_lo];
hi_vals = [npv_D_hi; npv_k_hi; npv_M_hi];


for i = 1:3
    if lo_vals(i) > hi_vals(i)
        tmp = lo_vals(i); lo_vals(i) = hi_vals(i); hi_vals(i) = tmp;
    end
end
ranges = hi_vals - lo_vals;


[~, ord] = sort(ranges, 'descend');
pnames   = pnames(ord);
lo_vals  = lo_vals(ord);
hi_vals  = hi_vals(ord);
ranges   = ranges(ord);

fprintf('\nTORNADO  (NPV in $M)\n');
fprintf('%-42s  %9s  %9s  %9s\n','Parameter','Lo','Hi','Range');
for i = 1:3
    fprintf('%-42s  %9.1f  %9.1f  %9.1f\n', pnames{i}, lo_vals(i), hi_vals(i), ranges(i));
end


figSA1 = figure('Position',[60 60 1080 420],'Color','w');
ax1    = axes(figSA1); hold(ax1,'on');

% Draw each tornado bar as two blocks around the baseline NPV.
bar_clrs = [C_TOTAL; C_WAVE; C_TOUR];


y_pos = [3, 2, 1];
bht   = 0.32;

for i = 1:3
    lo = lo_vals(i); hi = hi_vals(i); nb = NPV_base;

    fill(ax1, [lo, nb, nb, lo], ...
         [y_pos(i)-bht, y_pos(i)-bht, y_pos(i)+bht, y_pos(i)+bht], ...
         bar_clrs(i,:), 'FaceAlpha',0.40, 'EdgeColor','none');

    fill(ax1, [nb, hi, hi, nb], ...
         [y_pos(i)-bht, y_pos(i)-bht, y_pos(i)+bht, y_pos(i)+bht], ...
         bar_clrs(i,:), 'FaceAlpha',0.88, 'EdgeColor','none');

    pad_t = ranges(1) * 0.015;
    text(ax1, lo - pad_t, y_pos(i), sprintf('$%.0fM', lo), ...
        'HorizontalAlignment','right','FontSize',9,'Color','k','FontName',FONT);
    text(ax1, hi + pad_t, y_pos(i)+0.11, sprintf('$%.0fM', hi), ...
        'HorizontalAlignment','left','FontSize',9,'Color','k','FontName',FONT);
    text(ax1, hi + 2.2*pad_t, y_pos(i)-0.11, sprintf('range  $%.0fM', ranges(i)), ...
        'HorizontalAlignment','left','FontSize',8.5,'Color','k','FontName',FONT);
end


xline(ax1, NPV_base, '--', 'Color',C_REF, 'LineWidth',1.5);
baseline_dx = ranges(1) * 0.08;
text(ax1, NPV_base + 2.2*baseline_dx, 3+bht+0.20, sprintf('Baseline\n$%.0fM', NPV_base), ...
    'HorizontalAlignment','left','FontSize',9,'Color',C_REF,'FontName',FONT);


text(ax1, max(hi_vals)*0.995, 0.55, ...
    sprintf('Wave NPV: $%.0fM (%.0f%%)\nTourism NPV: $%.0fM (%.0f%%)', ...
    NPV_wave_base, 100*NPV_wave_base/NPV_base, ...
    NPV_tour_base, 100*NPV_tour_base/NPV_base), ...
    'HorizontalAlignment','right','FontSize',8.5,'Color',C_REF,'FontName',FONT, ...
    'VerticalAlignment','bottom');


set(ax1, 'YTick', [1 2 3], ...
         'YTickLabel', {pnames{3}, pnames{2}, pnames{1}}, ...
         'FontSize', 10.5, 'XColor','k','YColor','k', ...
         'Box','off','TickDir','out','Color','w');

xlabel(ax1,'NPV of Total Economic Losses ($M)','FontSize',11,'Color','k','FontName',FONT);
title(ax1,'Tornado Sensitivity: NPV Impact of $\pm$50\% Parameter Variation', ...
    'FontSize',12,'FontWeight','normal','Color','k','FontName',FONT);

pad_x = ranges(1) * 0.25;
xlim(ax1, [min(lo_vals) - pad_x, max(hi_vals) + pad_x]);
ylim(ax1, [0.4, 3.8]);
sa_cleanFig(figSA1);


r_vals     = [0.015, 0.03, 0.07];
r_labels   = {'r = 1.5%', 'r = 3%', 'r = 7%'};
scen_names = {'Status Quo', 'Hold Cover', 'Restoration'};

cover_sq   = fcst_cover;
cover_hold = repmat(obs_cover(end), 1, n_fcst);
cover_rest = linspace(obs_cover(end), A0, n_fcst);
scenarios  = {cover_sq, cover_hold, cover_rest};

npv_disc = zeros(3,3);
for ri = 1:3
    for sc = 1:3
        npv_disc(ri,sc) = calc_npv(scenarios{sc}, k_base, D_max_base, M0_param, r_vals(ri));
    end
end


fprintf('\nDISCOUNT RATE SENSITIVITY  (NPV $M)\n');
fprintf('%-12s  %-12s  %-12s  %-12s  %-12s  %-12s\n', ...
    'Rate','Status Quo','Hold Cover','Restoration','Save:Hold','Save:Rest');
for ri = 1:3
    fprintf('%-12s  %-12.1f  %-12.1f  %-12.1f  %-12.1f  %-12.1f\n', ...
        r_labels{ri}, npv_disc(ri,1), npv_disc(ri,2), npv_disc(ri,3), ...
        npv_disc(ri,1)-npv_disc(ri,2), npv_disc(ri,1)-npv_disc(ri,3));
end


Dw_sq  = Dw_fn(fcst_cover, k_base, D_max_base);
dM_sq  = dM_fn(fcst_cover, M0_param);
cum_npv = zeros(3, n_fcst);
for ri = 1:3
    cum_npv(ri,:) = cumsum((Dw_sq + dM_sq) ./ (1+r_vals(ri)).^t_disc_v)' / 1e6;
end

r_clrs  = [C_TOTAL; C_CORAL; C_WAVE];
figSA2  = figure('Position',[60 60 1320 560],'Color','w');


ax2 = subplot(1,2,1); hold(ax2,'on');
for ri = 1:3
    plot(ax2, fcst_years, cum_npv(ri,:), 'o-', ...
        'Color',r_clrs(ri,:),'MarkerFaceColor',r_clrs(ri,:), ...
        'LineWidth',LW,'MarkerSize',MS-1,'DisplayName',r_labels{ri});
end
xlabel(ax2,'Year','FontSize',11,'Color','k','FontName',FONT);
ylabel(ax2,'Cumulative Discounted Loss ($M)','FontSize',11,'Color','k','FontName',FONT);
title(ax2,'Cumulative NPV by Discount Rate (Status Quo)', ...
    'FontSize',12,'FontWeight','normal','Color','k','FontName',FONT);
legend(ax2,'Location','northwest');
xlim(ax2,[fcst_years(1)-0.5, fcst_years(end)+0.5]);
ylim(ax2,[0, max(cum_npv(:))*1.22]);

for ri = 1:3
    text(ax2, fcst_years(end)+1.2, cum_npv(ri,end), ...
        sprintf('$%.0fM', cum_npv(ri,end)), ...
        'FontSize',8,'Color',r_clrs(ri,:),'FontName',FONT,'VerticalAlignment','middle');
end
xlim(ax2,[fcst_years(1)-0.5, fcst_years(end)+2.5]);


ax3   = subplot(1,2,2); hold(ax3,'on');
bar_w = 0.22;
x_grp = 1:3;
offs  = [-bar_w, 0, bar_w];
for ri = 1:3
    bh2 = bar(ax3, x_grp+offs(ri), npv_disc(ri,:), bar_w*0.92);
    bh2.FaceColor   = r_clrs(ri,:);
    bh2.FaceAlpha   = 0.80;
    bh2.EdgeColor   = 'none';
    bh2.DisplayName = r_labels{ri};
end
for ri = 1:3
    for sc = 1:3
        text(ax3, x_grp(sc)+offs(ri), npv_disc(ri,sc)+max(npv_disc(:))*0.025, ...
            sprintf('%.0f', npv_disc(ri,sc)), ...
            'HorizontalAlignment','center','FontSize',8,'Color','k','FontName',FONT);
    end
end

sv_rest = npv_disc(1,1) - npv_disc(1,3);
annotation(figSA2,'textarrow', ...
    [0.745 0.745],[0.68 0.52], ...
    'String', sprintf('Save\n$%.0fM', sv_rest), ...
    'FontSize',8,'FontName',FONT,'Color',C_OPT, ...
    'HeadStyle','vback2','HeadWidth',6,'HeadLength',6,'LineWidth',1);

set(ax3,'XTick',x_grp,'XTickLabel',scen_names,'XTickLabelRotation',8, ...
    'Box','off','TickDir','out','XColor','k','YColor','k','FontName',FONT,'FontSize',11);
ylabel(ax3,'NPV ($M)','FontSize',11,'Color','k','FontName',FONT);
title(ax3,'NPV by Discount Rate and Mitigation Scenario', ...
    'FontSize',12,'FontWeight','normal','Color','k','FontName',FONT);
legend(ax3,'Location','northeast');
ylim(ax3,[0, max(npv_disc(:))*1.22]);

sg2 = sgtitle('Sensitivity to Social Discount Rate (OMB A-4: 1.5%, 3%, 7%)');
sg2.FontName = FONT; sg2.FontSize = 13; sg2.FontWeight = 'normal'; sg2.Color = 'k';
sa_cleanFig(figSA2);


slope_opt  = trend_slope_base + trend_slope_sd;
slope_pess = trend_slope_base - trend_slope_sd;

make_cov   = @(sl) max(0, min(100, obs_cover(end) + sl*(1:n_fcst)));
cover_opt  = make_cov(slope_opt);
cover_pess = make_cov(slope_pess);
cover_base = fcst_cover;

npv_opt  = calc_npv(cover_opt,  k_base, D_max_base, M0_param, r_base);
npv_bsl  = calc_npv(cover_base, k_base, D_max_base, M0_param, r_base);
npv_pess = calc_npv(cover_pess, k_base, D_max_base, M0_param, r_base);


ann_wave = @(cov) Dw_fn(cov, k_base, D_max_base)' / 1e6;
ann_tour = @(cov) dM_fn(cov, M0_param)' / 1e6;
ann_total = @(cov) ann_wave(cov) + ann_tour(cov);

wave_opt  = ann_wave(cover_opt);   tour_opt  = ann_tour(cover_opt);
wave_base = ann_wave(cover_base);  tour_base = ann_tour(cover_base);
wave_pess = ann_wave(cover_pess);  tour_pess = ann_tour(cover_pess);

loss_opt  = ann_total(cover_opt);
loss_base = ann_total(cover_base);
loss_pess = ann_total(cover_pess);

fprintf('\nTREND SLOPE SENSITIVITY\n');
fprintf('%-28s  %-12s  %-16s  %-14s  %-14s  %-10s\n', ...
    'Scenario','Slope (%/yr)',sprintf('Cover %d (%%)',target_forecast_year),'Wave NPV ($M)','Tour NPV ($M)','NPV ($M)');
fprintf('%-28s  %-12.4f  %-16.2f  %-14.1f  %-14.1f  %-10.1f\n', ...
    'Optimistic (+1 SD)',  slope_opt,  cover_opt(end), ...
    calc_npv_wave(cover_opt,k_base,D_max_base,r_base), ...
    calc_npv_tour(cover_opt,M0_param,r_base), npv_opt);
fprintf('%-28s  %-12.4f  %-16.2f  %-14.1f  %-14.1f  %-10.1f\n', ...
    'Baseline (ARIMA)', trend_slope_base, cover_base(end), ...
    NPV_wave_base, NPV_tour_base, npv_bsl);
fprintf('%-28s  %-12.4f  %-16.2f  %-14.1f  %-14.1f  %-10.1f\n', ...
    'Pessimistic (-1 SD)', slope_pess, cover_pess(end), ...
    calc_npv_wave(cover_pess,k_base,D_max_base,r_base), ...
    calc_npv_tour(cover_pess,M0_param,r_base), npv_pess);

fprintf('\n%s\n  DONE\n', repmat('=',1,60));
fprintf('  SA1 Tornado (total NPV + channel breakdown)\n');
fprintf('  SA2 Discount Rate (cumulative NPV + grouped bar)\n');
fprintf('%s\n', repmat('=',1,60));


function sa_cleanFig(fig)
    set(fig,'Color','w');
    for ax = findall(fig,'Type','axes')'
        set(ax,'Color','w','XColor','k','YColor','k','ZColor','k', ...
            'Box','off','XGrid','off','YGrid','off','TickDir','out', ...
            'LineWidth',0.9,'FontName','Helvetica','FontSize',11);
        ax.Title.Color    = 'k'; ax.Title.FontWeight = 'normal';
        ax.Title.FontSize = 12;  ax.Title.FontName   = 'Helvetica';
        ax.XLabel.Color   = 'k'; ax.XLabel.FontName  = 'Helvetica';
        ax.YLabel.Color   = 'k'; ax.YLabel.FontName  = 'Helvetica';
    end
    set(findall(fig,'Type','text'),'Color','k','FontName','Helvetica');
    for lg = findall(fig,'Type','legend')'
        set(lg,'Box','off','TextColor','k','Color','none','EdgeColor','none');
    end
end