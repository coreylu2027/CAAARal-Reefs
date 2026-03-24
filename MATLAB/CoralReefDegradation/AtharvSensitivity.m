% ═══════════════════════════════════════════════════════════════════════════
%  SENSITIVITY ANALYSIS — Coral Reef Economic Loss Model  (FIXED + IMPROVED)
%  Hardcoded from pipeline output. Standalone — no workspace needed.
%
%  SA1  Tornado chart:    NPV sensitivity to ±50% on k, D_max, M0
%  SA2  Discount rate:    r = 1.5%, 3%, 7% × 3 mitigation scenarios
%  SA3  Trend slope:      baseline ± 1 SD
% ═══════════════════════════════════════════════════════════════════════════

clc; close all;

% ── Palette ───────────────────────────────────────────────────────────────
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
    'defaultLegendBox','off','defaultLegendTextColor','k');

% ═══════════════════════════════════════════════════════════════════════════
%  A. DATA — hardcoded from CORIS/ARIMA pipeline output
% ═══════════════════════════════════════════════════════════════════════════
obs_years  = 2013:2023;
obs_cover  = [6.39,5.88,6.88,7.82,7.94,6.36,6.92,1.03,6.78,5.56,6.22];

fcst_years = 2024:2049;
fcst_cover = [5.32,6.18,5.58,4.74,2.97,3.95,3.48,4.27,5.62,4.38, ...
              4.53,3.44,1.79,2.71,2.26,3.03,4.36,3.14,3.28,2.20, ...
              0.58,1.48,1.04,1.79,3.09,1.89];

% Panel regression year coefficient SE (from bint of regress output)
trend_slope_base = -0.6397;
trend_slope_sd   =  0.2606;

% ═══════════════════════════════════════════════════════════════════════════
%  B. PARAMETERS  (must match main economic script)
% ═══════════════════════════════════════════════════════════════════════════
T0         = 2.1e6;          % Baseline tourists [persons/yr]
M0_param   = 105;            % Economic output per tourist [$/person]
M_base     = T0 * M0_param;  % Baseline tourism output [$]
E0         = 1.0e6;          % Offshore wave energy [kJ/m/yr]
k_base     = 0.004;          % Wave attenuation coefficient
A_reef_km  = 9.0;            % Reef area [km²]
A0         = obs_cover(1);   % Baseline cover = 6.39%
D_max_base = 850e6;          % Max coastal defence cost [$]
lam        = 2.5;            % Damage saturation rate
r_base     = 0.015;          % Social discount rate (OMB A-4)

n_fcst   = length(fcst_years);
t_disc_v = (1:n_fcst)';      % Column vector for discounting

% ── NPV helper functions ──────────────────────────────────────────────────
E_fn     = @(c, k)     E0 .* exp(-k .* c(:) .* A_reef_km);
Dw_fn    = @(c, k, Dm) Dm .* (1 - exp(-lam .* max(0, ...
               (E_fn(c,k) - E_fn(A0,k)) ./ E_fn(A0,k))));
dM_fn    = @(c, m0)    max(0, M_base - T0 .* m0 .* (c(:) ./ A0));
calc_npv = @(c, k, Dm, m0, r) ...
               sum((Dw_fn(c,k,Dm) + dM_fn(c,m0)) ./ (1+r).^t_disc_v) / 1e6;

% Also compute channel-split NPV for annotation
calc_npv_wave = @(c, k, Dm, r) ...
               sum(Dw_fn(c,k,Dm) ./ (1+r).^t_disc_v) / 1e6;
calc_npv_tour = @(c, m0, r) ...
               sum(dM_fn(c,m0) ./ (1+r).^t_disc_v) / 1e6;

% ═══════════════════════════════════════════════════════════════════════════
%  C. TORNADO: ±50% on k, D_max, M0
% ═══════════════════════════════════════════════════════════════════════════
NPV_base = calc_npv(fcst_cover, k_base, D_max_base, M0_param, r_base);
NPV_wave_base = calc_npv_wave(fcst_cover, k_base, D_max_base, r_base);
NPV_tour_base = calc_npv_tour(fcst_cover, M0_param, r_base);
fprintf('Baseline NPV: $%.1fM  (Wave: $%.1fM  Tourism: $%.1fM)\n', ...
    NPV_base, NPV_wave_base, NPV_tour_base);

% Compute lo/hi NPVs for each parameter
npv_k_lo = calc_npv(fcst_cover, k_base*0.5,      D_max_base,     M0_param,     r_base);
npv_k_hi = calc_npv(fcst_cover, k_base*1.5,      D_max_base,     M0_param,     r_base);
npv_D_lo = calc_npv(fcst_cover, k_base,           D_max_base*0.5, M0_param,     r_base);
npv_D_hi = calc_npv(fcst_cover, k_base,           D_max_base*1.5, M0_param,     r_base);
npv_M_lo = calc_npv(fcst_cover, k_base,           D_max_base,     M0_param*0.5, r_base);
npv_M_hi = calc_npv(fcst_cover, k_base,           D_max_base,     M0_param*1.5, r_base);

% Parameter names, lo/hi values, and ranges
pnames  = {'D_{max}  (coastal damage ceiling)', ...
           'k  (wave attenuation coeff.)', ...
           'M_0  (per-tourist output)'};
lo_vals = [npv_D_lo; npv_k_lo; npv_M_lo];
hi_vals = [npv_D_hi; npv_k_hi; npv_M_hi];

% Ensure lo < hi for each row (swap if needed — handles negative sensitivity)
for i = 1:3
    if lo_vals(i) > hi_vals(i)
        tmp = lo_vals(i); lo_vals(i) = hi_vals(i); hi_vals(i) = tmp;
    end
end
ranges = hi_vals - lo_vals;

% Sort descending by range
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

% ── Figure SA1 — Tornado ──────────────────────────────────────────────────
figSA1 = figure('Position',[60 60 1080 420],'Color','w');
ax1    = axes(figSA1); hold(ax1,'on');

bar_clrs = [C_TOTAL; C_WAVE; C_TOUR];

% FIX: y_pos must be INCREASING for set(ax,'YTick',...)
% Use y_pos = [1,2,3] with highest-impact bar at top (y=3)
% Sort order is already descending by range, so bar 1 = most impact = top
y_pos = [3, 2, 1];   % visual order: bar i drawn at y_pos(i)
bht   = 0.32;

for i = 1:3
    lo = lo_vals(i); hi = hi_vals(i); nb = NPV_base;
    % Left segment: lo → baseline (lighter)
    fill(ax1, [lo, nb, nb, lo], ...
         [y_pos(i)-bht, y_pos(i)-bht, y_pos(i)+bht, y_pos(i)+bht], ...
         bar_clrs(i,:), 'FaceAlpha',0.40, 'EdgeColor','none');
    % Right segment: baseline → hi (darker)
    fill(ax1, [nb, hi, hi, nb], ...
         [y_pos(i)-bht, y_pos(i)-bht, y_pos(i)+bht, y_pos(i)+bht], ...
         bar_clrs(i,:), 'FaceAlpha',0.88, 'EdgeColor','none');
    % Value labels
    pad_t = ranges(1) * 0.015;
    text(ax1, lo - pad_t, y_pos(i), sprintf('$%.0fM', lo), ...
        'HorizontalAlignment','right','FontSize',9,'Color','k','FontName',FONT);
    text(ax1, hi + pad_t, y_pos(i), sprintf('$%.0fM', hi), ...
        'HorizontalAlignment','left','FontSize',9,'Color','k','FontName',FONT);
    text(ax1, (lo+hi)/2, y_pos(i)+bht+0.10, sprintf('range  $%.0fM', ranges(i)), ...
        'HorizontalAlignment','center','FontSize',8.5,'Color','k','FontName',FONT);
end

% Baseline reference line
xline(ax1, NPV_base, '--', 'Color',C_REF, 'LineWidth',1.5);
text(ax1, NPV_base, 3+bht+0.38, sprintf('Baseline\n$%.0fM', NPV_base), ...
    'HorizontalAlignment','center','FontSize',9,'Color',C_REF,'FontName',FONT);

% Channel split annotation in bottom-right corner
text(ax1, max(hi_vals)*0.995, 0.55, ...
    sprintf('Wave NPV: $%.0fM (%.0f%%)\nTourism NPV: $%.0fM (%.0f%%)', ...
    NPV_wave_base, 100*NPV_wave_base/NPV_base, ...
    NPV_tour_base, 100*NPV_tour_base/NPV_base), ...
    'HorizontalAlignment','right','FontSize',8.5,'Color',C_REF,'FontName',FONT, ...
    'VerticalAlignment','bottom');

% FIX: YTick must be a strictly increasing vector
% We pass [1 2 3] and map labels in reverse so highest bar gets top label
set(ax1, 'YTick', [1 2 3], ...
         'YTickLabel', {pnames{3}, pnames{2}, pnames{1}}, ...
         'FontSize', 10.5, 'XColor','k','YColor','k', ...
         'Box','off','TickDir','out','Color','w');

xlabel(ax1,'NPV of Total Economic Losses ($M)','FontSize',11,'Color','k','FontName',FONT);
title(ax1,'Parameter Sensitivity — NPV Impact of \pm50% Variation', ...
    'FontSize',12,'FontWeight','normal','Color','k','FontName',FONT);

pad_x = ranges(1) * 0.25;
xlim(ax1, [min(lo_vals) - pad_x, max(hi_vals) + pad_x]);
ylim(ax1, [0.4, 3.8]);
sa_cleanFig(figSA1);

% ═══════════════════════════════════════════════════════════════════════════
%  D. DISCOUNT RATE SENSITIVITY
% ═══════════════════════════════════════════════════════════════════════════
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

% Savings relative to Status Quo per discount rate
fprintf('\nDISCOUNT RATE SENSITIVITY  (NPV $M)\n');
fprintf('%-12s  %-12s  %-12s  %-12s  %-12s  %-12s\n', ...
    'Rate','Status Quo','Hold Cover','Restoration','Save:Hold','Save:Rest');
for ri = 1:3
    fprintf('%-12s  %-12.1f  %-12.1f  %-12.1f  %-12.1f  %-12.1f\n', ...
        r_labels{ri}, npv_disc(ri,1), npv_disc(ri,2), npv_disc(ri,3), ...
        npv_disc(ri,1)-npv_disc(ri,2), npv_disc(ri,1)-npv_disc(ri,3));
end

% Cumulative NPV for status-quo at each rate
Dw_sq  = Dw_fn(fcst_cover, k_base, D_max_base);
dM_sq  = dM_fn(fcst_cover, M0_param);
cum_npv = zeros(3, n_fcst);
for ri = 1:3
    cum_npv(ri,:) = cumsum((Dw_sq + dM_sq) ./ (1+r_vals(ri)).^t_disc_v)' / 1e6;
end

r_clrs  = [C_TOTAL; C_CORAL; C_WAVE];
figSA2  = figure('Position',[60 60 1320 560],'Color','w');

% Left: cumulative NPV trajectories
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
% Annotate terminal values
for ri = 1:3
    text(ax2, fcst_years(end)+0.3, cum_npv(ri,end), ...
        sprintf('$%.0fM', cum_npv(ri,end)), ...
        'FontSize',8,'Color',r_clrs(ri,:),'FontName',FONT,'VerticalAlignment','middle');
end

% Right: grouped bar chart
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
% Savings callout arrows between Status Quo and Restoration for r=1.5%
sv_rest = npv_disc(1,1) - npv_disc(1,3);
annotation(figSA2,'textarrow', ...
    [0.745 0.745],[0.68 0.52], ...  % normalized figure coords (approx)
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

% ═══════════════════════════════════════════════════════════════════════════
%  E. TREND SLOPE SENSITIVITY  (benthic cover + economic loss)
% ═══════════════════════════════════════════════════════════════════════════
slope_opt  = trend_slope_base + trend_slope_sd;
slope_pess = trend_slope_base - trend_slope_sd;

make_cov   = @(sl) max(0, min(100, obs_cover(end) + sl*(1:n_fcst)));
cover_opt  = make_cov(slope_opt);
cover_pess = make_cov(slope_pess);
cover_base = fcst_cover;

npv_opt  = calc_npv(cover_opt,  k_base, D_max_base, M0_param, r_base);
npv_bsl  = calc_npv(cover_base, k_base, D_max_base, M0_param, r_base);
npv_pess = calc_npv(cover_pess, k_base, D_max_base, M0_param, r_base);

% Annual loss decomposition by scenario
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
    'Scenario','Slope (%/yr)','Cover 2049 (%)','Wave NPV ($M)','Tour NPV ($M)','NPV ($M)');
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

% ── Figure SA3 ─────────────────────────────────────────────────────────────
figSA3 = figure('Position',[60 60 1600 560],'Color','w');

% Left: benthic cover trajectories
ax4 = subplot(1,3,1); hold(ax4,'on');
fill([fcst_years, fliplr(fcst_years)], [cover_opt, fliplr(cover_pess)], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax4, obs_years,  obs_cover,  'o-', ...
    'Color',C_CORAL,'MarkerFaceColor',C_CORAL,'LineWidth',LW,'MarkerSize',MS, ...
    'DisplayName','Observed');
plot(ax4, fcst_years, cover_base, 'o--', ...
    'Color',C_CORAL,'MarkerFaceColor','w','LineWidth',LW,'MarkerSize',MS, ...
    'DisplayName',sprintf('Baseline (%.4f%%/yr)', trend_slope_base));
plot(ax4, fcst_years, cover_opt,  '-', 'Color',C_OPT,'LineWidth',LW, ...
    'DisplayName',sprintf('Optimistic (%.4f%%/yr)', slope_opt));
plot(ax4, fcst_years, cover_pess, '-', 'Color',C_PESS,'LineWidth',LW, ...
    'DisplayName',sprintf('Pessimistic (%.4f%%/yr)', slope_pess));
xline(ax4, fcst_years(1)-0.5,':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax4, fcst_years(1)+0.4, max(obs_cover)*1.05, '\leftarrow forecast', ...
    'FontSize',9,'Color',C_REF,'FontName',FONT);
xlabel(ax4,'Year','FontSize',11,'Color','k','FontName',FONT);
ylabel(ax4,'Hard Coral Cover (%)','FontSize',11,'Color','k','FontName',FONT);
title(ax4,'Cover Trajectories by Trend Scenario', ...
    'FontSize',12,'FontWeight','normal','Color','k','FontName',FONT);
legend(ax4,'Location','southwest','FontSize',9);
xlim(ax4,[obs_years(1)-1, fcst_years(end)+0.5]);
ylim(ax4,[0, max(obs_cover)*1.55]);

% Middle: wave damage channel
ax5 = subplot(1,3,2); hold(ax5,'on');
fill([fcst_years, fliplr(fcst_years)], [wave_opt, fliplr(wave_pess)], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax5, fcst_years, wave_base, 'o--', 'Color',C_WAVE,'MarkerFaceColor','w', ...
    'LineWidth',LW,'MarkerSize',MS, ...
    'DisplayName',sprintf('Baseline  NPV $%.0fM', calc_npv_wave(cover_base,k_base,D_max_base,r_base)));
plot(ax5, fcst_years, wave_opt,  '-', 'Color',C_OPT,'LineWidth',LW, ...
    'DisplayName',sprintf('Optimistic  NPV $%.0fM', calc_npv_wave(cover_opt,k_base,D_max_base,r_base)));
plot(ax5, fcst_years, wave_pess, '-', 'Color',C_PESS,'LineWidth',LW, ...
    'DisplayName',sprintf('Pessimistic  NPV $%.0fM', calc_npv_wave(cover_pess,k_base,D_max_base,r_base)));
xlabel(ax5,'Year','FontSize',11,'Color','k','FontName',FONT);
ylabel(ax5,'Wave Damage Loss ($M/yr)','FontSize',11,'Color','k','FontName',FONT);
title(ax5,'Wave Damage Channel by Trend Scenario', ...
    'FontSize',12,'FontWeight','normal','Color','k','FontName',FONT);
legend(ax5,'Location','northwest','FontSize',9);
xlim(ax5,[fcst_years(1)-0.5, fcst_years(end)+0.5]);
ylim(ax5,[0, max([wave_opt(:); wave_pess(:)])*1.30]);

% Right: tourism channel + total
ax6 = subplot(1,3,3); hold(ax6,'on');
fill([fcst_years, fliplr(fcst_years)], [loss_opt, fliplr(loss_pess)], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax6, fcst_years, loss_base, 'o--', 'Color',C_TOTAL,'MarkerFaceColor','w', ...
    'LineWidth',LW,'MarkerSize',MS, ...
    'DisplayName',sprintf('Baseline  NPV $%.0fM', npv_bsl));
plot(ax6, fcst_years, loss_opt,  '-', 'Color',C_OPT,'LineWidth',LW, ...
    'DisplayName',sprintf('Optimistic  NPV $%.0fM', npv_opt));
plot(ax6, fcst_years, loss_pess, '-', 'Color',C_PESS,'LineWidth',LW, ...
    'DisplayName',sprintf('Pessimistic  NPV $%.0fM', npv_pess));
% Overlay tourism-only dashed
plot(ax6, fcst_years, tour_base, '--', 'Color',C_TOUR,'LineWidth',1.2, ...
    'DisplayName','Tourism only (baseline)');
xlabel(ax6,'Year','FontSize',11,'Color','k','FontName',FONT);
ylabel(ax6,'Total Annual Loss ($M)','FontSize',11,'Color','k','FontName',FONT);
title(ax6,'Total Loss (Wave + Tourism) by Trend Scenario', ...
    'FontSize',12,'FontWeight','normal','Color','k','FontName',FONT);
legend(ax6,'Location','northwest','FontSize',9);
xlim(ax6,[fcst_years(1)-0.5, fcst_years(end)+0.5]);
ylim(ax6,[0, max([loss_opt(:); loss_pess(:)])*1.28]);

sg3 = sgtitle('Sensitivity to Hard Coral Cover Trend Slope (\pm1 SD)  |  Wave + Tourism Channels');
sg3.FontName = FONT; sg3.FontSize = 13; sg3.FontWeight = 'normal'; sg3.Color = 'k';
sa_cleanFig(figSA3);

% ═══════════════════════════════════════════════════════════════════════════
%  F. BONUS — Channel Split Sensitivity Comparison
%     Shows how each ±50% parameter swing affects Wave vs Tourism NPV
%     independently, giving intuition for which channel drives each parameter.
% ═══════════════════════════════════════════════════════════════════════════
figSA4 = figure('Position',[60 60 1100 480],'Color','w');
ax7 = axes(figSA4); hold(ax7,'on');

param_labels = {'D_{max} \times0.5', 'D_{max} \times1.5', ...
                'k \times0.5',       'k \times1.5', ...
                'M_0 \times0.5',     'M_0 \times1.5'};

wave_vals = [ ...
    calc_npv_wave(fcst_cover, k_base, D_max_base*0.5, r_base), ...
    calc_npv_wave(fcst_cover, k_base, D_max_base*1.5, r_base), ...
    calc_npv_wave(fcst_cover, k_base*0.5, D_max_base, r_base), ...
    calc_npv_wave(fcst_cover, k_base*1.5, D_max_base, r_base), ...
    calc_npv_wave(fcst_cover, k_base, D_max_base, r_base), ...   % M0 doesn't affect wave
    calc_npv_wave(fcst_cover, k_base, D_max_base, r_base)];

tour_vals = [ ...
    calc_npv_tour(fcst_cover, M0_param, r_base), ...  % D_max doesn't affect tourism
    calc_npv_tour(fcst_cover, M0_param, r_base), ...
    calc_npv_tour(fcst_cover, M0_param, r_base), ...  % k doesn't affect tourism
    calc_npv_tour(fcst_cover, M0_param, r_base), ...
    calc_npv_tour(fcst_cover, M0_param*0.5, r_base), ...
    calc_npv_tour(fcst_cover, M0_param*1.5, r_base)];

x_pos = 1:6;
bw = 0.32;

bw1 = bar(ax7, x_pos - bw/2, wave_vals, bw);
bw1.FaceColor = C_WAVE;  bw1.FaceAlpha = 0.85; bw1.EdgeColor = 'none';
bw1.DisplayName = 'Wave Damage NPV';

bw2 = bar(ax7, x_pos + bw/2, tour_vals, bw);
bw2.FaceColor = C_TOUR;  bw2.FaceAlpha = 0.85; bw2.EdgeColor = 'none';
bw2.DisplayName = 'Tourism Loss NPV';

% Baseline reference lines
yline(ax7, NPV_wave_base, '--', 'Color',C_WAVE, 'LineWidth',1.2, ...
    'Label',sprintf('Wave base $%.0fM', NPV_wave_base), ...
    'LabelHorizontalAlignment','left','HandleVisibility','off');
yline(ax7, NPV_tour_base, '--', 'Color',C_TOUR, 'LineWidth',1.2, ...
    'Label',sprintf('Tourism base $%.0fM', NPV_tour_base), ...
    'LabelHorizontalAlignment','left','HandleVisibility','off');

% Value labels on bars
for i = 1:6
    if abs(wave_vals(i) - NPV_wave_base) > 1
        text(ax7, x_pos(i)-bw/2, wave_vals(i)+NPV_wave_base*0.02, ...
            sprintf('$%.0f', wave_vals(i)), ...
            'HorizontalAlignment','center','FontSize',8,'Color','k','FontName',FONT);
    end
    text(ax7, x_pos(i)+bw/2, tour_vals(i)+NPV_wave_base*0.02, ...
        sprintf('$%.0f', tour_vals(i)), ...
        'HorizontalAlignment','center','FontSize',8,'Color','k','FontName',FONT);
end

set(ax7,'XTick',x_pos,'XTickLabel',param_labels,'FontSize',10.5, ...
    'Box','off','TickDir','out','XColor','k','YColor','k','Color','w', ...
    'FontName',FONT);
ylabel(ax7,'NPV of Losses ($M)','FontSize',11,'Color','k','FontName',FONT);
title(ax7,'Channel-Split Sensitivity: Wave vs Tourism NPV under ±50% Parameter Variation', ...
    'FontSize',12,'FontWeight','normal','Color','k','FontName',FONT);
legend(ax7,'Location','northeast');
ylim(ax7,[0, max([wave_vals, tour_vals])*1.22]);
sa_cleanFig(figSA4);

fprintf('\n%s\n  DONE\n', repmat('=',1,60));
fprintf('  SA1 Tornado (total NPV + channel breakdown)\n');
fprintf('  SA2 Discount Rate (cumulative NPV + grouped bar)\n');
fprintf('  SA3 Trend Slope (cover + wave channel + total loss)\n');
fprintf('  SA4 Channel Split (wave vs tourism under each parameter swing)\n');
fprintf('%s\n', repmat('=',1,60));

% ═══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  (must stay at end of file)
% ═══════════════════════════════════════════════════════════════════════════
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