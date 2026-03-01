% ═══════════════════════════════════════════════════════════════════════════
%  CORAL REEF ECONOMIC LOSS MODEL
%  Wave-Energy Coastal Protection Valuation  +  Linear Tourism Revenue Loss
%
%  Run AFTER the main benthic pipeline so that workspace variables
%  (years_all, annual_mean_hc, future_years_vec, yF, yF_lower, yF_upper)
%  are available.  Falls back to a synthetic FL Keys series if not found.
%
% ─────────────────────────────────────────────────────────────────────────
%  MATHEMATICAL FRAMEWORK
% ─────────────────────────────────────────────────────────────────────────
%
%  1. WAVE ENERGY ATTENUATION  (Ferrario et al. 2014; Storlazzi et al. 2019)
%
%       E_shore(t) = E0 · exp( −k · C(t) · A_reef )
%
%     E0      = incoming offshore wave energy  [kJ m⁻¹ yr⁻¹]
%     k       = bulk attenuation coefficient   [per % cover per km²]
%     C(t)    = hard coral cover  [%]
%     A_reef  = reef plan-form area  [km²]
%
%     Fractional excess shore energy above healthy-reef baseline:
%       ΔE_frac(t) = max(0, [E_shore(t) − E_shore(C_ref)] / E_shore(C_ref))
%
%     Annual coastal-damage cost  (saturating damage function):
%       D_wave(t) = D_max · [1 − exp(−λ · ΔE_frac(t))]
%
%  2. LINEAR TOURISM MODEL
%
%       T(A) = T0 · (A / A0)            tourists proportional to cover
%       M(T) = T  · M0                  output proportional to tourists
%
%     Combined:  M(A) = T0 · M0 · (A / A0)
%     Baseline:  M_base = T0 · M0       (at A = A0, zero loss)
%
%     Annual tourism economic loss:
%       ΔM(t) = max(0,  M_base · [1 − C(t)/A0])
%
%  3. TOTAL ANNUAL ECONOMIC LOSS:
%       L_total(t) = D_wave(t) + ΔM(t)
%
%  4. NET PRESENT VALUE  (US OMB A-4, r = 3 %):
%       NPV = Σ_{τ=1}^{T}  L_total(t_τ) / (1 + r)^τ
%
%  5. MITIGATION SCENARIOS:
%       (a) Status Quo   — cover follows ARIMA forecast
%       (b) Hold Cover   — no further decline from current level
%       (c) Restoration  — linear recovery to A0 over forecast window
%
%  6. SENSITIVITY: ±50 % on k, D_max, M0
%
% ═══════════════════════════════════════════════════════════════════════════

clc;
clear;
close all;

% ── Palette ───────────────────────────────────────────────────────────────
C_CORAL  = [0.13 0.39 0.68];
C_WAVE   = [0.20 0.60 0.86];
C_TOUR   = [0.90 0.45 0.00];
C_TOTAL  = [0.80 0.17 0.17];
C_REF    = [0.40 0.40 0.40];
C_FCST   = [0.65 0.65 0.65];
C_CI     = [0.85 0.85 0.85];

FONT     = 'Helvetica';
FS_AX    = 11;
FS_TTL   = 13;
FS_LBL   = 11;
FS_TXT   =  9;
LW_MAIN  = 1.8;
LW_REF   = 1.2;
MS       = 6;

% ── Root defaults (first pass — helps plots created inside subfunctions) ──
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
set(groot,'defaultLegendBox',       'off');
set(groot,'defaultLegendTextColor', 'k');
set(groot,'defaultLegendFontSize',  FS_TXT+1);

% ── Hardened cleanFig ─────────────────────────────────────────────────────
% Explicitly resets every property that MATLAB can override — including
% yyaxis secondary axes, colorbar labels, and sgtitle objects.
function cleanFig(fig)
    set(fig, 'Color','w');

    % All axes (including yyaxis left/right and colorbar axes)
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
            'FontSize',       11);
        ax.Title.Color      = 'k';
        ax.Title.FontName   = 'Helvetica';
        ax.Title.FontSize   = 12;
        ax.Title.FontWeight = 'normal';
        ax.XLabel.Color     = 'k';
        ax.XLabel.FontName  = 'Helvetica';
        ax.YLabel.Color     = 'k';
        ax.YLabel.FontName  = 'Helvetica';
    end

    % All text objects (annotations, data labels, etc.)
    all_txt = findall(fig, 'Type','text');
    set(all_txt, 'Color','k', 'FontName','Helvetica');

    % Legends
    all_lg = findall(fig, 'Type','legend');
    for lg = all_lg'
        set(lg, 'Box','off', 'TextColor','k', ...
            'Color','none', 'EdgeColor','none');
    end

    % Colorbars
    all_cb = findall(fig, 'Type','colorbar');
    for cb = all_cb'
        cb.Color     = 'k';
        cb.FontName  = 'Helvetica';
        cb.FontSize  = 10;
        if ~isempty(cb.Label)
            cb.Label.Color   = 'k';
            cb.Label.FontName = 'Helvetica';
        end
    end

    % sgtitle / subtitle objects live as text children of the figure
    all_fig_txt = findall(fig, 'Type','text');
    set(all_fig_txt, 'Color','k', 'FontName','Helvetica');
end

% ── Axis styling helper (call after all plot commands on an axis) ─────────
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
    ax.FontSize = 11;
    if nargin >= 2 && ~isempty(xl)
        xlabel(ax, xl, 'Color','k', 'FontName','Helvetica', 'FontSize',11);
    end
    if nargin >= 3 && ~isempty(yl)
        ylabel(ax, yl, 'Color','k', 'FontName','Helvetica', 'FontSize',11);
    end
    if nargin >= 4 && ~isempty(ttl)
        title(ax, ttl, 'Color','k', 'FontName','Helvetica', ...
            'FontSize',12, 'FontWeight','normal');
    end
end

function sg = makeSgtitle(txt)
    sg = sgtitle(txt);
    sg.FontName   = 'Helvetica';
    sg.FontSize   = 13;
    sg.FontWeight = 'normal';
    sg.Color      = 'k';
end

% ═══════════════════════════════════════════════════════════════════════════
%  A. LOAD OR SYNTHESISE TIME SERIES
% ═══════════════════════════════════════════════════════════════════════════

if exist('years_all','var')        && exist('annual_mean_hc','var') && ...
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
    
    % Change the end year to 2050
    fcst_years = 2024:2050; 
    n_f        = length(fcst_years); % Calculate length dynamically
    
    % Ensure the decline projects for the full duration (n_f)
    fcst_cover = max(1, obs_cover(end) - 0.55*(1:n_f)); 
    fcst_lo    = max(0,   fcst_cover - 3.5);
    fcst_hi    = min(100, fcst_cover + 3.5);
    source_tag = 'Synthetic (standalone)';
end

all_years = [obs_years,  fcst_years];
all_cover = [obs_cover,  fcst_cover];
n_obs     = length(obs_years);
n_fcst    = length(fcst_years);

% ═══════════════════════════════════════════════════════════════════════════
%  B. PARAMETERS
% ═══════════════════════════════════════════════════════════════════════════
A0         = obs_cover(1);      % Baseline coral cover [%]
T0         = 2.1e6;             % Baseline tourists [persons/yr]
M0         = 105;               % Economic output per tourist [$/person]
M_base     = T0 * M0;           % Baseline annual tourism output [$]

E0         = 1.0e6;             % Offshore wave energy [kJ/m/yr]
k_atten    = 0.004;             % Attenuation coefficient
A_reef_km2 = 9.0;               % Reef area [km²]
A_reef_m2  = A_reef_km2 * 1e6;
C_ref      = A0;
D_max      = 850e6;             % Max coastal defence cost [$]
lambda_dmg = 2.5;

r_disc     = 0.015;              % Social discount rate

% ═══════════════════════════════════════════════════════════════════════════
%  C. COMPUTE LOSSES
% ═══════════════════════════════════════════════════════════════════════════
E_shore   = @(C) E0 .* exp(-k_atten .* C .* (A_reef_m2/1e6));
E_ref_val = E_shore(C_ref);

% Wave
E_t      = E_shore(all_cover);
dE_frac  = max(0, (E_t - E_ref_val) / E_ref_val);
D_wave   = D_max .* (1 - exp(-lambda_dmg .* dE_frac));

% Linear tourism
T_t    = T0 .* (all_cover / A0);
M_t    = T_t .* M0;
dM_t   = max(0, M_base - M_t);

% Total
L_total = D_wave + dM_t;

% Split
D_wave_obs  = D_wave(1:n_obs);    dM_obs  = dM_t(1:n_obs);
D_wave_fcst = D_wave(n_obs+1:end);dM_fcst = dM_t(n_obs+1:end);
T_obs       = T_t(1:n_obs);       T_fcst  = T_t(n_obs+1:end);
L_obs       = L_total(1:n_obs);   L_fcst  = L_total(n_obs+1:end);
dT_obs      = max(0, T0-T_obs);   dT_fcst = max(0, T0-T_fcst);

% ── Uncertainty bounds (physically consistent ordering) ──
Dw_from_lo = D_max .* (1 - exp(-lambda_dmg .* ...
              max(0,(E_shore(fcst_lo)-E_ref_val)/E_ref_val)));
Dw_from_hi = D_max .* (1 - exp(-lambda_dmg .* ...
              max(0,(E_shore(fcst_hi)-E_ref_val)/E_ref_val)));

% Tourism (linear so simpler)
dM_from_lo = max(0, M_base - T0.*(fcst_lo/A0).*M0);
dM_from_hi = max(0, M_base - T0.*(fcst_hi/A0).*M0);

% Total loss bounds (explicit min/max ensures no crossing)
L_lo = min(Dw_from_lo + dM_from_lo, Dw_from_hi + dM_from_hi);
L_hi = max(Dw_from_lo + dM_from_lo, Dw_from_hi + dM_from_hi);

% Also fix wave-only bounds for NPV calc
D_wave_lo = min(Dw_from_lo, Dw_from_hi);
D_wave_hi = max(Dw_from_lo, Dw_from_hi);
dM_lo     = min(dM_from_lo, dM_from_hi);
dM_hi     = max(dM_from_lo, dM_from_hi);

% NPV
t_disc    = (1:n_fcst)';
NPV_wave  = sum(D_wave_fcst(:) ./ (1+r_disc).^t_disc);
NPV_tour  = sum(dM_fcst(:)     ./ (1+r_disc).^t_disc);
NPV_total = NPV_wave + NPV_tour;
NPV_lo    = sum((D_wave_lo(:)+dM_lo(:)) ./ (1+r_disc).^t_disc);
NPV_hi    = sum((D_wave_hi(:)+dM_hi(:)) ./ (1+r_disc).^t_disc);

% ═══════════════════════════════════════════════════════════════════════════
%  D. MITIGATION SCENARIOS
% ═══════════════════════════════════════════════════════════════════════════
cover_sq   = fcst_cover;
cover_hold = repmat(obs_cover(end), 1, n_fcst);
cover_rest = linspace(obs_cover(end), A0, n_fcst);

scenarios      = {cover_sq, cover_hold, cover_rest};
scenario_names = {'Status Quo', 'Hold Cover', 'Restoration'};
npv_scenarios  = zeros(1,3);
for sc = 1:3
    C_sc  = scenarios{sc};
    dE_sc = max(0, (E_shore(C_sc) - E_ref_val) / E_ref_val);
    Dw_sc = D_max .* (1 - exp(-lambda_dmg .* dE_sc));
    dM_sc = max(0, M_base - T0.*(C_sc/A0).*M0);
    npv_scenarios(sc) = sum((Dw_sc(:)+dM_sc(:))./(1+r_disc).^t_disc)/1e6;
end
npv_saved = npv_scenarios(1) - npv_scenarios;

% ═══════════════════════════════════════════════════════════════════════════
%  E. SENSITIVITY  (±50 % on k, D_max, M0)
% ═══════════════════════════════════════════════════════════════════════════
sens_labels = {'k  (attenuation)', 'D_{max}  ($M/yr)', 'M_0  ($/tourist)'};
mults       = [0.5, 1.0, 1.5];
sens_npv    = zeros(3,3);   % rows=low/base/high  cols=k/D_max/M0
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

% ═══════════════════════════════════════════════════════════════════════════
%  F. PRINT
% ═══════════════════════════════════════════════════════════════════════════
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
fprintf('\nMITIGATION\n');
for sc=1:3, fprintf('  %-20s  NPV $%.1fM  (save $%.1fM)\n',...
    scenario_names{sc},npv_scenarios(sc),npv_saved(sc)); end
fprintf('%s\n\n',DIV);

% ═══════════════════════════════════════════════════════════════════════════
%  G. FIGURES  — publication-ready
%
%  Every property is set explicitly after plot(), area(), bar() calls.
%  cleanFig() then does a final sweep to catch anything missed.
%  No style is left to MATLAB defaults or figure inheritance.
% ═══════════════════════════════════════════════════════════════════════════
xlim_all  = [obs_years(1)-0.5, fcst_years(end)+0.5];
xbreak    = fcst_years(1) - 0.5;

% ─────────────────────────────────────────────────────────────────────────
%  Fig E1 — Linear tourism: tourist count and economic output
% ─────────────────────────────────────────────────────────────────────────
figE1 = figure('Position',[60 60 1300 540],'Color','w');

% Left panel — tourists
ax1 = subplot(1,2,1); hold(ax1,'on');
fill([fcst_years, fliplr(fcst_years)], ...
     [T0.*(fcst_hi/A0)/1e6, fliplr(T0.*(fcst_lo/A0)/1e6)], ...
     C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax1, all_years, T_t/1e6, 'o-', ...
    'Color',C_TOUR,'MarkerFaceColor',C_TOUR, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','T(A) = T_0(A/A_0)');
yline(ax1, T0/1e6, '--', 'Color',C_REF,'LineWidth',LW_REF, ...
    'Label',sprintf('Baseline T_0 = %.1fM',T0/1e6), ...
    'LabelHorizontalAlignment','left','HandleVisibility','off');
xline(ax1, xbreak, ':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax1, fcst_years(1)+0.1, T0/1e6*0.97, '\leftarrow forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT);
styleAx(ax1,'Year','Tourists (millions / yr)','T(A) = T_0 \cdot (A / A_0)');
legend(ax1,'Location','southwest');
xlim(ax1, xlim_all); ylim(ax1, [0, T0/1e6*1.18]);

% Right panel — economic output
ax2 = subplot(1,2,2); hold(ax2,'on');
fill([fcst_years, fliplr(fcst_years)], ...
     [T0.*(fcst_hi/A0).*M0/1e6, fliplr(T0.*(fcst_lo/A0).*M0/1e6)], ...
     C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax2, all_years, M_t/1e6, 'o-', ...
    'Color',C_TOUR,'MarkerFaceColor',C_TOUR, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','M(A) = T_0 M_0(A/A_0)');
yline(ax2, M_base/1e6, '--','Color',C_REF,'LineWidth',LW_REF, ...
    'Label',sprintf('M_{base} = $%.1fM',M_base/1e6), ...
    'LabelHorizontalAlignment','left','HandleVisibility','off');
xline(ax2, xbreak, ':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax2, fcst_years(1)+0.1, M_base/1e6*0.97, '\leftarrow forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT);
styleAx(ax2,'Year','Economic Output ($M / yr)','M(T) = T \cdot M_0');
legend(ax2,'Location','southwest');
xlim(ax2, xlim_all); ylim(ax2, [0, M_base/1e6*1.18]);

makeSgtitle('Linear Tourism Model: Visitor Count and Economic Output');
cleanFig(figE1);

% ─────────────────────────────────────────────────────────────────────────
%  Fig E2 — Annual losses: stacked area + component lines
% ─────────────────────────────────────────────────────────────────────────
figE2 = figure('Position',[60 60 1300 560],'Color','w');
y_ceil = max(L_hi/1e6)*1.35 + 2;

% Left panel — stacked area
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
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','Total (forecast)');
xline(ax3, xbreak,':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax3, fcst_years(1)+0.1, y_ceil*0.90, '\leftarrow forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT);
styleAx(ax3,'Year','Annual Economic Loss ($M)','Stacked Annual Losses');
legend(ax3, [ah(1),ah(2)], {'Wave Damage','Tourism Loss'}, 'Location','northwest');
xlim(ax3, xlim_all); ylim(ax3, [0, y_ceil]);

% Right panel — all components as lines
ax4 = subplot(1,2,2); hold(ax4,'on');
fill([fcst_years, fliplr(fcst_years)], [L_hi/1e6, fliplr(L_lo/1e6)], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'DisplayName','95% CI');
plot(ax4, all_years, D_wave/1e6, 's-', ...
    'Color',C_WAVE,'MarkerFaceColor',C_WAVE, ...
    'LineWidth',1.4,'MarkerSize',MS-1,'DisplayName','Wave damage');
plot(ax4, all_years, dM_t/1e6, 'd-', ...
    'Color',C_TOUR,'MarkerFaceColor',C_TOUR, ...
    'LineWidth',1.4,'MarkerSize',MS-1,'DisplayName','Tourism loss');
plot(ax4, all_years, L_total/1e6, 'o-', ...
    'Color',C_TOTAL,'MarkerFaceColor',C_TOTAL, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','Total loss');
xline(ax4, xbreak,':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax4, fcst_years(1)+0.1, y_ceil*0.90, '\leftarrow forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT);
styleAx(ax4,'Year','Annual Economic Loss ($M)','Loss Components');
legend(ax4,'Location','northwest');
xlim(ax4, xlim_all); ylim(ax4, [0, y_ceil]);

makeSgtitle('Coral Reef Economic Losses: Wave Damage + Tourism Revenue');
cleanFig(figE2);

% ─────────────────────────────────────────────────────────────────────────
%  Fig E3 — Transfer function curves
% ─────────────────────────────────────────────────────────────────────────
figE3 = figure('Position',[60 60 1300 520],'Color','w');
cover_range = linspace(0, A0*1.1, 500);
dot_obs  = repmat(C_CORAL, n_obs,  1);
dot_fcst = repmat(C_FCST,  n_fcst, 1);

% T(A)
ax5 = subplot(1,3,1); hold(ax5,'on');
plot(ax5, cover_range, T0.*(cover_range/A0)/1e6, '-', ...
    'Color',C_TOUR,'LineWidth',LW_MAIN+0.2);
xline(ax5, A0,'--','Color',C_REF,'LineWidth',LW_REF);
scatter(ax5, all_cover, T_t/1e6, 22, [dot_obs; dot_fcst], ...
    'filled','MarkerFaceAlpha',0.65);
text(ax5, A0*0.08, T0/1e6*0.78, ...
    sprintf('Slope = T_0/A_0 = %.3fM/%%', T0/A0/1e6), ...
    'FontSize',FS_TXT,'Color',C_TOUR,'FontName',FONT);
styleAx(ax5,'Hard Coral Cover (%)','Tourists (millions/yr)', ...
    'T(A) = T_0 \cdot A / A_0');
ylim(ax5,[0, T0/1e6*1.15]);

% ΔM(A)
ax6 = subplot(1,3,2); hold(ax6,'on');
dM_range = max(0, M_base*(1 - cover_range/A0)) / 1e6;
plot(ax6, cover_range, dM_range, '-', ...
    'Color',C_TOUR,'LineWidth',LW_MAIN+0.2);
xline(ax6, A0,'--','Color',C_REF,'LineWidth',LW_REF);
scatter(ax6, all_cover, dM_t/1e6, 22, [dot_obs; dot_fcst], ...
    'filled','MarkerFaceAlpha',0.65);
text(ax6, A0*0.08, max(dM_range)*0.72, ...
    sprintf('Slope = -M_{base}/A_0\n= -$%.2fM/%%', M_base/A0/1e6), ...
    'FontSize',FS_TXT,'Color',C_TOUR,'FontName',FONT);
styleAx(ax6,'Hard Coral Cover (%)','Tourism Loss ($M/yr)', ...
    '\DeltaM = M_{base}(1 - A/A_0)');

% D_wave(C)
ax7 = subplot(1,3,3); hold(ax7,'on');
dE_range = max(0, (E_shore(cover_range) - E_ref_val) / E_ref_val);
Dw_range = D_max .* (1 - exp(-lambda_dmg .* dE_range)) / 1e6;
plot(ax7, cover_range, Dw_range, '-', ...
    'Color',C_WAVE,'LineWidth',LW_MAIN+0.2);
xline(ax7, A0,'--','Color',C_REF,'LineWidth',LW_REF);
scatter(ax7, all_cover, D_wave/1e6, 22, [dot_obs; dot_fcst], ...
    'filled','MarkerFaceAlpha',0.65);
text(ax7, A0*0.08, max(Dw_range)*0.78, ...
    sprintf('D_{max} = $%.0fM', D_max/1e6), ...
    'FontSize',FS_TXT,'Color',C_WAVE,'FontName',FONT);
styleAx(ax7,'Hard Coral Cover (%)','Wave Damage ($M/yr)', ...
    'D_{wave} = D_{max}[1 - e^{-\lambda\DeltaE}]');

makeSgtitle('Transfer Functions: Coral Cover to Economic Loss');
cleanFig(figE3);

% ─────────────────────────────────────────────────────────────────────────
%  Fig E4 — Cumulative NPV + Mitigation bar chart
% ─────────────────────────────────────────────────────────────────────────
figE4 = figure('Position',[60 60 1300 560],'Color','w');

% Left — cumulative NPV waterfall
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
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','Total NPV');
for i = 1:n_fcst
    text(ax8, fcst_years(i), cum_Dw(i)+cum_dM(i)+max(cum_hi)*0.045, ...
        sprintf('$%.0fM', cum_Dw(i)+cum_dM(i)), ...
        'HorizontalAlignment','center','FontSize',FS_TXT,'Color',C_TOTAL,'FontName',FONT);
end
styleAx(ax8,'Year','Cumulative NPV ($M)', ...
    sprintf('Cumulative Discounted Losses  (r = %.0f%%)',r_disc*100));
legend(ax8, [ac(1),ac(2)], {'Wave NPV','Tourism NPV'},'Location','northwest');
xlim(ax8, [fcst_years(1)-0.5, fcst_years(end)+0.5]);
ylim(ax8, [0, max(cum_hi)*1.25]);

% Right — mitigation bar chart
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
        sprintf('$%.1fM', npv_scenarios(sc)), ...
        'HorizontalAlignment','center','FontSize',FS_TXT,'Color','k','FontName',FONT);
    if sc > 1 && npv_saved(sc) > 0.5
        text(ax9, sc, npv_scenarios(sc)/2, ...
            sprintf('Save $%.1fM', npv_saved(sc)), ...
            'HorizontalAlignment','center','FontSize',FS_TXT-1, ...
            'Color','w','FontWeight','bold','FontName',FONT);
    end
end
yline(ax9, npv_scenarios(1),'--','Color',C_REF,'LineWidth',LW_REF);
set(ax9,'XTick',1:3,'XTickLabel',scenario_names, ...
    'XTickLabelRotation',8,'TickDir','out', ...
    'Box','off','XColor','k','YColor','k', ...
    'FontName',FONT,'FontSize',FS_AX,'Color','w');
ylabel(ax9,'NPV of Total Losses ($M)', ...
    'Color','k','FontName',FONT,'FontSize',FS_LBL);
title(ax9,'Mitigation Scenarios — NPV Comparison', ...
    'Color','k','FontName',FONT,'FontSize',12,'FontWeight','normal');
ylim(ax9, [0, max(npv_scenarios)*1.22]);

makeSgtitle('NPV of Reef-Degradation Losses + Mitigation Value');
cleanFig(figE4);

% ─────────────────────────────────────────────────────────────────────────
%  Fig E5 — Dual-axis: coral cover vs. total annual loss
% ─────────────────────────────────────────────────────────────────────────
figE5 = figure('Position',[60 60 1100 500],'Color','w');
ax10 = axes(figE5); hold(ax10,'on');

yyaxis(ax10,'left');
fill([fcst_years, fliplr(fcst_years)], [fcst_hi, fliplr(fcst_lo)], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax10, obs_years,  obs_cover,'o-', ...
    'Color',C_CORAL,'MarkerFaceColor',C_CORAL, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','Cover (obs)');
plot(ax10, fcst_years, fcst_cover,'o--', ...
    'Color',C_FCST,'MarkerFaceColor','w', ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','Cover (forecast)');
yline(ax10, A0,':','Color',C_REF,'LineWidth',LW_REF);
text(ax10, obs_years(1)+0.2, A0+0.6, ...
    sprintf('A_0 = %.1f%%', A0),'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT);
ax10.YAxis(1).Color = 'k';
ylabel(ax10,'Hard Coral Cover (%)','Color','k','FontName',FONT,'FontSize',FS_LBL);
ylim(ax10, [0, A0*1.65]);

yyaxis(ax10,'right');
fill([fcst_years, fliplr(fcst_years)], [L_hi/1e6, fliplr(L_lo/1e6)], ...
    C_CI,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
plot(ax10, obs_years,  L_obs/1e6,'s-', ...
    'Color',C_TOTAL,'MarkerFaceColor',C_TOTAL, ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','Total loss (obs)');
plot(ax10, fcst_years, L_fcst/1e6,'s--', ...
    'Color',C_TOTAL,'MarkerFaceColor','w', ...
    'LineWidth',LW_MAIN,'MarkerSize',MS,'DisplayName','Total loss (forecast)');
ax10.YAxis(2).Color = 'k';
ylabel(ax10,'Total Annual Loss ($M)','Color','k','FontName',FONT,'FontSize',FS_LBL);

xline(ax10, xbreak,':', 'Color',C_REF,'LineWidth',1,'HandleVisibility','off');
text(ax10, fcst_years(1)+0.1, max(L_hi/1e6)*0.95, '\leftarrow forecast', ...
    'FontSize',FS_TXT,'Color',C_REF,'FontName',FONT);
xlabel(ax10,'Year','Color','k','FontName',FONT,'FontSize',FS_LBL);
lg5 = legend(ax10, ...
    {'Cover (obs)','Cover (forecast)', ...
     'Total loss (obs)','Total loss (forecast)'}, ...
     'Location','northwest');
lg5.Box = 'off'; lg5.TextColor = 'k';
xlim(ax10, xlim_all);

% Explicitly clean both yyaxis sides
ax10.Box      = 'off';
ax10.XGrid    = 'off';  ax10.YGrid = 'off';
ax10.Color    = 'w';
ax10.XColor   = 'k';
ax10.TickDir  = 'out';
ax10.FontName = FONT;
ax10.FontSize = FS_AX;
ax10.LineWidth = 0.9;

% Shift axes down slightly so sgtitle has clear space above
ax10.Position(2) = ax10.Position(2) - 0.04;   % nudge axes down 4%
ax10.Position(4) = ax10.Position(4) - 0.04;   % shrink height to match

makeSgtitle('Coral Cover Decline and Associated Economic Losses');
cleanFig(figE5);

% ═══════════════════════════════════════════════════════════════════════════
%  H. EXPORT CSVs
% ═══════════════════════════════════════════════════════════════════════════
writetable(table(obs_years',obs_cover',T_obs'/1e6,dT_obs', ...
    D_wave_obs'/1e6,dM_obs'/1e6,L_obs'/1e6, ...
    'VariableNames',{'Year','Cover_Pct','Tourists_M','Tourist_Loss', ...
        'WaveDamage_$M','TourismLoss_$M','TotalLoss_$M'}), ...
    'reef_economic_loss_observed.csv');

writetable(table(fcst_years',fcst_cover',T_fcst'/1e6,dT_fcst', ...
    D_wave_fcst'/1e6,dM_fcst'/1e6,L_fcst'/1e6,L_lo'/1e6,L_hi'/1e6, ...
    'VariableNames',{'Year','Cover_Pct','Tourists_M','Tourist_Loss', ...
        'WaveDamage_$M','TourismLoss_$M','TotalLoss_$M', ...
        'TotalLoss_Lower_$M','TotalLoss_Upper_$M'}), ...
    'reef_economic_loss_forecast.csv');

writetable(table({'Wave Protection';'Tourism Revenue';'TOTAL';'CI Low';'CI High'}, ...
    [NPV_wave;NPV_tour;NPV_total;NPV_lo;NPV_hi]/1e6, ...
    'VariableNames',{'Component','NPV_$M'}), 'reef_economic_npv.csv');

writetable(table(scenario_names',npv_scenarios',npv_saved', ...
    'VariableNames',{'Scenario','NPV_Loss_$M','NPV_Saved_$M'}), ...
    'reef_economic_mitigation.csv');

fprintf('CSVs written:\n');
fprintf('  reef_economic_loss_observed.csv\n');
fprintf('  reef_economic_loss_forecast.csv\n');
fprintf('  reef_economic_npv.csv\n');
fprintf('  reef_economic_mitigation.csv\n\n');
fprintf('Figures:  E1 Tourism Model  |  E2 Annual Losses  |  E3 Transfer Functions\n');
fprintf('          E4 NPV + Mitigation  |  E5 Cover vs Loss\n');