% Pharmaceutical Option Value Model for Reef Biodiversity Loss
% Based on Simpson et al. (1996), Costello & Ward (2006)

clear; clc; close all;

set(groot, 'defaultTextInterpreter',              'latex');
set(groot, 'defaultAxesTickLabelInterpreter',     'latex');
set(groot, 'defaultLegendInterpreter',            'latex');
set(groot, 'defaultAxesFontSize', 11);

% Region parameters

regions = {'Florida Keys', 'Puerto Rico', 'U.S.\ Virgin Islands', ...
           'American Samoa', 'Mariana Islands', 'Hawaii', ...
           'Pacific Islands (PRIA)'};
nRegions = numel(regions);

reef_area_km2 = [9.0; 6.3; 2.8; 8.1; 5.9; 11.2; 41.0];

cover_baseline = [13.49; 8.20; 8.10; 14.50; 10.80; 12.30; 20.00];

% Species scaling constants calibrated to NOAA RVC richness surveys
c_scaling = [310; 280; 190; 340; 295; 325; 410];

% Time grid 2024-2100

t_start = 2024;
t_end   = 2100;
t_vec   = (t_start : t_end)';
nT      = numel(t_vec);

% SARIMAX coral cover projections from Table 9 anchor points

anchor_yr  = [2024, 2025, 2026, 2027, 2028, 2029, 2030, 2032, ...
              2036, 2040, 2044, 2048, 2049];
anchor_cov = [5.32, 6.18, 5.58, 4.74, 2.97, 3.95, 3.48, 5.62, ...
              1.79, 4.36, 0.58, 3.09, 1.89];
anchor_lo  = [3.67, 3.41, 2.69, 1.84, 0.07, 1.04, 0.58, 2.72, ...
              0.00, 1.45, 0.00, 0.19, 0.00];
anchor_hi  = [6.97, 8.96, 8.47, 7.64, 5.88, 6.85, 6.38, 8.52, ...
              4.69, 7.26, 3.49, 6.00, 4.80];

cover_mid_2049 = interp1(anchor_yr, anchor_cov, t_start:2049, 'pchip');
cover_lo_2049  = interp1(anchor_yr, anchor_lo,  t_start:2049, 'pchip');
cover_hi_2049  = interp1(anchor_yr, anchor_hi,  t_start:2049, 'pchip');

trend        = -0.6397;
t_extra      = (2050:t_end)';
cov_2049_end = cover_mid_2049(end);

cover_mid_extra = max(0, cov_2049_end + trend*(t_extra - 2049));
cover_lo_extra  = max(0, cover_lo_2049(end)  + trend*(t_extra - 2049));
cover_hi_extra  = max(0, cover_hi_2049(end)  + trend*(t_extra - 2049));

cover_mid = max(0, [cover_mid_2049, cover_mid_extra'])';
cover_lo  = max(0, [cover_lo_2049,  cover_lo_extra']');
cover_hi  = max(0, [cover_hi_2049,  cover_hi_extra']');

% Species-area curve S(A,t) = c * A(t)^z

z     = 0.25;
C_max = 27.52;

S_mat = zeros(nT, nRegions);
for rr = 1:nRegions
    A_t = reef_area_km2(rr) .* max(0, cover_mid ./ C_max);
    S_mat(:, rr) = c_scaling(rr) .* max(0, A_t).^z;
end
S_global = sum(S_mat, 2);

% Bioeconomic prospecting — lifetime NPV drug values

p_discovery = 1e-4;

rho      = 0.015;
pi_k     = [0.05,  0.20,  0.25,  0.50];
R_annual = [400e6, 75e6,  175e6, 30e6];
L_k      = [20,    15,    15,    30  ];
tau_k    = [12,    15,    14,    10  ];

% Annuity factor: PV of $1/yr for L_k years discounted back tau_k years to discovery
annuity_k    = (1 - exp(-rho .* L_k)) ./ rho .* exp(-rho .* tau_k);
v_lifetime_k = R_annual .* annuity_k;
v_compound   = sum(pi_k .* v_lifetime_k);
pv_per_species = p_discovery * v_compound;

fprintf('=== Bioeconomic Parameters ===\n');
fprintf('  Discovery probability p             : %.2e\n', p_discovery);
fprintf('\n  Outcome-level lifetime NPVs:\n');
outcome_labels = {'Blockbuster (AZT-like)', ...
                  'Niche (Ziconotide)', ...
                  'Moderate (Trabectedin)', ...
                  'Generic/societal (Cytarabine)'};
for kk = 1:4
    fprintf('    k=%d  %-28s pi=%.2f  annual=$%.0fM  L=%d yr  tau=%d yr  -> lifetime NPV=$%.1fB\n', ...
        kk, outcome_labels{kk}, pi_k(kk), R_annual(kk)/1e6, L_k(kk), tau_k(kk), ...
        v_lifetime_k(kk)/1e9);
end
fprintf('\n  E[lifetime NPV per compound] v      : $%.3fB\n', v_compound/1e9);
fprintf('  E[V(s)] = p * v per species         : $%.4fM\n', pv_per_species/1e6);

% Cumulative pharmaceutical option value loss (triple integral)

disc = exp(-rho * (t_vec - t_start));

pharma_value_undiscounted = zeros(nT, nRegions);
for rr = 1:nRegions
    A_t = reef_area_km2(rr) .* max(0, cover_mid ./ C_max);
    pharma_value_undiscounted(:, rr) = p_discovery * v_compound * ...
        c_scaling(rr) / (z+1) .* max(0, A_t).^(z+1);
end

pharma_baseline   = pharma_value_undiscounted(1, :);
pharma_loss_annual = max(0, pharma_baseline - pharma_value_undiscounted);

pharma_loss_discounted = pharma_loss_annual .* disc;
cumulative_NPV_loss = zeros(nT, nRegions);
for rr = 1:nRegions
    for tt = 2:nT
        cumulative_NPV_loss(tt, rr) = trapz(t_vec(1:tt), ...
            pharma_loss_discounted(1:tt, rr));
    end
end
total_NPV_loss = sum(cumulative_NPV_loss, 2);

fprintf('\n=== Pharmaceutical NPV Loss (Status Quo, 2024--2100) ===\n');
for rr = 1:nRegions
    fprintf('  %-32s: $%.3fB\n', regions{rr}, cumulative_NPV_loss(end, rr)/1e9);
end
fprintf('  %-32s: $%.3fB\n', 'GLOBAL TOTAL', total_NPV_loss(end)/1e9);

% Episodic bleaching — permanent loss at trough years

bleach_years = [2028, 2036, 2044, 2049];
bleach_idx   = bleach_years - t_start + 1;

species_at_trough  = S_global(bleach_idx);
species_prior      = S_global(max(1, bleach_idx - 1));
species_lost_event = max(0, species_prior - species_at_trough);
permanent_loss     = species_lost_event * pv_per_species;

fprintf('\n=== Permanent Option Value at Bleaching Troughs ===\n');
for kk = 1:numel(bleach_years)
    fprintf('  %d: ~%.0f species lost,  permanent loss ~$%.3fB\n', ...
            bleach_years(kk), species_lost_event(kk), permanent_loss(kk)/1e9);
end

% Sensitivity analysis — discount rates 1.5%, 3%, 7%

disc_rates  = [0.015, 0.03, 0.07];
NPV_by_rate = zeros(1, 3);
for dd = 1:3
    disc_d = exp(-disc_rates(dd) * (t_vec - t_start));
    for rr = 1:nRegions
        A_t = reef_area_km2(rr) .* max(0, cover_mid ./ C_max);
        loss_ann = max(0, pharma_baseline(rr) - ...
            p_discovery * v_compound * c_scaling(rr) / (z+1) .* ...
            max(0, A_t).^(z+1));
        NPV_by_rate(dd) = NPV_by_rate(dd) + trapz(t_vec, loss_ann .* disc_d);
    end
end

fprintf('\n=== Sensitivity: Total NPV Loss by Discount Rate ===\n');
disc_label_str = {'1.5%', '3.0%', '7.0%'};
for dd = 1:3
    fprintf('  r = %s:  $%.3fB\n', disc_label_str{dd}, NPV_by_rate(dd)/1e9);
end

% Color definitions

col_coral  = [0.85, 0.33, 0.10];
col_algae  = [0.18, 0.55, 0.34];
col_pharma = [0.27, 0.45, 0.77];
col_ci     = [0.70, 0.80, 0.93];
col_episod = [0.93, 0.69, 0.13];
col_perm   = [0.64, 0.08, 0.18];
colors_r   = lines(nRegions);

% Figure 1: Hard coral cover trajectory

fig1 = figure('Units','inches','Position',[1 1 7 4]);
hold on;
fill([t_vec; flipud(t_vec)], [cover_hi; flipud(cover_lo)], ...
     col_ci, 'EdgeColor','none','FaceAlpha',0.5);
plot(t_vec, cover_mid, '-', 'Color', col_coral, 'LineWidth', 1.8);
for kk = 1:numel(bleach_years)
    xline(bleach_years(kk), '--', 'Color', col_episod, 'LineWidth', 1.2);
end
yline(10, ':', 'Color',[0.3 0.3 0.3], 'LineWidth', 1.0);
xlabel('Year');
ylabel('Hard Coral Cover (\%)');
title({'SARMA$(1,1)\times(1,1)_4$ Hard Coral Cover Forecast', ...
       'with 95\% Prediction Interval (2024--2100)'});
legend({'95\% PI', 'SARMA Forecast', 'Episodic Trough Years', ...
        'Critical Threshold (10\%)'}, 'Location','northeast','FontSize',9);
xlim([t_start t_end]); ylim([0 12]);
grid on; box on;
hold off;

% Figure 2: Global species count stacked by region

fig2 = figure('Units','inches','Position',[1 1 7 4]);
hold on;
h_area = area(t_vec, S_mat);
for rr = 1:nRegions
    h_area(rr).FaceColor = colors_r(rr,:);
    h_area(rr).EdgeColor = 'none';
    h_area(rr).FaceAlpha = 0.75;
end
for kk = 1:numel(bleach_years)
    xline(bleach_years(kk), '--', 'Color', col_episod, 'LineWidth', 1.1);
end
xlabel('Year');
ylabel('Estimated Species Count $S(A,t)$');
title({'Reef-Associated Species Count via Species--Area Relationship', ...
       '$S(A,t) = c \cdot A(t)^{z}$, $z = 0.25$ (stacked by region)'});
legend(regions, 'Location','northeast','FontSize',7.5);
xlim([t_start t_end]); grid on; box on;
hold off;

% Figure 3: Annual pharmaceutical option value by region

fig3 = figure('Units','inches','Position',[1 1 7 4]);
hold on;
for rr = 1:nRegions
    plot(t_vec, pharma_value_undiscounted(:,rr)/1e9, ...
         'LineWidth', 1.4, 'Color', colors_r(rr,:));
end
for kk = 1:numel(bleach_years)
    xline(bleach_years(kk), '--', 'Color', col_episod, 'LineWidth', 1.0);
end
xlabel('Year');
ylabel('Annual Pharma Option Value (\$B/yr)');
title({'Annual Pharmaceutical Option Value by Region', ...
       '$\mathbb{E}[V(A,t)] = p \cdot v_{\mathrm{lifetime}} \cdot c/(z+1) \cdot A(t)^{z+1}$'});
legend(regions, 'Location','northeast','FontSize',7.5);
xlim([t_start t_end]); grid on; box on;
hold off;

% Figure 4: Cumulative NPV loss and discount rate sensitivity

fig4 = figure('Units','inches','Position',[1 1 12 5]);

subplot(1,2,1); hold on;
for rr = 1:nRegions
    plot(t_vec, cumulative_NPV_loss(:,rr)/1e9, ...
         'LineWidth', 1.4, 'Color', colors_r(rr,:));
end
plot(t_vec, total_NPV_loss/1e9, 'k-', 'LineWidth', 2.2);
for kk = 1:numel(bleach_years)
    xline(bleach_years(kk), '--', 'Color', col_episod, 'LineWidth', 1.0);
end
xlabel('Year');
ylabel('Cumulative NPV Loss (\$B)');
title({'Cumulative Discounted Pharmaceutical', ...
       'Option Value Loss ($r=1.5\%$)'});
legend([regions, {'Global Total'}], 'Location','northwest','FontSize',7.0);
xlim([t_start t_end]); grid on; box on; hold off;

subplot(1,2,2); hold on;
disc_labels = {'$r = 1.5\%$','$r = 3.0\%$','$r = 7.0\%$'};
disc_styles = {'-','--',':'};
disc_colors = {col_pharma, col_coral, col_algae};
for dd = 1:3
    disc_d = exp(-disc_rates(dd) * (t_vec - t_start));
    total_loss_d = zeros(nT, 1);
    for rr = 1:nRegions
        A_t = reef_area_km2(rr) .* max(0, cover_mid ./ C_max);
        loss_ann = max(0, pharma_baseline(rr) - ...
            p_discovery * v_compound * c_scaling(rr) / (z+1) .* ...
            max(0, A_t).^(z+1));
        cum_d = zeros(nT,1);
        for tt = 2:nT
            cum_d(tt) = trapz(t_vec(1:tt), loss_ann(1:tt) .* disc_d(1:tt));
        end
        total_loss_d = total_loss_d + cum_d;
    end
    plot(t_vec, total_loss_d/1e9, disc_styles{dd}, ...
         'Color', disc_colors{dd}, 'LineWidth', 2.0, ...
         'DisplayName', disc_labels{dd});
end
for kk = 1:numel(bleach_years)
    xline(bleach_years(kk), '--', 'Color', col_episod, 'LineWidth', 1.0);
end
xlabel('Year');
ylabel('Global Cumulative NPV Loss (\$B)');
title({'Sensitivity: Global NPV Loss', 'by Discount Rate'});
legend('Location','northwest','FontSize',9);
xlim([t_start t_end]); grid on; box on; hold off;

sgtitle({'\textbf{Pharmaceutical Option Value -- Cumulative Loss}', ...
         'Lifetime NPV Drug Values Applied'}, ...
         'Interpreter','latex','FontSize',12);

% Figure 5: Permanent loss bar chart at bleaching troughs

fig5 = figure('Units','inches','Position',[1 1 7 4.5]);
hold on;
perm_frac = 0.15;
annual_loss_trough = zeros(1, numel(bleach_years));
for kk = 1:numel(bleach_years)
    tt = bleach_idx(kk);
    annual_loss_trough(kk) = sum(pharma_loss_annual(tt,:));
end
perm_component  = perm_frac * annual_loss_trough / 1e9;
recov_component = (1 - perm_frac) * annual_loss_trough / 1e9;

b = bar(bleach_years, [recov_component; perm_component]', 'stacked');
b(1).FaceColor = col_pharma;  b(1).EdgeColor = 'none';
b(2).FaceColor = col_perm;    b(2).EdgeColor = 'none';
xlabel('Bleaching Trough Year');
ylabel('Annual Option Value Loss (\$B/yr)');
title({'\textbf{Option Value Lost at Episodic Bleaching Troughs}', ...
       'Permanent vs.\ Potentially Recoverable'});
legend({'Potentially Recoverable','Permanent (pre-discovery extinction)'}, ...
       'Location','northwest','FontSize',9);
xticks(bleach_years);
for kk = 1:numel(bleach_years)
    text(bleach_years(kk), annual_loss_trough(kk)/1e9 + 0.01, ...
         sprintf('\\$%.3fB', annual_loss_trough(kk)/1e9), ...
         'HorizontalAlignment','center','FontSize',8.5,'Interpreter','latex');
end
grid on; box on;
hold off;

% Figure 6: Triple-integral integrand surface

fig6 = figure('Units','inches','Position',[1 1 7 4]);
t_mesh = linspace(t_start, t_end, 80);
a_frac = linspace(0, 1, 60);
[T_m, A_m] = meshgrid(t_mesh, a_frac);

c_global = sum(c_scaling);
A_total  = sum(reef_area_km2);
C_t_mesh = max(0, interp1(t_vec, cover_mid, t_mesh, 'linear','extrap'));
A_t_frac = max(0, C_t_mesh / C_max);

A_eff          = A_m .* A_total .* repmat(A_t_frac, size(a_frac,2), 1);
integrand      = p_discovery * v_compound * c_global .* max(0, A_eff).^z;
disc_mesh      = exp(-rho * (T_m - t_start));
integrand_disc = integrand .* disc_mesh / 1e9;

surf(T_m, A_m, integrand_disc, 'EdgeColor','none');
colormap(parula);
cb = colorbar;
cb.Label.String      = '\$B / (yr $\cdot$ area fraction)';
cb.Label.Interpreter = 'latex';
xlabel('Year $t$');
ylabel('Reef Area Fraction $A / A_{\max}$');
zlabel('Discounted Integrand (\$B)');
title({'\textbf{Triple-Integral Integrand Surface}', ...
       '$\mathcal{L}_{\mathrm{pharma}} = \int_0^T \int_0^{A(t)} \int_0^{S(A,t)} p(s)\,v_{\mathrm{lifetime}}(s)\,ds\,dA\,dt$'});
view(45, 30); grid on; box on;

% Export figures

exportgraphics(fig1, 'pharma_fig1_cover_forecast.png',  'Resolution', 200);
exportgraphics(fig2, 'pharma_fig2_species_count.png',   'Resolution', 200);
exportgraphics(fig3, 'pharma_fig3_annual_value.png',    'Resolution', 200);
exportgraphics(fig4, 'pharma_fig4_cumulative_npv.png',  'Resolution', 200);
exportgraphics(fig5, 'pharma_fig5_trough_loss.png',     'Resolution', 200);
exportgraphics(fig6, 'pharma_fig6_triple_integral.png', 'Resolution', 200);

fprintf('\nAll figures exported.\n');
fprintf('\n=== FINAL SUMMARY ===\n');
fprintf('  Global pharmaceutical NPV loss (r=1.5%%): $%.2fB\n', total_NPV_loss(end)/1e9);
fprintf('  Global pharmaceutical NPV loss (r=3.0%%): $%.2fB\n', NPV_by_rate(2)/1e9);
fprintf('  Global pharmaceutical NPV loss (r=7.0%%): $%.2fB\n', NPV_by_rate(3)/1e9);
