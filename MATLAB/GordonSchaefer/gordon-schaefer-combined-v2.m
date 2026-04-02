% Stochastic Gordon-Schaefer Coral Reef Model
% Figure 1: Fish biomass (X/K) — all 7 regions
% Figure 3: NPV degradation premium bar chart — all 7 regions

clear; clc; close all;
rng(42);

set(groot, 'defaultTextInterpreter',              'latex');
set(groot, 'defaultAxesTickLabelInterpreter',     'latex');
set(groot, 'defaultLegendInterpreter',            'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');

% Region definitions

K_FL_base    = 19687072.72 / 1000;
reef_area_FL = 2900;

regions = struct();

regions(1).name            = 'American Samoa';
regions(1).label           = 'Am.\ Samoa';
regions(1).reef_area       = 520;
K_AS                       = K_FL_base * (520 / reef_area_FL);
regions(1).K               = K_AS;
regions(1).X0              = K_AS * (4943.6017 / 11991.3717);
regions(1).r               = 0.35;
regions(1).sigma_X         = 0.3857;
regions(1).Rmax            = 30.9197;
regions(1).R0              = 30.9197;
regions(1).g               = 0.022;
regions(1).delta           = 0.009;
regions(1).sigma_R         = 0.1390;
regions(1).p_shock         = 0.12;
regions(1).reef_shock_mean = 0.28;
regions(1).reef_shock_sd   = 0.10;
regions(1).fish_shock_mean = 0.26;
regions(1).fish_shock_sd   = 0.10;
regions(1).discount_rate   = 0.015;
regions(1).scale_to_global = 284000 / 520;

regions(2).name            = 'CNMI';
regions(2).label           = 'CNMI';
regions(2).reef_area       = 1000;
K_CNMI                     = K_FL_base * (1000 / reef_area_FL);
regions(2).K               = K_CNMI;
regions(2).X0              = K_CNMI * (9043.1394 / 12870.4682);
regions(2).r               = 0.35;
regions(2).sigma_X         = 0.1758;
regions(2).Rmax            = 13.6757;
regions(2).R0              = 11.5481;
regions(2).g               = 0.005;
regions(2).delta           = 0.013;
regions(2).sigma_R         = 0.1460;
regions(2).p_shock         = 0.15;
regions(2).reef_shock_mean = 0.28;
regions(2).reef_shock_sd   = 0.10;
regions(2).fish_shock_mean = 0.25;
regions(2).fish_shock_sd   = 0.10;
regions(2).discount_rate   = 0.015;
regions(2).scale_to_global = 284000 / 1000;

regions(3).name            = 'Florida';
regions(3).label           = 'Florida';
regions(3).reef_area       = reef_area_FL;
regions(3).K               = K_FL_base;
regions(3).X0              = 7731440.87 / 1000;
regions(3).r               = 0.4632;
regions(3).sigma_X         = 0.2161;
regions(3).Rmax            = 10.6953;
regions(3).R0              = 8.1234;
regions(3).g               = 0.05;
regions(3).delta           = 0.04;
regions(3).sigma_R         = 0.3883;
regions(3).p_shock         = 0.04;
regions(3).reef_shock_mean = 0.25;
regions(3).reef_shock_sd   = 0.08;
regions(3).fish_shock_mean = 0.25;
regions(3).fish_shock_sd   = 0.10;
regions(3).discount_rate   = 0.03;
regions(3).scale_to_global = 284000 / 2900;

regions(4).name            = 'Hawaii';
regions(4).label           = 'Hawaii';
regions(4).reef_area       = 1500;
K_HI                       = K_FL_base * (1500 / reef_area_FL);
regions(4).K               = K_HI;
regions(4).X0              = K_HI * (5672.6445 / 19921.8763);
regions(4).r               = 0.35;
regions(4).sigma_X         = 0.5402;
regions(4).Rmax            = 12.1962;
regions(4).R0              = 9.9176;
regions(4).g               = 0.010;
regions(4).delta           = 0.023;
regions(4).sigma_R         = 0.1136;
regions(4).p_shock         = 0.10;
regions(4).reef_shock_mean = 0.22;
regions(4).reef_shock_sd   = 0.08;
regions(4).fish_shock_mean = 0.20;
regions(4).fish_shock_sd   = 0.08;
regions(4).discount_rate   = 0.015;
regions(4).scale_to_global = 284000 / 1500;

regions(5).name            = 'PRIA';
regions(5).label           = 'PRIA';
regions(5).reef_area       = 3000;
K_PRIA                     = K_FL_base * (3000 / reef_area_FL);
regions(5).K               = K_PRIA;
regions(5).X0              = K_PRIA * (13186.4637 / 40989.3037);
regions(5).r               = 0.35;
regions(5).sigma_X         = 0.40;
regions(5).Rmax            = 39.1591;
regions(5).R0              = 39.1591;
regions(5).g               = 0.109;
regions(5).delta           = 0.047;
regions(5).sigma_R         = 0.20;
regions(5).p_shock         = 0.12;
regions(5).reef_shock_mean = 0.30;
regions(5).reef_shock_sd   = 0.10;
regions(5).fish_shock_mean = 0.22;
regions(5).fish_shock_sd   = 0.08;
regions(5).discount_rate   = 0.015;
regions(5).scale_to_global = 284000 / 3000;

regions(6).name            = 'Puerto Rico';
regions(6).label           = 'Puerto Rico';
regions(6).reef_area       = 1790;
K_PR                       = K_FL_base * (1790 / reef_area_FL);
regions(6).K               = K_PR;
regions(6).X0              = K_PR * (485.9264 / 644.5919);
regions(6).r               = 0.35;
regions(6).sigma_X         = 0.3177;
regions(6).Rmax            = 7.6324;
regions(6).R0              = 7.6324;
regions(6).g               = 0.02;
regions(6).delta           = 0.07;
regions(6).sigma_R         = 0.0847;
regions(6).p_shock         = 0.20;
regions(6).reef_shock_mean = 0.30;
regions(6).reef_shock_sd   = 0.10;
regions(6).fish_shock_mean = 0.28;
regions(6).fish_shock_sd   = 0.10;
regions(6).discount_rate   = 0.015;
regions(6).scale_to_global = 284000 / 1790;

regions(7).name            = 'USVI';
regions(7).label           = 'USVI';
regions(7).reef_area       = 230;
K_USVI                     = K_FL_base * (230 / reef_area_FL);
regions(7).K               = K_USVI;
regions(7).X0              = K_USVI * (5527.14 / 7639.84);
regions(7).r               = 0.35;
regions(7).sigma_X         = 0.3652;
regions(7).Rmax            = 9.556;
regions(7).R0              = 9.556;
regions(7).g               = 0.05;
regions(7).delta           = 0.02;
regions(7).sigma_R         = 0.1076;
regions(7).p_shock         = 0.06;
regions(7).reef_shock_mean = 0.25;
regions(7).reef_shock_sd   = 0.09;
regions(7).fish_shock_mean = 0.25;
regions(7).fish_shock_sd   = 0.10;
regions(7).discount_rate   = 0.015;
regions(7).scale_to_global = 284000 / 230;

% Shared simulation parameters

q      = 0.005;
E      = 20;
p_fish = 8;
MC     = 500;
T_sim  = 80;
dt     = 0.05;
N      = round(T_sim / dt);
t      = linspace(0, T_sim, N);
years  = 2023 + t;
nR     = numel(regions);

% Monte Carlo simulation across all regions

results = struct();

for ri = 1:nR
    reg    = regions(ri);
    K      = reg.K;
    X0     = reg.X0;
    Rmax   = reg.Rmax;
    R0     = reg.R0;
    alpha  = K / Rmax;
    collapse_threshold = 0.10 * K;

    NPV_loss_with_deg = zeros(MC,1);
    NPV_loss_no_deg   = zeros(MC,1);
    collapse_count    = 0;
    X_all_deg         = zeros(MC,N);
    X_all_nodeg       = zeros(MC,N);
    R_all             = zeros(MC,N);

    for m = 1:MC
        for scenario = 1:2
            delta_use = reg.delta * (scenario == 1);
            X = zeros(1,N);  R = zeros(1,N);
            X(1) = X0;       R(1) = R0;
            for i = 1:N-1
                K_t         = max(alpha * R(i), 0.02 * K);
                fish_growth = reg.r * X(i) * (1 - X(i)/K_t) - q * E * X(i);
                fish_noise  = reg.sigma_X * X(i) * sqrt(dt) * randn;
                reef_growth = reg.g * R(i) * (1 - R(i)/Rmax) - delta_use * R(i);
                reef_noise  = reg.sigma_R * R(i) * sqrt(dt) * randn;
                fish_shock  = 0;  reef_shock = 0;
                if rand < reg.p_shock * dt
                    fish_shock = max(0, reg.fish_shock_mean + reg.fish_shock_sd * randn) * X(i);
                    reef_shock = max(0, reg.reef_shock_mean + reg.reef_shock_sd * randn) * R(i);
                end
                X(i+1) = min(max(X(i) + fish_growth*dt + fish_noise - fish_shock, 0), K_t);
                R(i+1) = min(max(R(i) + reef_growth*dt + reef_noise - reef_shock, 0), Rmax);
            end
            biomass_loss = max(0, K - X);
            NPV_loss     = sum(p_fish .* biomass_loss .* exp(-reg.discount_rate * t) * dt);
            if scenario == 1
                NPV_loss_with_deg(m) = NPV_loss;
                X_all_deg(m,:)       = X;
                R_all(m,:)           = R;
                if min(X) < collapse_threshold
                    collapse_count = collapse_count + 1;
                end
            else
                NPV_loss_no_deg(m) = NPV_loss;
                X_all_nodeg(m,:)   = X;
            end
        end
    end

    results(ri).X_mean_deg = mean(X_all_deg,   1);
    results(ri).X_mean_nd  = mean(X_all_nodeg, 1);
    results(ri).X_p05      = prctile(X_all_deg,  5, 1);
    results(ri).X_p95      = prctile(X_all_deg, 95, 1);
    results(ri).R_mean     = mean(R_all, 1);
    results(ri).R_p05      = prctile(R_all,  5, 1);
    results(ri).R_p95      = prctile(R_all, 95, 1);
    results(ri).K          = K;
    results(ri).Rmax       = Rmax;
    results(ri).collapse_pct = collapse_count / MC * 100;
    results(ri).npv_premium  = mean(NPV_loss_with_deg - NPV_loss_no_deg);
    results(ri).npv_prem_sd  = std( NPV_loss_with_deg - NPV_loss_no_deg);

    fprintf('%-14s  Collapse: %5.1f%%   NPV premium: $%.3f M\n', ...
        reg.name, collapse_count/MC*100, results(ri).npv_premium/1e6);
end

% Color palette

c_red   = [0.84 0.19 0.15];
c_shade = [0.92 0.92 0.92];

reg_colors = [
    0.13 0.47 0.71;   % AS    — steel blue
    0.55 0.34 0.65;   % CNMI  — purple
    0.17 0.63 0.43;   % FL    — green
    0.95 0.55 0.15;   % HI    — orange
    0.84 0.19 0.15;   % PRIA  — red
    0.10 0.60 0.75;   % PR    — teal
    0.65 0.40 0.20;   % USVI  — brown
];

axisStyle = {'Color','w','FontSize',10, ...
    'XColor',[0.20 0.20 0.20],'YColor',[0.20 0.20 0.20], ...
    'LineWidth',0.6,'TickDir','out','Box','off', ...
    'TickLength',[0.007 0.007],'XTick',2030:10:2100, ...
    'XGrid','off','YGrid','off','TickLabelInterpreter','latex'};

% Figure 1: Fish biomass panel

fig1 = figure('Color','w','Position',[40 40 1400 520]);
tl1  = tiledlayout(1, 1, 'TileSpacing','compact', 'Padding','loose');

axA = nexttile(1);
hold(axA, 'on');

xregion(axA, 2026, 2026.5, 'FaceColor',c_shade,'FaceAlpha',0.7,'EdgeColor','none');

npv_prems   = [results.npv_premium];
[~, idx_hi] = max(npv_prems);
[~, idx_lo] = min(npv_prems);

% Draw 90% CI bands for all regions
for ri = 1:nR
    is_extreme = (ri == idx_hi || ri == idx_lo);
    alpha_fill = 0.18 * is_extreme + 0.07 * ~is_extreme;
    fill(axA, [years, fliplr(years)], ...
        [results(ri).X_p95 ./ results(ri).K, ...
         fliplr(results(ri).X_p05 ./ results(ri).K)], ...
        reg_colors(ri,:), 'FaceAlpha', alpha_fill, 'EdgeColor','none');
end

% Draw mean trajectories
hLines = gobjects(nR,1);
for ri = 1:nR
    is_extreme = (ri == idx_hi || ri == idx_lo);
    lw = 2.4 * is_extreme + 1.6 * ~is_extreme;
    hLines(ri) = plot(axA, years, results(ri).X_mean_deg ./ results(ri).K, ...
        '-', 'Color', reg_colors(ri,:), 'LineWidth', lw);
end

yline(axA, 0.10, ':', 'Color', [c_red 0.65], 'LineWidth', 1.3);
text(axA, 2104, 0.115, '\textit{Collapse} ($10\%\,K$)', ...
    'FontSize',8.5,'Interpreter','latex','Color',c_red,'VerticalAlignment','bottom');

set(axA, axisStyle{:});
axA.XAxis.MinorTick       = 'on';
axA.XAxis.MinorTickValues = 2025:5:2100;
xlim(axA, [2022 2107]);
ylim(axA, [0 1.15]);
ylabel(axA, 'Fish Biomass  $(X / K)$', ...
    'FontSize',11,'Interpreter','latex','Color',[0.15 0.15 0.15]);
xlabel(axA, 'Year', ...
    'FontSize',11,'Interpreter','latex','Color',[0.15 0.15 0.15]);

legend_labels = {regions.name};
hShock = fill(axA, NaN, NaN, c_shade, 'FaceAlpha',0.7,'EdgeColor','none');
hCI    = fill(axA, NaN, NaN, [0.6 0.6 0.6], 'FaceAlpha',0.15,'EdgeColor','none');
legend(axA, [hShock; hCI; hLines], ...
    ['Shock window (2026)'; '$90\%$ CI'; legend_labels(:)], ...
    'FontSize',8.5,'Interpreter','latex', ...
    'EdgeColor',[0.82 0.82 0.82],'Color','w', ...
    'Location','northeast','NumColumns',2);

title(axA, '\textbf{A}\ \ Fish Biomass', ...
    'FontSize',12,'Interpreter','latex','Color',[0.15 0.15 0.15], ...
    'HorizontalAlignment','center');

exportgraphics(fig1, 'figure1_comparative_trajectories.pdf', ...
    'ContentType','vector','BackgroundColor','white');
disp('Saved: figure1_comparative_trajectories.pdf');

% Figure 3: NPV degradation premium bar chart

npv_means     = [results.npv_premium] ./ 1e6;
npv_sds       = [results.npv_prem_sd] ./ 1e6;
collapse_pcts = [results.collapse_pct];
region_labels = {regions.label};

[npv_sorted, sort_idx] = sort(npv_means, 'descend');
npv_sd_sorted   = npv_sds(sort_idx);
collapse_sorted = collapse_pcts(sort_idx);
labels_sorted   = region_labels(sort_idx);
colors_sorted   = reg_colors(sort_idx,:);

fig3 = figure('Color','w','Position',[60 60 1050 560]);
tiledlayout(1,1,'Padding','loose');

axC = nexttile;
hold(axC,'on');

bh = bar(axC, 1:nR, npv_sorted, 0.62, 'FaceColor','flat','EdgeColor','none','BaseValue',0);
for ri = 1:nR
    bh.CData(ri,:) = colors_sorted(ri,:);
end

errorbar(axC, 1:nR, npv_sorted, npv_sd_sorted, ...
    'k.','LineWidth',1.2,'CapSize',6,'MarkerSize',1);

% Collapse probability labels above each bar
for ri = 1:nR
    bar_top = npv_sorted(ri) + npv_sd_sorted(ri);
    text(axC, ri, bar_top + max(npv_sorted)*0.025, ...
        sprintf('$%.0f\\%%$', collapse_sorted(ri)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',9,'Interpreter','latex', ...
        'Color',[0.25 0.25 0.25],'FontWeight','bold');
end

text(axC, 0.97, 0.97, '\textit{\% = collapse probability}', ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',8.5,'Interpreter','latex','Color',[0.40 0.40 0.40]);

yline(axC, 0, '-','Color',[0.6 0.6 0.6],'LineWidth',0.8);

set(axC,'Color','w','FontSize',11, ...
    'XColor',[0.20 0.20 0.20],'YColor',[0.20 0.20 0.20], ...
    'LineWidth',0.6,'TickDir','out','Box','off', ...
    'XTick',1:nR,'XTickLabel',labels_sorted, ...
    'TickLabelInterpreter','latex', ...
    'TickLength',[0.007 0.007], ...
    'XGrid','off','YGrid','on', ...
    'GridColor',[0.88 0.88 0.88],'GridAlpha',1);

ylabel(axC,'NPV Degradation Premium  (\$M)', ...
    'FontSize',12,'Interpreter','latex','Color',[0.15 0.15 0.15]);
title(axC, ...
    ['Economic Cost of Reef Degradation by Region' newline ...
     'NPV Premium (with vs.\ without degradation, $n = 500$)'], ...
    'FontSize',12,'FontWeight','bold','Color',[0.10 0.10 0.10],'Interpreter','latex');

ylim(axC, [0, max(npv_sorted + npv_sd_sorted) * 1.20]);
xlim(axC, [0.35 nR+0.65]);

exportgraphics(fig3, 'figure3_npv_premium.pdf', ...
    'ContentType','vector','BackgroundColor','white');
disp('Saved: figure3_npv_premium.pdf');
