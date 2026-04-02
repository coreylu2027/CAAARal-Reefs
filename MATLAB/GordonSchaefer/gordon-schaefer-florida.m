clear; clc; close all;
rng(42);

set(groot, ...
    'defaultTextInterpreter',              'latex', ...
    'defaultAxesTickLabelInterpreter',     'latex', ...
    'defaultLegendInterpreter',            'latex', ...
    'defaultColorbarTickLabelInterpreter', 'latex', ...
    'defaultAxesFontSize',                 13, ...
    'defaultTextFontSize',                 13, ...
    'defaultAxesLineWidth',                0.8, ...
    'defaultLineLineWidth',                1.8, ...
    'defaultAxesBox',                      'off', ...
    'defaultAxesTickDir',                  'out', ...
    'defaultFigureColor',                  'w');

K       = 19687.1;
X0      = 7731.4;
r       = 0.46;
q       = 0.0003;
E       = 120;
sigma_X = 0.2161;

Rmax    = 10.6953;
R0      = 8.1234;
g       = 0.05;
delta   = 0.04;
sigma_R = 0.3883;
alpha   = K / Rmax;

p_shock         = 0.04;
reef_shock_mean = 0.25;
reef_shock_sd   = 0.08;
fish_shock_mean = 0.25;
fish_shock_sd   = 0.10;

p_fish        = 8;
discount_rate = 0.015;
scale_to_global = 284000 / 2900;
Global_Pop      = 8e9;

T_sim = 80;
dt    = 0.05;
N     = round(T_sim / dt);
t     = linspace(0, T_sim, N);
years = 2022 + t;
MC    = 500;

collapse_threshold = 0.10 * K;

NPV_loss_with_deg = zeros(MC, 1);
NPV_loss_no_deg   = zeros(MC, 1);
collapse_count    = 0;
X_all_deg         = zeros(MC, N);
X_all_nodeg       = zeros(MC, N);
R_all             = zeros(MC, N);

% scenario 1 = with degradation, scenario 2 = no degradation
for m = 1:MC
    for scenario = 1:2
        delta_use = delta * (scenario == 1);

        X    = zeros(1, N);
        R    = zeros(1, N);
        X(1) = X0;
        R(1) = R0;
        npv  = 0;

        for i = 1:N-1
            K_t = max(alpha * R(i), 0.02 * K);

            fish_growth = r * X(i) * (1 - X(i)/K_t) - q*E*X(i);
            fish_noise  = sigma_X * X(i) * sqrt(dt) * randn;

            reef_growth = g * R(i) * (1 - R(i)/Rmax) - delta_use * R(i);
            reef_noise  = sigma_R * R(i) * sqrt(dt) * randn;

            fish_shock = 0;
            reef_shock = 0;
            if rand < p_shock * dt
                fish_shock = -abs(normrnd(fish_shock_mean, fish_shock_sd)) * X(i);
                reef_shock = -abs(normrnd(reef_shock_mean, reef_shock_sd)) * R(i);
            end

            X(i+1) = max(0, min(X(i) + fish_growth*dt + fish_noise + fish_shock, K_t));
            R(i+1) = max(0, min(R(i) + reef_growth*dt + reef_noise + reef_shock, Rmax));

            npv = npv + p_fish * max(0, K - X(i+1)) * exp(-discount_rate * t(i)) * dt;
        end

        if scenario == 1
            X_all_deg(m,:)       = X;
            R_all(m,:)           = R;
            NPV_loss_with_deg(m) = npv;
            if X(end) < collapse_threshold
                collapse_count = collapse_count + 1;
            end
        else
            X_all_nodeg(m,:)   = X;
            NPV_loss_no_deg(m) = npv;
        end
    end
end

X_mean_deg = mean(X_all_deg,   1);
X_mean_nd  = mean(X_all_nodeg, 1);

% CI over non-collapsed trajectories only
non_collapsed = X_all_deg(:,end) > 0.01 * K;
X_ci = X_all_deg(non_collapsed, :);
if isempty(X_ci), X_ci = X_all_deg; end
X_upper = prctile(X_ci, 95, 1);
X_lower = prctile(X_ci,  5, 1);

collapse_by_year = mean(X_all_deg < collapse_threshold, 1);

obs_years   = [2014, 2016, 2018, 2021, 2022] - 2022;
obs_biomass = [19687072.72, 15833018.77, 10143101.85, 7288779.06, 7731440.87] / 1000;

fprintf('\n============ RESULTS ============\n');
fprintf('Mean PV loss (w/ degradation):    $%.2f\n',  mean(NPV_loss_with_deg));
fprintf('Mean PV loss (no degradation):    $%.2f\n',  mean(NPV_loss_no_deg));
fprintf('Additional loss from degradation: $%.2f\n',  mean(NPV_loss_with_deg - NPV_loss_no_deg));
fprintf('Degradation share of total loss:  %.2f%%\n', ...
    mean(NPV_loss_with_deg - NPV_loss_no_deg) / mean(NPV_loss_with_deg) * 100);
fprintf('Collapse probability:             %.2f%%\n', collapse_count/MC*100);
fprintf('Global loss from degradation:     $%.2f M\n', ...
    mean(NPV_loss_with_deg - NPV_loss_no_deg) * scale_to_global / 1e6);
fprintf('Global total biomass loss value:  $%.2f M\n', ...
    mean(NPV_loss_with_deg) * scale_to_global / 1e6);
fprintf('Per capita global impact:         $%.4f\n', ...
    mean(NPV_loss_with_deg - NPV_loss_no_deg) * scale_to_global / Global_Pop);

c_blue   = [0.17 0.45 0.72];
c_black  = [0.12 0.12 0.12];
c_green  = [0.10 0.55 0.25];
c_red    = [0.80 0.10 0.10];
c_orange = [0.90 0.40 0.10];
c_grey   = [0.55 0.55 0.55];

fig = figure('Position', [80 80 1120 940], 'Color', 'w');

% Panel 1: Fish Biomass
ax1 = subplot(2, 2, [1 2]);
fill(ax1, [t, fliplr(t)], [X_upper, fliplr(X_lower)], ...
    c_blue, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', '$90\%$ CI');
hold(ax1, 'on');
plot(ax1, t, X_mean_deg, '-',  'Color', c_blue,  'LineWidth', 2.2, 'DisplayName', 'Mean --- with degradation');
plot(ax1, t, X_mean_nd,  '--', 'Color', c_black, 'LineWidth', 1.8, 'DisplayName', 'Mean --- no degradation');
plot(ax1, obs_years, obs_biomass, 'o', 'Color', c_red, 'MarkerFaceColor', c_red, ...
    'MarkerSize', 7, 'DisplayName', 'Observed (2014--2022)');
yline(ax1, K, '--', 'Color', c_green, 'LineWidth', 1.4, ...
    'Label', '$K$', 'LabelHorizontalAlignment', 'right', 'Interpreter', 'latex', 'FontSize', 11);
yline(ax1, collapse_threshold, ':', 'Color', c_red, 'LineWidth', 1.4, ...
    'Label', 'Collapse ($10\%\,K$)', 'LabelHorizontalAlignment', 'right', ...
    'Interpreter', 'latex', 'FontSize', 11);
ylabel(ax1, 'Fish Biomass (kg)');
xlabel(ax1, 'Years from 2022');
title(ax1, '\textbf{Fish Biomass Projection} --- Florida Keys / Dry Tortugas ($n = 500$)', 'FontSize', 14);
legend(ax1, 'Location', 'northeast', 'FontSize', 11, 'Box', 'off');
grid(ax1, 'on');  ax1.GridAlpha = 0.15;  ax1.GridLineStyle = ':';

% Panel 2: Collapse Probability Over Time
ax2 = subplot(2, 2, 3);
fill(ax2, [years, fliplr(years)], [collapse_by_year*100, zeros(1,length(years))], ...
    c_red, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
hold(ax2, 'on');
plot(ax2, years, collapse_by_year*100, '-', 'Color', c_red, 'LineWidth', 2.2);
yline(ax2, 50, '--', 'Color', c_grey, 'LineWidth', 1.2, ...
    'Label', '$P = 0.50$', 'LabelHorizontalAlignment', 'right', ...
    'Interpreter', 'latex', 'FontSize', 10);
xlabel(ax2, 'Year');
ylabel(ax2, '$P(\mathrm{collapse})$ (\%)');
title(ax2, '\textbf{Collapse Probability Over Time}', 'FontSize', 14);
ylim(ax2, [0 100]);
xlim(ax2, [years(1) years(end)]);
grid(ax2, 'on');  ax2.GridAlpha = 0.15;  ax2.GridLineStyle = ':';

% Panel 3: NPV Loss Distribution
ax3 = subplot(2, 2, 4);
deg_loss = (NPV_loss_with_deg - NPV_loss_no_deg) / 1e3;
histogram(ax3, deg_loss, 40, 'FaceColor', c_orange, 'EdgeColor', 'w', 'FaceAlpha', 0.88);
hold(ax3, 'on');
mu = mean(deg_loss);
xline(ax3, mu, '--', 'Color', c_black, 'LineWidth', 2.0, ...
    'Label', sprintf('$\\bar{x} = %.1f$k', mu), ...
    'Interpreter', 'latex', 'FontSize', 11, 'LabelHorizontalAlignment', 'right');
xlabel(ax3, 'Degradation Loss (thousand USD)');
ylabel(ax3, 'Frequency');
title(ax3, '\textbf{PV Loss from Reef Degradation}', 'FontSize', 14);
grid(ax3, 'on');  ax3.GridAlpha = 0.15;  ax3.GridLineStyle = ':';

sgtitle(fig, ...
    ['\textbf{Stochastic Gordon--Schaefer Model:} ' ...
     'Florida Keys / Dry Tortugas ($n = 500$ Monte Carlo runs)'], ...
    'FontSize', 16, 'Interpreter', 'latex');

exportgraphics(fig, 'florida_gs.pdf', 'ContentType', 'vector', 'BackgroundColor', 'white');
