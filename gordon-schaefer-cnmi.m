clear; clc; close all;
rng(42);

% K scaled from Florida anchor by reef area ratio
reef_area_MA = 1000;
reef_area_FL = 2900;
K_FL         = 19687072.72/1000;
K            = K_FL * (reef_area_MA / reef_area_FL);
X0_raw = 9043.1394; K_raw = 12870.4682;
X0     = K * (X0_raw / K_raw);

r       = 0.35;
q       = 0.005;
sigma_X = 0.1758;
E       = 20;

Rmax    = 13.6757;
R0      = 11.5481;
g       = 0.005;
delta   = 0.013;
sigma_R = 0.1460;
alpha   = K / Rmax;

p_shock         = 0.15;
reef_shock_mean = 0.28;
reef_shock_sd   = 0.10;
fish_shock_mean = 0.25;
fish_shock_sd   = 0.10;

p_fish        = 8;
discount_rate = 0.015;
scale_to_global = 284000 / reef_area_MA;
Global_Pop      = 8e9;

T_sim = 80; dt = 0.05; N = round(T_sim/dt);
t = linspace(0,T_sim,N); years = 2022+t; MC = 500;
collapse_threshold = 0.10*K;

NPV_loss_with_deg = zeros(MC,1);
NPV_loss_no_deg   = zeros(MC,1);
collapse_count    = 0;
X_all_deg         = zeros(MC,N);
X_all_nodeg       = zeros(MC,N);
R_all             = zeros(MC,N);

% scenario 1 = with degradation, scenario 2 = no degradation
for m = 1:MC
    for scenario = 1:2
        delta_use = delta * (scenario == 1);
        X = zeros(1,N); R = zeros(1,N); X(1) = X0; R(1) = R0;
        for i = 1:N-1
            K_t = max(alpha*R(i), 0.02*K);
            fish_growth = r*X(i)*(1-X(i)/K_t) - q*E*X(i);
            fish_noise  = sigma_X*X(i)*sqrt(dt)*randn;
            reef_growth = g*R(i)*(1-R(i)/Rmax) - delta_use*R(i);
            reef_noise  = sigma_R*R(i)*sqrt(dt)*randn;
            fish_shock = 0; reef_shock = 0;
            if rand < p_shock*dt
                fish_shock = max(0, fish_shock_mean + fish_shock_sd*randn)*X(i);
                reef_shock = max(0, reef_shock_mean + reef_shock_sd*randn)*R(i);
            end
            X(i+1) = max(X(i) + fish_growth*dt + fish_noise - fish_shock, 0);
            R(i+1) = max(R(i) + reef_growth*dt + reef_noise - reef_shock, 0);
        end
        biomass_loss = max(0, K-X);
        NPV_loss = sum(p_fish .* biomass_loss .* exp(-discount_rate*t) * dt);
        if scenario == 1
            NPV_loss_with_deg(m) = NPV_loss; X_all_deg(m,:) = X; R_all(m,:) = R;
            if min(X) < collapse_threshold, collapse_count = collapse_count+1; end
        else
            NPV_loss_no_deg(m) = NPV_loss; X_all_nodeg(m,:) = X;
        end
    end
end

disp('--- RESULTS ---')
disp(['K (kg):                    ' num2str(K)])
disp(['X0 (kg):                   ' num2str(X0)])
disp(['Total NPV loss (w/ deg):   $' num2str(mean(NPV_loss_with_deg)/1e6, '%.2f') ' M'])
disp(['Total NPV loss (no deg):   $' num2str(mean(NPV_loss_no_deg)/1e6,   '%.2f') ' M'])
disp(['Degradation loss:          $' num2str(mean(NPV_loss_with_deg - NPV_loss_no_deg)/1e6, '%.2f') ' M'])
disp(['Global degradation loss:   $' num2str(mean(NPV_loss_with_deg - NPV_loss_no_deg)*scale_to_global/1e6, '%.2f') ' M'])
disp(['Collapse probability:      ' num2str(collapse_count/MC*100, '%.1f') '%'])

X_mean_deg  = mean(X_all_deg,   1);
X_upper_deg = prctile(X_all_deg, 95, 1);
X_lower_deg = prctile(X_all_deg,  5, 1);
X_mean_nd   = mean(X_all_nodeg, 1);

figure('Color','w','Position',[100 100 1400 600]);
fill([years,fliplr(years)],[X_upper_deg,fliplr(X_lower_deg)], ...
    [0.7 0.85 1],'EdgeColor','none','DisplayName','90% CI'); hold on
plot(years, X_mean_deg, 'b-',  'LineWidth',2,'DisplayName','Mean (with degradation)')
plot(years, X_mean_nd,  'k--', 'LineWidth',2,'DisplayName','Mean (no degradation)')
plot(years(1), X0, 'ro','MarkerSize',5,'MarkerFaceColor','r','DisplayName','Observed 2022')
yline(K,'--g','LineWidth',1.5,'DisplayName','K')
yline(collapse_threshold,'--r','LineWidth',1.5,'DisplayName','Collapse Threshold')
xline(2026,'-r','LineWidth',2,'DisplayName','Year 2026')
ylabel('Biomass (kg)')
xlabel('Year')
title('Fish Biomass Projection — Mariana Archipelago (CNMI)  (n = 500)')
ylim([0 10000]); xlim([2022 2103])
legend('Location','northeast')
grid on

exportgraphics(gcf,'cnmi_projection.png','Resolution',600)
