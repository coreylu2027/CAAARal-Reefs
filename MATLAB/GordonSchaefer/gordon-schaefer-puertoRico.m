clear; clc; close all;
rng(42);

% K scaled from Florida anchor by reef area ratio
reef_area_PR = 1790;
reef_area_FL = 2900;
K_FL         = 19687072.72/1000;
K            = K_FL * (reef_area_PR / reef_area_FL);
X0_raw = 485.9264; K_raw = 644.5919;
X0     = K * (X0_raw / K_raw);

r       = 0.35;
q       = 0.005;
sigma_X = 0.3177;
E       = 20;

Rmax    = 7.6324;
R0      = 7.6324;
g       = 0.02;
delta   = 0.07;
sigma_R = 0.0847;
alpha   = K / Rmax;

p_shock         = 0.20;
reef_shock_mean = 0.30;
reef_shock_sd   = 0.10;
fish_shock_mean = 0.28;
fish_shock_sd   = 0.10;

p_fish        = 8;
discount_rate = 0.015;
scale_to_global = 284000 / reef_area_PR;
Global_Pop      = 8e9;

T_sim = 80; dt = 0.05; N = round(T_sim/dt);
t = linspace(0,T_sim,N); years = 2023+t; MC = 500;
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

c_blue  = [0.13 0.47 0.71];
c_grey  = [0.50 0.50 0.50];
c_red   = [0.84 0.19 0.15];
c_gold  = [0.80 0.58 0.10];
c_shade = [0.94 0.94 0.94];

fig = figure('Color','w','Position',[80 80 1400 600]);
ax1 = axes;
hold on;
xregion(2026, 2026.5,'FaceColor',c_shade,'FaceAlpha',0.6,'EdgeColor','none')
fill([years,fliplr(years)],[X_upper_deg,fliplr(X_lower_deg)], ...
    c_blue,'FaceAlpha',0.12,'EdgeColor','none')
plot(years, X_mean_nd,  '--','Color',c_grey,'LineWidth',1.4)
plot(years, X_mean_deg, '-', 'Color',c_blue,'LineWidth',2.2)
yline(K, 'Color',c_gold,'LineStyle','--','LineWidth',1.6, ...
    'Label','K  ','LabelHorizontalAlignment','right','FontAngle','italic','FontSize',9)
yline(collapse_threshold, 'Color',c_red,'LineStyle',':','LineWidth',1.4, ...
    'Label','Collapse  ','LabelHorizontalAlignment','right','FontAngle','italic','FontSize',9)
scatter(years(1), X0, 55, c_red,'filled','MarkerEdgeColor','w','LineWidth',0.8)
set(ax1,'Color','w','FontName','Helvetica','FontSize',11,'XColor',[0.2 0.2 0.2],'YColor',[0.2 0.2 0.2], ...
    'LineWidth',0.6,'TickDir','out','Box','off','XTick',2030:10:2100)
ax1.XAxis.MinorTick = 'on'; ax1.XAxis.MinorTickValues = 2025:5:2100;
xlim([2022 2103]); ylim([0 16000])
ylabel('Fish Biomass  (kg)','FontSize',12,'Color',[0.15 0.15 0.15])
xlabel('Year','FontSize',12,'Color',[0.15 0.15 0.15])
title('Fish Biomass Projection — Puerto Rico  (n = 500)','FontSize',13,'FontWeight','bold','Color',[0.1 0.1 0.1])
legend({'Year 2026','90% CI','No degradation','With degradation','Observed (2023)'}, ...
    'FontSize',9.5,'EdgeColor',[0.82 0.82 0.82],'Color','w','Location','northeast')

exportgraphics(fig,'puerto_rico_projection.png','Resolution',600,'BackgroundColor','white')
