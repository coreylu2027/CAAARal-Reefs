clear; clc; close all;
rng(42);

% K scaled from Florida anchor by reef area ratio
reef_area_AS = 520;
reef_area_FL = 2900;
K_FL         = 19687072.72/1000;
K            = K_FL * (reef_area_AS / reef_area_FL);
X0_raw = 4943.6017; K_raw = 11991.3717;
X0     = K * (X0_raw / K_raw);

r       = 0.35;
q       = 0.005;
sigma_X = 0.3857;
E       = 20;

Rmax    = 30.9197;
R0      = 30.9197;
g       = 0.022;
delta   = 0.009;
sigma_R = 0.1390;
alpha   = K / Rmax;

p_shock         = 0.12;
reef_shock_mean = 0.28;
reef_shock_sd   = 0.10;
fish_shock_mean = 0.26;
fish_shock_sd   = 0.10;

p_fish        = 8;
discount_rate = 0.015;
scale_to_global = 284000 / reef_area_AS;

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
        delta_use = delta*(scenario==1);
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
        NPV_loss = sum(max(0,K-X) * p_fish .* exp(-discount_rate*t) * dt);
        if scenario == 1
            NPV_loss_with_deg(m) = NPV_loss; X_all_deg(m,:) = X; R_all(m,:) = R;
            if min(X) < collapse_threshold, collapse_count = collapse_count+1; end
        else
            NPV_loss_no_deg(m) = NPV_loss; X_all_nodeg(m,:) = X;
        end
    end
end

disp(['Collapse probability: ' num2str(collapse_count/MC*100,'%.1f') '%'])
disp(['Degradation premium:  $' num2str((mean(NPV_loss_with_deg)-mean(NPV_loss_no_deg))/1e6,'%.3f') ' M'])

X_p05      = prctile(X_all_deg, 5,  1);
X_p95      = prctile(X_all_deg, 95, 1);
X_mean_deg = mean(X_all_deg,  1);
X_mean_nd  = mean(X_all_nodeg,1);
inst_coll  = mean(X_all_deg < collapse_threshold, 1)*100;
cum_coll   = mean(cummax(X_all_deg < collapse_threshold, 2), 1)*100;

c_blue  = [0.13 0.47 0.71];
c_grey  = [0.45 0.45 0.45];
c_red   = [0.84 0.19 0.15];
c_gold  = [0.72 0.50 0.05];
c_shade = [0.92 0.92 0.92];

fig = figure('Color','w','Position',[60 60 1380 600]);

ax1 = axes('Position',[0.08 0.12 0.86 0.76]);
hold(ax1,'on');
xregion(ax1, 2026, 2026.5, 'FaceColor',c_shade,'FaceAlpha',0.7,'EdgeColor','none')
fill(ax1, [years,fliplr(years)], [X_p95,fliplr(X_p05)], c_blue,'FaceAlpha',0.18,'EdgeColor','none')
fill(ax1, [years,fliplr(years)], [X_mean_nd,fliplr(X_mean_deg)], c_red,'FaceAlpha',0.07,'EdgeColor','none')
plot(ax1, years, X_mean_nd,  '--','Color',c_grey,'LineWidth',1.5)
plot(ax1, years, X_mean_deg, '-', 'Color',c_blue,'LineWidth',2.4)
yline(ax1, K, '--', 'Color',[c_gold 0.75],'LineWidth',1.5)
text(ax1, 2103.2, K, 'K','FontSize',9,'FontAngle','italic','Color',c_gold,'VerticalAlignment','middle')
yline(ax1, collapse_threshold, ':', 'Color',[c_red 0.70],'LineWidth',1.4)
text(ax1, 2103.2, collapse_threshold, 'Collapse','FontSize',8.5,'FontAngle','italic','Color',c_red,'VerticalAlignment','middle')
scatter(ax1, years(1), X0, 60, c_red,'filled','MarkerEdgeColor','w','LineWidth',0.9)
set(ax1,'Color','w','FontName','Helvetica','FontSize',11,'XColor',[0.2 0.2 0.2],'YColor',[0.2 0.2 0.2], ...
    'LineWidth',0.6,'TickDir','out','Box','off','XTick',2030:10:2100)
ax1.XAxis.MinorTick = 'on'; ax1.XAxis.MinorTickValues = 2025:5:2100;
xlim(ax1,[2022 2107]); ylim(ax1,[0 5000])
ylabel(ax1,'Fish Biomass  (kg)','FontSize',12,'Color',[0.15 0.15 0.15])
xlabel(ax1,'Year','FontSize',12,'Color',[0.15 0.15 0.15])
title(ax1,'Fish Biomass Projection — American Samoa  (n = 500)','FontSize',13,'FontWeight','bold','Color',[0.1 0.1 0.1])
legend(ax1, {'Event (2026)','95% CI','Degradation gap','No degradation','With degradation','Observed (2023)'}, ...
    'FontSize',9.5,'EdgeColor',[0.82 0.82 0.82],'Color','w','Location','northeast')

% Inset: collapse probability
ax1_pos = ax1.Position;
ax_ins = axes('Position',[ax1_pos(1)+0.025, ax1_pos(2)+0.30, 0.175, 0.22]);
hold(ax_ins,'on');
fill(ax_ins,[years,fliplr(years)],[cum_coll,fliplr(inst_coll)],c_red,'FaceAlpha',0.12,'EdgeColor','none')
plot(ax_ins,years,inst_coll,'-','Color',[c_red 0.55],'LineWidth',1.1)
plot(ax_ins,years,cum_coll, '-','Color',c_red,'LineWidth',1.7)
set(ax_ins,'Color',[0.99 0.99 0.99],'FontName','Helvetica','FontSize',7.5, ...
    'XColor',[0.35 0.35 0.35],'YColor',[0.35 0.35 0.35],'TickDir','out','Box','on', ...
    'LineWidth',0.5,'XTick',2040:20:2100)
xlim(ax_ins,[2023 2103]); ylim(ax_ins,[0 105])
ylabel(ax_ins,'%','FontSize',7)
title(ax_ins,'Collapse probability','FontSize',7.5,'FontWeight','bold','Color',[0.25 0.25 0.25])
legend(ax_ins,{'Range','Instantaneous','Cumulative'},'FontSize',6.5,'Location','northwest','EdgeColor',[0.80 0.80 0.80])

exportgraphics(fig,'american_samoa_projection.png','Resolution',600,'BackgroundColor','white')
disp('Saved: american_samoa_projection.png')
