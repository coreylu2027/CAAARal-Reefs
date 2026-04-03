clc;
clear;
close all;


set(groot, 'defaultAxesColor',       'w');
set(groot, 'defaultFigureColor',     'w');
set(groot, 'defaultAxesFontName',    'Helvetica');
set(groot, 'defaultAxesFontSize',    11);
set(groot, 'defaultTextFontName',    'Helvetica');
set(groot, 'defaultTextColor',       'k');

CLR_CORAL  = [0.13 0.39 0.68];
CLR_ALGAE  = [0.13 0.60 0.32];
CLR_SPONGE = [0.90 0.62 0.00];
CLR_TREND  = [0.30 0.30 0.30];


function cleanFig(fig)
    for ax = findall(fig, 'Type', 'axes')'
        set(ax, 'Box','off', 'TickDir','out', 'LineWidth',0.9, ...
            'XColor','k', 'YColor','k', 'ZColor','k', ...
            'FontName','Helvetica', 'FontSize',11, 'Color','w');
        if isprop(ax,'Title')
            ax.Title.Color      = 'k';
            ax.Title.FontWeight = 'normal';
            ax.Title.FontSize   = 12;
        end
    end
    set(findall(fig,'Type','text'), 'Color','k', 'FontName','Helvetica');
    for lg = findall(fig,'Type','legend')'
        set(lg, 'Box','off', 'TextColor','k');
    end
    set(fig, 'Color','w');
end


folderPath = '../../data/raw/CORIS_DATA/';
files = dir(fullfile(folderPath, '**', '*Benthic*'));
if isempty(files)
    error('No benthic files found under %s', folderPath);
end
fprintf('%d file(s) found\n', length(files));

raw_tables   = {};
all_col_sets = {};
for f = 1:length(files)
    fpath = fullfile(files(f).folder, files(f).name);
    try
        tbl = readtable(fpath, 'VariableNamingRule','preserve', 'TextType','string');
        if height(tbl) > 0 && ismember('YEAR', tbl.Properties.VariableNames)
            if isnan(str2double(tbl.YEAR(1))), tbl(1,:) = []; end
        end
        raw_tables{end+1}   = tbl;
        all_col_sets{end+1} = tbl.Properties.VariableNames;
    catch ME
        fprintf('Skipped %s: %s\n', files(f).name, ME.message);
    end
end

all_cols = unique([all_col_sets{:}], 'stable');
for f = 1:length(raw_tables)
    for col = setdiff(all_cols, raw_tables{f}.Properties.VariableNames)
        raw_tables{f}.(col{1}) = repmat("", height(raw_tables{f}), 1);
    end
    raw_tables{f} = raw_tables{f}(:, all_cols);
end
data_all  = vertcat(raw_tables{:});
col_names = data_all.Properties.VariableNames;


num_cols = {'YEAR','HARDBOTTOM_P','latitude','longitude','MIN_DEPTH', ...
            'MAX_DEPTH','PROT','WTD_RUG','STATION_NR'};
n_data = height(data_all);
for nc = 1:length(num_cols)
    col = num_cols{nc};
    if ~ismember(col, col_names), continue; end
    v = data_all.(col);
    if isnumeric(v), continue; end
    v = v(:); v(ismissing(v)) = "";
    data_all.(col) = str2double(v);
end


hard_coral_codes = [ ...
    "ACR CERV","ACR PALM","AGA AGAR","AGA FRAG","AGA HUMI","AGA LAMA","AGA SPE.", ...
    "COL NATA","COL SPE.","DEN CYLI","DIC STOK","DIP LABY","DIP SPE.","EUS FAST", ...
    "FAV FRAG","HEL CUCU","ISO SINU","MAD AURE","MAD CARM","MAD DECA","MAD FORM", ...
    "MAD SPE.","MAN AREO","MEA MEAN","MEA SPE.","MIL SPE.","MON CAVE","MUS ANGU", ...
    "MYC ALIC","MYC FERO","MYC LAMA","OCU DIFF","OCU ROBU","ORB ANNU","ORB FAVE", ...
    "ORB FRAN","OTH CORA","POR ASTR","POR DIVA","POR FURC","POR PORI","POR SPE.", ...
    "PSE CLIV","PSE STRI","SCO CUBE","SCO SPE.","SID RADI","SID SIDE","SID SPE.", ...
    "SOL BOUR","SOL SPE.","STE INTE"];

algae_codes = [ ...
    "DIC SPE.","HAL SPE.","LOB SPE.","MAC FLES","MAC CALC","RHO CRUS", ...
    "CYA SPE.","TUR SEDI","TUR FREE","MAG SPE.","PEY SPE.","RAM SPE."];

sponge_codes = ["SPO OTHE","CLI SPE."];


data_all = data_all(~isnan(data_all.YEAR), :);
data_all.COVER_CAT_CD  = strtrim(data_all.COVER_CAT_CD);
if ~isnumeric(data_all.HARDBOTTOM_P)
    data_all.HARDBOTTOM_P = str2double(data_all.HARDBOTTOM_P);
end

col_names = data_all.Properties.VariableNames;
n_rows    = height(data_all);

psu_col = repmat("0", n_rows, 1);
if ismember('PRIMARY_SAMPLE_UNIT', col_names)
    psu_col = data_all.PRIMARY_SAMPLE_UNIT;
    psu_col(ismissing(psu_col)) = "UNKNOWN";
end

stn_col = repmat("0", n_rows, 1);
if ismember('STATION_NR', col_names)
    stn_raw = data_all.STATION_NR;
    stn_raw(isnan(stn_raw)) = -1;
    stn_col = string(stn_raw);
end

year_col_num = data_all.YEAR;
transect_id  = string(year_col_num) + "_" + psu_col + "_" + stn_col;
bad          = contains(transect_id,"NaN") | ismissing(transect_id);
data_all(bad,:)     = [];
transect_id(bad)    = [];
psu_col(bad)        = [];
stn_col(bad)        = [];
year_col_num(bad)   = [];

unique_transects = unique(transect_id);
n_transects      = length(unique_transects);
cat_col          = data_all.COVER_CAT_CD;
hb_col           = data_all.HARDBOTTOM_P;

T = struct( ...
    'year',       zeros(n_transects,1), ...
    'hard_coral', zeros(n_transects,1), ...
    'algae',      zeros(n_transects,1), ...
    'sponge',     zeros(n_transects,1), ...
    'depth_m',    NaN(n_transects,1), ...
    'sub_region', strings(n_transects,1));

for i = 1:n_transects
    rows = find(transect_id == unique_transects(i));
    T.year(i) = year_col_num(rows(1));
    if ismember('MIN_DEPTH',col_names) && ismember('MAX_DEPTH',col_names)
        T.depth_m(i) = (data_all.MIN_DEPTH(rows(1)) + data_all.MAX_DEPTH(rows(1))) / 2;
    end
    if ismember('SUB_REGION_NAME',col_names)
        T.sub_region(i) = data_all.SUB_REGION_NAME(rows(1));
    end
    cats = cat_col(rows);
    hb   = hb_col(rows);
    T.hard_coral(i) = sum(hb(ismember(cats, hard_coral_codes)), 'omitnan');
    T.algae(i)      = sum(hb(ismember(cats, algae_codes)),      'omitnan');
    T.sponge(i)     = sum(hb(ismember(cats, sponge_codes)),     'omitnan');
end

TT = struct2table(T);
TT = TT(~isnan(TT.depth_m) & TT.depth_m > 0, :);

years_all = sort(unique(TT.year));
years_all = years_all(~isnan(years_all));


ann_hc  = arrayfun(@(y) mean(TT.hard_coral(TT.year==y),'omitnan'), years_all);
ann_alg = arrayfun(@(y) mean(TT.algae(TT.year==y),     'omitnan'), years_all);
ann_spo = arrayfun(@(y) mean(TT.sponge(TT.year==y),    'omitnan'), years_all);


figA = figure('Position',[30 30 1500 500], 'Name','3D Benthic Cover — Global');


axA1 = subplot(1,3,1);
scatter3(TT.year, TT.algae, TT.hard_coral, 12, TT.depth_m, ...
    'filled', 'MarkerFaceAlpha', 0.35);
colormap(axA1, flipud(turbo));
cb1 = colorbar; cb1.Label.String = 'Depth (m)'; cb1.Color = 'k';
xlabel('Year'); ylabel('Algae Cover (%)'); zlabel('Hard Coral Cover (%)');
title('Transect Scatter');
view(-35, 25);
grid on;


axA2 = subplot(1,3,2);


year_grid_vec = linspace(min(TT.year), max(TT.year), 60);
alg_grid_vec  = linspace(0, prctile(TT.algae, 97), 60);
[YR_grid, ALG_grid] = meshgrid(year_grid_vec, alg_grid_vec);


valid_pts = ~isnan(TT.hard_coral) & ~isnan(TT.algae);
F_interp  = scatteredInterpolant(TT.year(valid_pts), TT.algae(valid_pts), ...
                                  TT.hard_coral(valid_pts), ...
                                  'natural', 'nearest');
HC_grid = F_interp(YR_grid, ALG_grid);
HC_grid = max(0, HC_grid);

surf(axA2, YR_grid, ALG_grid, HC_grid, 'EdgeAlpha', 0.08, 'FaceAlpha', 0.85);
colormap(axA2, parula);
cb2 = colorbar; cb2.Label.String = 'Hard Coral Cover (%)'; cb2.Color = 'k';
xlabel('Year'); ylabel('Algae Cover (%)'); zlabel('Hard Coral Cover (%)');
title('Interpolated Surface');
view(-40, 30);
grid on;
shading interp;


axA3 = subplot(1,3,3);


ribbon_data = [ann_hc(:), ann_alg(:), ann_spo(:)];
Y_ribbon    = repmat(years_all(:), 1, 3);

ribbon_colors = {CLR_CORAL, CLR_ALGAE, CLR_SPONGE};
ribbon_labels = {'Hard Coral','Algae','Sponge'};
hold(axA3, 'on');
for rb = 1:3

    y_offset = (rb-1) * 0.5;
    z_vals   = ribbon_data(:, rb);
    x_lo     = Y_ribbon(:,rb);
    x_hi     = Y_ribbon(:,rb);
    z_lo     = zeros(size(z_vals));
    z_hi     = z_vals;
    for k = 1:(length(years_all)-1)
        xpatch = [x_lo(k), x_hi(k+1), x_hi(k+1), x_lo(k)];
        ypatch = [y_offset, y_offset, y_offset, y_offset];
        zpatch = [z_lo(k), z_lo(k+1), z_hi(k+1), z_hi(k)];
        fill3(xpatch, ypatch, zpatch, ribbon_colors{rb}, ...
            'FaceAlpha', 0.75, 'EdgeColor', 'none');
    end

    plot3(years_all, repmat(y_offset,size(years_all)), z_hi, ...
        '-', 'Color', ribbon_colors{rb}*0.7, 'LineWidth', 1.5, ...
        'DisplayName', ribbon_labels{rb});
end
xlabel('Year'); zlabel('Cover (%)');
yticks(0:0.5:1); yticklabels({'Coral','','Algae','','Sponge'});
title('Ribbon View');
legend('Location','northeast','FontSize',9);
view(-25, 30);
grid on;

sgtA = sgtitle('Benthic Cover 3D — Coral | Algae | Sponge × Time');
sgtA.FontSize = 13; sgtA.Color = 'k'; sgtA.FontWeight = 'normal';
cleanFig(figA);


sub_regions  = unique(TT.sub_region);
sub_regions  = sub_regions(strtrim(sub_regions) ~= "");
sub_yr_count = arrayfun(@(s) ...
    length(unique(TT.year(TT.sub_region==s))), sub_regions);
plot_regions = sub_regions(sub_yr_count >= 3);
n_reg        = length(plot_regions);

if n_reg > 0
    n_cols_r = min(3, n_reg);
    n_rows_r = ceil(n_reg / n_cols_r);

    figB = figure('Position',[50 50 500*n_cols_r 420*n_rows_r], ...
                  'Name','3D Benthic Cover — Per Region');

    for s = 1:n_reg
        rname  = plot_regions(s);
        mask_r = TT.sub_region == rname;
        sub_d  = TT(mask_r, :);

        axR = subplot(n_rows_r, n_cols_r, s);
        scatter3(sub_d.year, sub_d.algae, sub_d.hard_coral, 20, sub_d.sponge, ...
            'filled', 'MarkerFaceAlpha', 0.55);
        colormap(axR, cool);
        cb_r = colorbar; cb_r.Color = 'k';
        cb_r.Label.String = 'Sponge (%)';


        sub_years  = sort(unique(sub_d.year));
        mean_hc_r  = arrayfun(@(y) mean(sub_d.hard_coral(sub_d.year==y),'omitnan'), sub_years);
        mean_alg_r = arrayfun(@(y) mean(sub_d.algae(sub_d.year==y),     'omitnan'), sub_years);
        hold(axR, 'on');
        plot3(sub_years, mean_alg_r, mean_hc_r, 'k-o', ...
            'LineWidth', 1.8, 'MarkerSize', 5, 'MarkerFaceColor','k', ...
            'DisplayName','Annual mean');
        hold(axR, 'off');

        xlabel('Year'); ylabel('Algae (%)'); zlabel('Coral (%)');
        title(sprintf('%s', strtrim(rname)), 'Interpreter','none');
        view(-35, 25);
        grid on;
    end

    sgtB = sgtitle('Hard Coral vs. Algae vs. Time — By Region');
    sgtB.FontSize = 13; sgtB.Color = 'k'; sgtB.FontWeight = 'normal';
    cleanFig(figB);
end


figC = figure('Position',[80 80 900 650], 'Name','3D Trajectory — Rotating');

axC = axes(figC);
surf(axC, YR_grid, ALG_grid, HC_grid, 'EdgeAlpha',0.05, 'FaceAlpha',0.9);
colormap(axC, parula); shading interp;
cb3 = colorbar; cb3.Label.String = 'Hard Coral Cover (%)'; cb3.Color = 'k';
hold(axC, 'on');


plot3(years_all, ann_alg, ann_hc, 'ko-', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.85 0.12 0.12], 'DisplayName','Annual mean');

xlabel('Year'); ylabel('Algae Cover (%)'); zlabel('Hard Coral Cover (%)');
title('Hard Coral ~ f(Algae, Year) — Interpolated Surface', ...
    'FontWeight','normal', 'Color','k');
legend('Location','northeast','FontSize',10);
grid on;
view(-40, 28);
cleanFig(figC);


fprintf('\n%s\n', repmat('=',1,60));
fprintf('3D Visualisation Summary\n');
fprintf('%s\n', repmat('=',1,60));
fprintf('  Figure A: Global 3D — 3 panels\n');
fprintf('    [A1] Scatter:  Year × Algae × Coral, colored by depth\n');
fprintf('    [A2] Surface:  Interpolated Coral ~ f(Year, Algae)\n');
fprintf('    [A3] Ribbons:  Coral / Algae / Sponge stacked through time\n');
if n_reg > 0
    fprintf('  Figure B: Regional — %d subplots\n', n_reg);
    for s = 1:n_reg
        fprintf('    [%d] %s\n', s, strtrim(plot_regions(s)));
    end
end
fprintf('  Figure C: Rotating surface + annual mean trajectory\n');
fprintf('%s\n', repmat('=',1,60));