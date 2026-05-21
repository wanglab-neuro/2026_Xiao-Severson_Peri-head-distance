%% Fig S4, panel I
% Analyze probe map information

fig_dir = 'C:\Users\kssev\Desktop\Manuscript\revision\Figures\Figure S4';

% Load sst with probe map info
waveform_stats = 'E:\PrV_Wall_Recordings\WaveformStats_naive';
sst_file = 'PrV_spike_stats_table_naive_probe.mat';

% Load psth data table with units classified into map, proximity, etc.
dataTable = matfile('E:\PrV_Wall_Recordings\Analysis\IntermediateData\tuningTable_PrV_stats.mat').dataTable;

% Load sst
sst = matfile(fullfile(waveform_stats,sst_file)).sst;

% Plot histogram
population_names = {'Proximity', 'Map', 'Suppressed', 'Untuned'}; % 'Ambiguous', unique(dataTable.map)

% Get colormap
cmap = flipud(turbo(256));
colors = zeros(numel(population_names),3);
colors(strcmp(population_names,'Proximity'),:) = cmap(1,:); % proximity = red
colors(strcmp(population_names,'Untuned'),:) = [0.5 0.5 0.5]; % untuned = dark gray
colors(strcmp(population_names,'Map'),:) = [1 1 0]; % map = white (yellow)
colors(strcmp(population_names,'Suppressed'),:) = [0.4 0.4 1]; % suppressed = blue

edges = 0:100:800; % depth bins
binCenters = mean([edges(1:end-1); edges(2:end)]);

% mask = ~dataTable.low_firing & ~contains(sst.probe,'poly3'); % forgot "acceptable" mask
poly3_mask = ~dataTable.low_firing & sst.acceptable & dataTable.wallResp & strcmpi(sst.probe, 'NN poly3');
poly2_mask = ~dataTable.low_firing & sst.acceptable & dataTable.wallResp & strcmpi(sst.probe, 'NN poly2');
depths = sst.depth(poly2_mask);
map_class = dataTable.map(poly2_mask);

% disp(['Total poly2: ' num2str(sum(poly2_mask))])
% disp(['Total poly3: ' num2str(sum(poly3_mask))])

depths = sst.depth(poly2_mask);
map_class = dataTable.map(poly2_mask);

figure(1); clf
tl = tiledlayout(1, numel(population_names), 'Padding','tight', 'TileSpacing','compact');
for pop = 1:numel(population_names)
    nexttile
    pop_name = population_names{pop};
    set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',[5 5 10 4]);
    x = binCenters;
    y = histcounts(depths(strcmpi(map_class, pop_name)), edges);
    total = sum(y);
    y = y / total;
    barh(x, y, 1, 'FaceColor',colors(pop,:), 'EdgeColor','k')
    title([pop_name ' (N=' num2str(total) ')'], 'FontWeight','normal', 'FontSize',6)
    xlabel('Frac. units')
    set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)
    if pop==1
        ylabel('Depth from tip (um)')
        ylab = get(gca, 'YTick');
        ylab = arrayfun(@(x) num2str(-x), ylab, 'Uni',0);
        set(gca, 'YTickLabels',ylab)
    else
        set(gca, 'YTickLabels',{})
    end
    xlim([0 0.275])
    ylim([-25 825])
    pause(0.001)
end

% Save figure
saveas(gcf, fullfile(fig_dir, 'S4_I_probe_map_poly2_all'), 'pdf')
saveas(gcf, fullfile(fig_dir, 'S4_I_probe_map_poly2_all'), 'png')

% mask = ~dataTable.low_firing & dataTable.wallResp & ~contains(sst.probe,'poly3');
% depths = sst.depth(mask);
% map_class = dataTable.map(mask);
% 
% figure(2); clf
% tl = tiledlayout(1, numel(population_names), 'Padding','tight', 'TileSpacing','compact');
% for pop = 1:numel(population_names)
%     nexttile
%     pop_name = population_names{pop};
%     set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',[5 5 10 4]);
%     x = binCenters;
%     y = histcounts(depths(strcmpi(map_class, pop_name)), edges);
%     total = histcounts(depths(ismember(map_class, lower(population_names))), edges);
% %     total = sum(y);
%     y = y / total;
%     barh(x, y, 1, 'FaceColor',colors(pop,:), 'EdgeColor','k')
%     title([pop_name ' (N=' num2str(total(pop)) ')'], 'FontWeight','normal', 'FontSize',6)
%     xlabel('Frac. units')
%     set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)
%     if pop==1
%         ylabel('Depth from tip (um)')
%         ylab = get(gca, 'YTick');
%         ylab = arrayfun(@(x) num2str(-x), ylab, 'Uni',0);
%         set(gca, 'YTickLabels',ylab)
%     else
%         set(gca, 'YTickLabels',{})
%     end
%     xlim([0 0.6])
%     ylim([-25 825])
%     pause(0.001)
% end
% 
% % Save figure
% saveas(gcf, fullfile(fig_dir, 'S4_I_probe_map_poly2_resp'), 'pdf')
% saveas(gcf, fullfile(fig_dir, 'S4_I_probe_map_poly2_resp'), 'png')

response_types = {'Wall-responsive', 'Unresponsive'};

poly3_mask = ~dataTable.low_firing & sst.acceptable & strcmpi(sst.probe, 'NN poly3');
poly2_mask = ~dataTable.low_firing & sst.acceptable & strcmpi(sst.probe, 'NN poly2');

depths = sst.depth(poly2_mask);
map_class = dataTable.map(poly2_mask);

resp_class = repmat({'Unresponsive'},numel(depths),1);
resp_class(dataTable.wallResp(poly2_mask)) = {'Wall-responsive'};

resp_colors = [1 1 1; 0.25 0.25 0.25]; % responsive

xl = [0 0.2];
yl = [-25 825];

figure(4); clf
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
for pop = 1:2
    nexttile
    pop_name = response_types{pop};
    set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',[5 5 4 4]);
    x = binCenters;
    y_resp = histcounts(depths(strcmpi(resp_class, 'Wall-responsive')), edges);
    y_unresp = histcounts(depths(strcmpi(resp_class, 'Unresponsive')), edges);
    total = y_resp + y_unresp; % sum over each depth reange
    if strcmp(response_types{pop}, 'Wall-responsive')
        N = sum(y_resp);
        y = y_resp / N;
    else
        N = sum(y_unresp);
        y = y_unresp / N;
    end
    barh(x, y, 1, 'FaceColor',resp_colors(pop,:), 'EdgeColor','k')
    title(['N=' num2str(N)], 'FontWeight','normal', 'FontSize',6)
    xlabel('Frac. units')
    set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)
    if pop==1
        ylabel('Depth from tip (um)')
        ylab = get(gca, 'YTick');
        set(gca, 'YTickLabels',ylab)
    else
        set(gca, 'YTickLabels',{})
    end
    xlim(xl)
    ylim(yl)
    pause(0.001)
end

% Save figure
filename = 'S4_I_probe_map_poly2_fraction_depth_resp';
saveas(gcf, fullfile(fig_dir, filename), 'pdf')
saveas(gcf, fullfile(fig_dir, filename), 'png')

