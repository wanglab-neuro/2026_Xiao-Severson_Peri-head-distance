% Figure 4:
% *** Requires MATLAB 2024 (?) to run mafdr (multiple comparisons correction)
% Plot pairwise baseline and wall pass firing rates for cutting condition
clear; clc; close all

% Settings
settings.region = 'PrV';
settings.fig_dir = 'E:\PrV_Wall_Recordings\Figures\One_Whisker_cut\Pairwise';
settings.save_plots = true; % false; % 
settings.dist_centers = matfile('E:\Code\Analysis\NeuralDecoder\dist_centers.mat').dist_centers;
settings.dist_edges = (settings.dist_centers(1) - 1) : 2 : (settings.dist_centers(end) + 1);

% Define manipulation condition names
conditions = {'Full','Cut'};

% Load spike stats
settings.wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_OneWhisker';
sst = matfile(fullfile(settings.wv_stats, 'PrV_spike_stats_table_onewhisker.mat')).sst;

% Load psth data
data_folder = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\One_Whisker_cut';
dataTable = matfile(fullfile(data_folder,'tuningTable_PrV_stats.mat')).dataTable;
dataTable_cond = matfile(fullfile(data_folder,'tuningTable_PrV_stats_condition.mat')).dataTable;

%% Load baseline firing rate between control (1=full array) and condition (2=one whisker)

% Set unit inclusion criteria
incl = sst.acceptable & ~dataTable.low_firing;
included_units = find(incl);

% Get mean firing rate and repetitions for baseline and each wall distance
clear fr_base fr_base_distance fr_base_all fr_wall fr_wall_distance fr_wall_all
for c = numel(conditions):-1:1
    disp(['Loading data from ' conditions{c}])
    if strcmp(conditions{c},'Cut') || strcmp(conditions{c},'Laser')
        cond = '_cond';
        dt = dataTable_cond;
    else
        cond = '';
        dt = dataTable;
    end

    % Get median wall firing rate
    fr_wall(:,c) = cellfun(@(x) median(x), dt.(['fr_mean' cond]), 'Uni',1);

    % Get mean firing rate and repetitions for each wall distance
    for u = numel(incl):-1:1
        try
            distances = dt.distanceMM{u};
            [~, ~, included_distances] = histcounts(distances, settings.dist_edges);
            unique_distances = unique(included_distances);
            % tmp = dt.(['fr_base' cond]){included_units(u)};
            % fr_base(u,c) = mean(cellfun(@(x) mean(x), tmp, 'Uni',1));
            for ud = 1:numel(unique_distances)
                d = unique_distances(ud);
                if d > 0
                    fr_base_distance(u,d,c) = mean(cellfun(@(x) mean(x), dt.(['fr_base' cond]){u}(included_distances==d), 'Uni',1));
                    fr_wall_distance(u,d,c) = mean(dt.(['fr_mean' cond]){u}(included_distances==d));
                end
            end
            fr_base(u,c) = mean(fr_base_distance(u,:,c),2);
            fr_base_all{u,c} = cell2mat(dt.(['fr_base' cond]){u}(included_distances)');
            fr_wall_all{u,c} = cell2mat(dt.(['fr_rep' cond]){u}(included_distances)');
        catch e
            disp(['Error on unit ' num2str(u) ': '])
            disp(e)
        end
    end
end

% Get max fr and index of max
num_units = size(fr_wall_distance,1);
[fr_wall_max, fr_wall_max_idx] = deal(NaN(num_units,2));
for u = 1:num_units
    for c = 1:2
        fr_base_max(u,c) = max(squeeze(fr_base_distance(u,:,c)),[],2);
        [fr_wall_max(u,c), fr_wall_max_idx(u,c)] = max(squeeze(fr_wall_distance(u,:,c)),[],2);
    end
end

% Set minimum firing rate effectively to 0 (below 1 Hz)
min_fr = 1; % Hz
wallTuned = incl & dataTable.wallTuned & dataTable.wallResp & fr_wall_max(:,1) > min_fr;
silenced_units = incl & dataTable.activated & wallTuned & fr_wall_max(:,1) > min_fr & fr_wall_max(:,2) < min_fr;
unit_names = dataTable.Properties.RowNames;
silenced_unit_names = unit_names(silenced_units);
disp([num2str(numel(silenced_unit_names)) '/' num2str(sum(wallTuned & dataTable.activated)) ...
    ' tuned & activated units reduced their firing rates below ' num2str(min_fr) ' Hz in manipulation trials'])
disp(shiftdim(silenced_unit_names))

map_incl = incl & strcmp(dataTable.map, 'map');
map_silenced = unit_names(map_incl & silenced_units);

prox_incl = incl & strcmp(dataTable.map, 'proximity');
prox_silenced = unit_names(prox_incl & silenced_units);

% Set inclusion criteria for activated, wall-tuned, acceptable (pre-cut)
mask = incl & dataTable.activated & wallTuned;

disp('Computing bootstrap confidence intervals ...')
ci = [2.5 97.5]; % 95%
nboot = 1e3;
modulation_index = nan(numel(mask),1); % effect size metric
mi_boot = nan(numel(mask),nboot); % effect size metric
for u = 1:numel(incl)
    if incl(u)
        ctrl_rep  = fr_wall_all{u,1};
        manip_rep = fr_wall_all{u,2};
    
        % ensure paired and same length
        if isempty(ctrl_rep) || isempty(manip_rep)
            continue
        end
        L = min(numel(ctrl_rep), numel(manip_rep));
        x = ctrl_rep(1:L);
        y = manip_rep(1:L);
    
        % Compute modulation index
        FR_on = mean(y);
        FR_off = mean(x);
        eps = 1e-9; % epsilon to avoid divide-by-zero errors
        modulation_index(u) = (FR_on - FR_off) ./ (FR_on + FR_off + eps);  % eps to avoid 0/0

        % Compute confidence intervals
        rng('Default')
        rng(42)
        for i = 1:nboot
            idx = datasample(1:numel(x), numel(x), 'Replace',true);
            xx = mean(x(idx));
            yy = mean(y(idx));
            mi_boot(u,i) = (yy - xx) ./ (yy + xx + eps);  % eps to avoid 0/0
        end
    end
end
modulation_index_ci = prctile(mi_boot, ci, 2);

% Use threshold
num_inhib = sum(modulation_index <= -0.5);
num_inc = sum(mask & modulation_index >= 0.05);
num_dec = sum(mask & modulation_index <= -0.05);
num_nochange = sum(mask & abs(modulation_index) < 0.05);

% disp('Increased units:')
% disp(unit_names(mask & modulation_index >= 0.05))
% disp('Decreased units:')
% disp(unit_names(mask & modulation_index <= -0.05))

% Use confidence intervals
sig = all(modulation_index_ci < 0,2) | all(modulation_index_ci > 0,2);
num_inhib = sum(sig & modulation_index <= -0.5);
num_nochange = sum(~sig & mask);
num_inc = sum(sig & mask & modulation_index > 0);
num_dec = sum(sig & mask & modulation_index < 0);

disp(['Strongly inhibited: ' num2str(num_inhib) '/' num2str(sum(incl))])
disp(['Decreased: ' num2str(num_dec) '/' num2str(sum(mask))])
disp(['Increased: ' num2str(num_inc) '/' num2str(sum(mask))])
disp(['No change: ' num2str(num_nochange) '/' num2str(sum(mask))])

disp('Increased units:')
disp(unit_names(sig & mask & modulation_index > 0))
disp('Decreased units:')
disp(unit_names(sig & mask & modulation_index < 0))

% --- Define significance and categories ---

% Significant modulation if CI doesn't cross zero
sig = all(modulation_index_ci < 0, 2) | all(modulation_index_ci > 0, 2);

%% Do stats to find units that significantly increased or decreased firing in manipulation condition
n_units = sum(mask);
units2include = find(mask);
unitnames_tuned = dataTable.Properties.RowNames(mask);
pvals = nan(n_units,1);
signed_diff = nan(n_units,1);   % mean or median diff for direction

for uu = 1:n_units
    u = units2include(uu);
    ctrl_rep  = fr_wall_all{u,1};
    manip_rep = fr_wall_all{u,2};

    % ensure paired and same length
    if isempty(ctrl_rep) || isempty(manip_rep)
        continue
    end
    L = min(numel(ctrl_rep), numel(manip_rep));
    x = ctrl_rep(1:L);
    y = manip_rep(1:L);
    
    % differences
    d = y - x;
    signed_diff(uu) = median(d);

    % Wilcoxon signed-rank test (paired test)
    % H1: median difference != 0
    pvals(uu) = signrank(x, y, 'method','exact');  % exact p-values
end

adj_p = mafdr(pvals, 'BHFDR', true);
alpha = 0.05;

direction = strings(n_units,1);
for u = 1:n_units
    if adj_p(u) < alpha
        if signed_diff(u) > 0
            direction(u) = "increase";
        elseif signed_diff(u) < 0
            direction(u) = "decrease";
        end
    else
        direction(u) = "nochange";
    end
end
num_inc = sum(direction == "increase");
num_dec = sum(direction == "decrease");
num_nochange = sum(direction == "nochange");

disp(['Decreased < 1 Hz: ' num2str(sum(mask & silenced_units)) '/' num2str(sum(mask))])
% disp(['Decreased: ' num2str(num_dec) '/' num2str(sum(mask))])
% disp(['Increased: ' num2str(num_inc) '/' num2str(sum(mask))])
% disp(['No change: ' num2str(num_nc) '/' num2str(sum(mask))])

disp('Increased units:')
disp(unitnames_tuned(direction == "increase"))
disp('Decreased units:')
disp(unitnames_tuned(direction == "decrease"))

% Count activated units tuned before and after
activated_tuned_mask = mask & ~silenced_units & dataTable_cond.activated & dataTable_cond.wallTuned; %  & ~strcmp(dataTable_cond.cut_sig,'significant')
disp([num2str(sum(activated_tuned_mask)) '/' num2str(sum(~silenced_units & mask)) ...
    ' activated units remained tuned in manipulation condition'])
disp(unit_names(activated_tuned_mask))

% Set inclusion criteria for suppressed, wall-tuned, acceptable (pre-cut)
suppressed_mask = incl & dataTable.suppressed & dataTable.wallTuned; % true(size(fr_base,1),1); % 

% Count suppressed units tuned before and after
suppressed_tuned_mask = suppressed_mask & ~silenced_units & dataTable_cond.wallTuned;
disp([num2str(sum(suppressed_tuned_mask)) '/' num2str(sum(~silenced_units & suppressed_mask)) ...
    ' suppressed units remained tuned in manipulation condition'])
disp(unit_names(suppressed_tuned_mask))

%% Run stats on activated population (two-sided Wilcoxon signed rank test for paired observations)

% Set inclusion criteria for activated, wall-tuned, acceptable (pre-cut)
mask = incl & dataTable.activated & dataTable.wallTuned; % true(size(fr_base,1),1); % 

% Expect baseline could change either way (two-sided Wilcoxon signed rank test)
disp(['Full vs. Cut Baseline Period (n=' num2str(sum(incl)) '):'])
base1 = fr_base(incl,1); % e.g., "full" condition
base2 = fr_base(incl,2); % e.g., "cut" condition
[mean_base1, ci_base1, mean_base2, ci_base2] = run_stats_cdf_wenxi(conditions{1}, conditions{2}, base1, base2, 'unequal');
p_base = signrank(base1, base2);

% Expect wall firing rates to change in condition 2 (two-sided Wilcoxon signed rank test)
disp(['Full vs. Cut Wall Period (n=' num2str(sum(incl)) '):'])
wall1 = fr_wall(mask,1); % e.g., "full" condition
wall2 = fr_wall(mask,2); % e.g., "cut" condition
[mean_wall1, ci_wall1, mean_wall2, ci_wall2] = run_stats_cdf_wenxi(conditions{1}, conditions{2}, wall1, wall2, 'unequal');
p_wall = signrank(wall1, wall2);

% Run stats for each unit (how many decreased, increased, no change)
% Wilcoxon signed rank test (paired, non-parametric test)
alpha = 0.05;
units = find(mask);
num_units = numel(units);
alpha_corrected = alpha / (num_units * 2); % correct for multiple comparison (neurons * 2) left and right tail

clear pw pw_left pw_right
for neuron = num_units:-1:1
    % Run signed rank test for each unit to compare firing rates during wall period
    b1 = fr_base_all{neuron,1}(:); % average across distances
    b2 = fr_base_all{neuron,2}(:);
    if numel(b1) > numel(b2)
        b1 = b1(1:numel(b2));
    elseif numel(b2) > numel(b1)
        b2 = b2(1:numel(b1));
    end
    pb(neuron,1) = signrank(b1, b2);
    
    % Run rank sum test for each unit to compare firing rates during wall
    % period in each condition
    w1 = fr_wall_all{neuron,1}(:);
    w2 = fr_wall_all{neuron,2}(:);

    if numel(w1) > numel(w2)
        w1 = w1(1:numel(w2));
    elseif numel(w2) > numel(w1)
        w2 = w2(1:numel(w1));
    end
    pw(neuron,1) = signrank(w1, w2, 'Tail','both');

    % One-sided (right) = reject null hypothesis that w1 - w2 > 0
    % if p < 0.05, w2 > w1
    pw_right(neuron,1) = signrank(w1, w2, 'Tail','right');
    
    % One-sided (left) = reject null hypothesis that w1 - w2 < 0
    % if p < 0.05, w2 > w1
    pw_left(neuron,1) = signrank(w1, w2, 'Tail','left');
end

adj_pb = mafdr(pb, 'BHFDR',true);
adj_pw = mafdr(pw, 'BHFDR',true);

disp([num2str(sum(adj_pb > alpha)) '/' num2str(num_units) ' units show no significant change in baseline firing rates in manipulation condition'])
disp([num2str(sum(adj_pb < alpha)) '/' num2str(num_units) ' units have significantly different baseline firing rates across conditions'])
disp(newline)
disp([num2str(sum(adj_pw > alpha)) '/' num2str(num_units) ' units show no significant change in wall firing rates in manipulation condition'])
disp([num2str(sum(adj_pw < alpha)) '/' num2str(num_units) ' units have significantly different wall firing rates across conditions'])
disp(newline)

decreased_cells = adj_pw < alpha & wall2 < wall1;
increased_cells = adj_pw < alpha & wall2 > wall1;
disp([num2str(sum(decreased_cells)) '/' num2str(num_units) ' tuned units have significantly lower wall firing rates in manipulation condition'])
disp([num2str(sum(increased_cells)) '/' num2str(num_units) ' tuned units have significantly higher wall firing rates in manipulation condition'])

map_mask = map_incl(mask);
map_decreased = map_mask(decreased_cells);
map_increased = map_mask(increased_cells);

%% Panel B: Plot pairwise baseline and wall period firing rates
close all

% Set inclusion criteria for activated, wall-tuned, acceptable (pre-cut)
mask = incl & dataTable.activated & dataTable.wallTuned; % true(size(fr_base,1),1); % 

min_fr = 1e-3; % set ~0 Hz to 1e-3 for log plots
fr_base(fr_base(:,1) < min_fr, 1) = min_fr;
fr_wall(fr_wall(:,1) < min_fr, 1) = min_fr;
fr_wall_max(fr_wall_max(:,1) < min_fr, 1) = min_fr;

fr_base(fr_base(:,2) < min_fr, 2) = min_fr;
fr_wall(fr_wall(:,2) < min_fr, 2) = min_fr;
fr_wall_max(fr_wall_max(:,2) < min_fr, 2) = min_fr;

% Plot pairwise comparison for firing rates during baseline period for
% distance-tuned units (in control condition)
figname = ['Baseline_' conditions{1} '_vs_' conditions{2}];
sig = pb < alpha_corrected;
h3 = plot_pairwise(settings, figname, fr_base(mask,:), conditions, 'Baseline Period', sig);

% Plot pairwise comparison for firing rates during wall period for
% distance-tuned units (in control condition)
figname = ['Wall_' conditions{1} '_vs_' conditions{2}];
sig = pw < alpha_corrected;
h4 = plot_pairwise(settings, figname, fr_wall_max(mask,:), conditions, 'Wall Period', sig);

disp(['Plotting pairwise firing rate for n=' num2str(sum(mask)) ' units'])

%% Panel C: "String" plot comparing preferred distance in Full vs. Cut

distance_centers = 7:2:23;
tuned_units = unit_names(activated_tuned_mask); % 36 units that were tuned in both full and cut whisker conditions

min_fr = 1; % Hz
wallTuned = incl & dataTable.wallTuned & dataTable.wallResp & fr_wall_max(:,1) > min_fr;
silenced_units = incl & dataTable.activated & wallTuned & fr_wall_max(:,1) > min_fr & fr_wall_max(:,2) < min_fr;
mask = incl & dataTable.activated & wallTuned; % true(size(fr_base,1),1); % 
activated_tuned_mask = mask & ~silenced_units & dataTable_cond.activated & dataTable_cond.wallTuned; %  & ~strcmp(dataTable_cond.cut_sig,'significant')
disp([num2str(sum(activated_tuned_mask)) '/' num2str(sum(~silenced_units & mask)) ...
    ' activated units remained tuned in manipulation condition'])
disp(tuned_units)

% Get labels (map, ambiguous, proximity) from control and manipulation condition
label_pre = dataTable.map(activated_tuned_mask);
label_post = dataTable_cond.map(activated_tuned_mask);

labels = {'map','ambiguous','proximity'};
symbols = {'M','A','P'};
colors = [1.0 1.0 0  % yellow
          1.0 0.5 0; % orange
          0.6 0.0 0];% dark red

% Get the number of units
data = fr_wall_max_idx(activated_tuned_mask,:);
num_units = size(data, 1);

% Convert preferred distance index to actual mm values
d1_all = distance_centers(data(:,1)); % Pre-manipulation (control)
d2_all = distance_centers(data(:,2)); % Post-manipulation (manipulation)

% Get manipulation effects (same distance, closer, farther)
effect_labels = {'Same','Closer','Farther'};
for i = numel(d1_all):-1:1
    if d1_all(i) == d2_all(i)
        effects(i) = 1; % same across conditions
    elseif d1_all(i) > d2_all(i)
        effects(i) = 2; % closer in manipulation condition
    elseif d1_all(i) < d2_all(i)
        effects(i) = 3; % farther in manipulation
    end
end

% Sort units by pre-manipulation preferred distance, then post-manipulation
[~, sort_idx] = sortrows([shiftdim(effects), d1_all(:), d2_all(:)]);

% Apply sorting to all arrays
d1_sorted = d1_all(sort_idx);
d2_sorted = d2_all(sort_idx);
tuned_units_sorted = tuned_units(sort_idx);
label_pre = label_pre(sort_idx);
label_post = label_post(sort_idx);
effects = effects(sort_idx);

mark_big = 12;
mark_sm = 8;

% Create the plot
figure(1); clf
set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',[2 2 8 18])
hold on;
yticks(1:num_units);
% yticklabels(strrep(tuned_units_sorted,'_',' '));
xlabel('Preferred Distance (mm)');
ylabel('Unit (sorted by pre → post distance)');

% Set x and y limits
xlim([min(distance_centers)-2, max(distance_centers)+1]);
ylim([0.25, num_units + 0.75]);

% Plot each unit's before/after preferred distance
for i = 1:num_units
    x1 = d1_sorted(i); % Pre-manipulation (control)
    x2 = d2_sorted(i); % Post-manipulation (manipulation)
    y = i;

    c1 = strcmp(labels, label_pre{i});
    c2 = strcmp(labels, label_post{i});
    if x2 < x1
        offset = -1;
    else
        offset = 1;
    end

    % Line connecting the two points
    plot([x1, x2], [y, y], 'k-', 'LineWidth',0.5);

    % --- Cut condition (post-manipulation) as a colored square ---
    % Plot first so it sits "behind" the full-array circle when x1 == x2
    plot(x2, y, 's', ...                         % square marker
        'MarkerEdgeColor','k', ...
        'MarkerFaceColor',colors(c2,:), ...      % color by post-manip label
        'MarkerSize',mark_big, ...
        'LineWidth',0.5);                         % slightly larger than full-array circle

    % --- Full-array condition (pre-manipulation) as a colored circle ---
    % Plot second so it appears on top
    plot(x1, y, 'o', ...
        'MarkerEdgeColor','k', ...
        'MarkerFaceColor',colors(c1,:), ...
        'MarkerSize',mark_sm, ...
        'LineWidth',0.5);                         % slightly smaller than square
end

set(gca, 'YDir', 'reverse'); % Optional: unit 1 at top
set(gca, 'XTick',distance_centers)
set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)

% Create dummy plots for the legend
h_full = plot(nan, nan, 'ko', 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'MarkerSize',mark_big); % Full Array (closed)
h_cut = plot(nan, nan, 'ks', 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'MarkerSize',mark_big);
h_prox = plot(nan, nan, 'o', 'MarkerFaceColor',colors(strcmp(labels,'proximity'),:), 'MarkerSize',mark_big); 
h_amb = plot(nan, nan, 'o', 'MarkerFaceColor',colors(strcmp(labels,'ambiguous'),:), 'MarkerSize',mark_big); 
h_map = plot(nan, nan, 'o', 'MarkerFaceColor',colors(strcmp(labels,'map'),:), 'MarkerSize',mark_big); 

% Add the legend
leg = legend([h_full, h_cut, h_prox, h_amb, h_map], ...
    {'Full Array', 'Cut', 'Proximity', 'Ambiguous', 'Map'}, ...
    'Location','best');
set(leg, 'Box','off', 'Location','northeast', 'FontSize',8)

if settings.save_plots
    filename = fullfile(settings.fig_dir, ['pairwise_onewhisker_cut_distance_stacked']);
    saveas(gcf,filename,'pdf');
    exportgraphics(gcf,[filename '.png'], 'Resolution',600);
    pause(0.01)
end

% Count number of map cells that remained tuned after cut
map_before = strcmp(dataTable.map, 'map');
map_before_tuned_after = map_before & dataTable_cond.wallTuned & dataTable_cond.activated;
disp([num2str(sum(map_before_tuned_after)) ' out of ' num2str(sum(map_before)) ...
    ' map cells remained tuned after cutting'])

% Count number of proximity cells that remained tuned after cut
prox_before = strcmp(dataTable.map, 'proximity');
prox_before_tuned_after = prox_before & dataTable_cond.wallTuned & dataTable_cond.activated;
disp([num2str(sum(prox_before_tuned_after)) ' out of ' num2str(sum(prox_before)) ...
    ' proximity cells remained tuned after cutting'])

% Count number of tuned proximity cells that increased/decreased activity
prox_before_tuned_after_increased = prox_before_tuned_after(mask) & increased_cells;
prox_before_tuned_after_decreased = prox_before_tuned_after(mask) & decreased_cells;
disp([num2str(sum(prox_before_tuned_after_increased)) ' increased activity after cutting'])
disp([num2str(sum(prox_before_tuned_after_decreased)) ' decreased activity after cutting'])

%% Functions

function h = plot_pairwise(settings, figname, data, conditions, name, sig)
    alpha = 0.25;

    figure(randi(1000)); clf
    hold on
    w = 4; h = 6;
    figpos = [5 5 w h];
    MPlot.Paperize(gcf)
    set(gcf, 'Position',figpos)
    yy = data;
    xx = cat(2, zeros(size(yy,1),1), ones(size(yy,1),1) );
    plot(xx', yy', 'ko', 'MarkerSize',4) % plot points
    plot(xx', yy', 'Color',[0 0 0 alpha]) % plot lines
    text(xx(sig,2)+0.1, yy(sig,2), '*', 'FontSize',8, 'Color','k')
    MPlot.FormatAxes(gca)
    xlim([-0.5 1.5])
    ylim([3e-2 3e2])
    set(gca, 'YScale','log') % set y-axis to log scale
    ax = gca;
    ylabel('Firing Rate (Sp/s)')
    yticks = logspace(-1,2,4);
    set(gca, 'YTick',yticks)
    ytlabels = ax.YTickLabels;
    set(gca, 'YTickLabels',ytlabels)
    set(gca, 'XTick',[0 1], 'XTickLabels',conditions)
    title(name)

    if settings.save_plots
        filename = fullfile(settings.fig_dir, ['pairwise_' figname]);
        saveas(gcf,filename,'pdf');
        exportgraphics(gcf,[filename '.png'], 'Resolution',600);
        pause(0.01)
    end
end


function h = plot_pairwise_pref_distance(settings, figname, data, conditions, name, distance_centers)
    alpha = 0.25;

    figure(randi(1000)); clf
    hold on
    w = 4; h = 6;
    figpos = [5 5 w h];
    MPlot.Paperize(gcf)
    set(gcf, 'Position',figpos)
    yy = distance_centers(data); % convert index to mm
    xx = cat(2, zeros(size(yy,1),1), ones(size(yy,1),1) );
    plot(xx', yy', 'ko', 'MarkerSize',4) % plot points
    plot(xx', yy', 'Color',[0 0 0 alpha]) % plot lines
    MPlot.FormatAxes(gca)
    xlim([-0.5 1.5])
    ylim([min(distance_centers)-1 max(distance_centers)+1])
    ax = gca;
    ylabel('Pref. Wall Distance (mm)')
    set(gca, 'YTick',distance_centers)
    set(gca, 'XTick',[0 1], 'XTickLabels',conditions)
    title(name)

    if settings.save_plots
        filename = fullfile(settings.fig_dir, ['pairwise_pref_dist_' figname]);
        saveas(gcf,filename,'pdf');
        exportgraphics(gcf,[filename '.png'], 'Resolution',600);
        pause(0.01)
    end
end


function make_wall_cdf_wenxi(settings, figname, data1, data2, color1, color2, name1, name2)
    
    % Clamp data1 and data2 to range
    min_fr = 1e-3;
    data1(data1 < min_fr) = min_fr;
    data2(data2 < min_fr) = min_fr;

    % Make figure for cdf comparing two datasets
    figure(randi(1000)); clf;
    set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',[2 2 5 5])
    hold on
    h1 = cdfplot(data1);
    h2 = cdfplot(data2);
    set(h1, 'Color',color1, 'LineWidth',1);
    set(h2, 'Color',color2, 'LineWidth',1);
    MPlot.FormatAxes(gca)
    set(gca, 'XScale','log')
    axis square
    grid off
    xlim([min_fr 100])
    set(gca, 'XTick',[min_fr .01 .1 1 10 100])
    ax = gca;
    xtlabels = ax.XTickLabels;
    xtlabels{1} = '0';
    ylabel('Cumulative Fraction Units')
    xlabel('Firing Rate (Sp/s)')
    title(strrep(['CDF ' figname],'_',' '))
    legend({name1,name2}, 'Location','northwest', 'Box','off')

    if settings.save_plots
        filename = fullfile(settings.fig_dir, ['cdf_' figname]);
        saveas(gcf,filename,'pdf');
        exportgraphics(gcf,[filename '.png'], 'Resolution',600);
        pause(0.01)
    end
end


function [mean1, ci1, mean2, ci2, p, stats] = run_stats_cdf_wenxi(g1, g2, data1, data2, tail)
    if strcmp(tail,'smaller')
        tail = 'right';
    elseif strcmp(tail, 'larger')
        tail = 'left';
    else
        tail = 'both';
    end

    % Fix random seed
    rng('default')
    rng(42);

    % Compute mean +/- CI for Group 1 data
    mean1 = mean(data1);
    ci1 = diff(bootci(1e4, @(x) mean(x), data1)) / 2;
    disp(['Mean for ' g1 ' : ' num2str(mean1) ' +/- ' num2str(ci1)])
    
    % Fix random seed
    rng('default')
    rng(42);

    % Compute mean +/- CI for Group 2 data
    mean2 = mean(data2);
    ci2 = diff(bootci(1e4, @(x) mean(x), data2)) / 2;
    disp(['Mean for ' g2 ' : ' num2str(mean2) ' +/- ' num2str(ci2)])
    

    % Run Mann-Whitney U Test for significance
    [p,~,stats] = signrank(data1, data2, 'Tail',tail);
    disp(['P-value from Wilcoxon Signed Rank Test: ' num2str(p)])
    disp(newline)
end
