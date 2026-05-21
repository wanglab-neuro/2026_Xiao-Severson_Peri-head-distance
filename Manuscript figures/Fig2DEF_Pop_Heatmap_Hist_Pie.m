% Run stats to test significance of wall tuning
% Also compute tuning modulation index (wall - baseline) / (wall + baseline)
% and assign 'activated' vs. 'suppressed' label
clear; clc
settings.region = 'PrV';
settings.output_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData';
settings.fig_dir = 'E:\PrV_Wall_Recordings\Figures';

%% Naive -- 1 way stats - wall distance for naive experiment (naive + before cut)
settings.dataset = 'naive';
settings.data_dir = 'E:\PrV_Wall_Recordings';
settings.wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_naive';

% Also computes "Activation" and "Suppression"
clearvars -except settings

% Load tuning curve data combined from naive and whisker cutting
% experiments (1358 total units)
settings.data_folder = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData';
dataTable = matfile(fullfile(settings.data_folder, 'tuningTable_PrV.mat')).dataTable;
dataTable_cut = matfile(fullfile('E:\PrV_Wall_Recordings\Analysis\IntermediateData\Whisker_Cutting_naive\tuningTable_PrV.mat')).dataTable;
dataTable_one = matfile(fullfile('E:\PrV_Wall_Recordings\Analysis\IntermediateData\One_Whisker_naive\tuningTable_PrV.mat')).dataTable;
dataTable = cat(1, dataTable, dataTable_cut(:,1:size(dataTable,2)), dataTable_one(:,1:size(dataTable,2))); % only take columns without condition labels
clear dataTable_cut dataTable_one

% Load spike stats table (1358 units)
sst = matfile(fullfile(settings.wv_stats, 'PrV_spike_stats_table_naive.mat')).sst;

%% Prerequisite in Figure 2: Run stats for naive full array condition

rng('Default')
rng(42)

% Define units to be included
threshold_stats = 1; % 1.0; % mean baseline or wall fr, in Hz
dataTable.low_firing = cellfun(@(b,w) max(b) < threshold_stats & mean(w) < threshold_stats, ...
    dataTable.fr_mean, dataTable.fr_base_mean, 'Uni',1);

included_units = sst.acceptable & ~dataTable.low_firing;
num_units = sum(included_units);

% Define distance conditions
if strcmp(settings.dataset, 'tungsten')
    num_distances = []; % empirical based on each session
else
    num_distances = [];
end

% Define p-value with Bonferroni correction
alpha = 0.05;
alpha_corrected = alpha; %  / (nchoosek(10,2)); %  / (num_distances * 2); % Bonferroni correction for multiple independent comparisons (units)

% Pre-allocate
dataTable.fr_mean_part1 = cell(size(dataTable,1),1);
dataTable.fr_mean_part2 = cell(size(dataTable,1),1);
dataTable.wallResp = false(size(dataTable,1),1);
dataTable.wallTuned = false(size(dataTable,1),1);
dataTable.stats_p_tuned = cell(size(dataTable,1),1);
dataTable.stats_p_map = cell(size(dataTable,1),1);
dataTable.stats_err = false(size(dataTable,1),1);
dataTable.activated = false(size(dataTable,1),1);
dataTable.suppressed = false(size(dataTable,1),1);
dataTable.map = cell(size(dataTable,1),1);
dataTable.map_ambiguous = false(size(dataTable,1),1);
dataTable.Properties.RowNames = cellfun(@(x,u) strjoin({x,num2str(u)},'_'), dataTable.recID, num2cell(dataTable.unit), 'Uni',0);
dataTable.pref_dist = NaN(size(dataTable,1),1);

mod_index_thresh = 0;
show_suppressed = 0;

error_trials = [];
cnt = 1;

for neuron = size(dataTable,1):-1:1
    if included_units(neuron)
        try
            disp(neuron)

            % Get data from baseline and wall pass periods
            baseline_data = dataTable.fr_base{neuron};
            wallpass_data = dataTable.fr_rep{neuron};
            distance_data = unique(dataTable.distanceMM{neuron});

            % Handle missing data 
            if strcmp(dataTable.recID(neuron), '20221213_sp_session_2')
                % Remove last 3 loops to balance
                baseline_data = cellfun(@(x) x(1:17), baseline_data, 'Uni',0);
                wallpass_data = cellfun(@(x) x(1:17), wallpass_data, 'Uni',0);
            elseif strcmp(dataTable.recID(neuron), '20230131_sp_session_2')
                % Remove last loop to balance
                baseline_data = cellfun(@(x) x(1:19), baseline_data, 'Uni',0);
                wallpass_data = cellfun(@(x) x(1:19), wallpass_data, 'Uni',0);
            end
            
            % Compute activation/suppression scores, modulation index
            base_mean = cellfun(@(x) mean(x,'omitnan'), dataTable.fr_base{neuron}, 'Uni',1);
            wall_mean = cellfun(@(x) mean(x,'omitnan'), dataTable.fr_rep{neuron}, 'Uni',1); % dataTable.fr_mean{neuron};
            wall_mean_z = dataTable.zfr_mean{neuron};
            wall_min = min(dataTable.fr_mean{neuron}, [], 'omitnan');
            wall_max = max(dataTable.fr_mean{neuron}, [], 'omitnan');

            % If any remaining non-paired distances, remove them here
            num_passes_per_dist = cellfun(@(x) numel(x), wallpass_data, 'Uni',1);
            % d2remove = num_passes_per_dist ~= max(num_passes_per_dist);
    
            % Get wall distance groups
            reps = numel(baseline_data{1,1}); % count number of loops (wall distance replications);
            num_wd = numel(distance_data);
            wd = repelem(1:num_wd, reps)';

            % Partition data into 2 groups for cross-validation
            idx = true(1,reps);
            idx(datasample(1:reps, floor(reps/2), 'Replace',false)) = false;
            dataTable.fr_mean_part1{neuron,:} = cell2mat(cellfun(@(x) mean(x(idx),1,'omitnan'), wallpass_data, 'Uni',0));
            dataTable.fr_mean_part2{neuron,:} = cell2mat(cellfun(@(x) mean(x(~idx),1,'omitnan'), wallpass_data, 'Uni',0));
            
            % Reshape datatable for 1-way stats to test for responsiveness 
            % baseline different from wall firing rates
            bd = baseline_data;
            wpd = wallpass_data;

            % y = cell2mat(auc_data);
            if numel(wpd) ~= numel(wd)
                error_trials(cnt) = neuron;
                cnt = cnt + 1;
            end
            
            % ----------------------------------------------
            %             RESPONSIVENESS TEST 
            %           ACTIVATION / SUPPRESSION
            % ----------------------------------------------

            % Manually perform one-tailed stats on each distances to
            % compare firing rates versus baseline
            if isempty(num_distances)
                nd = numel(unique(wd));
            else
                nd = num_distances;
            end

            N = 0:nd;
            combos = cat(2, zeros(nd,1), (1:nd)');
            [p0_act, p0_sup] = deal(zeros(size(combos,1),1));
            d1 = cell(size(combos,1),1);
            d2 = cell(size(combos,1),1);
            for d = 1:nd
                d2{d} = wpd{d}; % wall
                d1{d} = bd{d}; % baseline
                
%                 [~, p0_act(d)] = ttest(d2{d}, d1{d}, 'Tail','right'); % one-tailed paired t-test
%                 [~, p0_sup(d)] = ttest(d2{d}, d1{d}, 'Tail','left'); % one-tailed paired t-test

                p0_act(d) = signrank(d2{d}, d1{d}, 'Tail','right'); % one-tailed paired non-parametric test
                p0_sup(d) = signrank(d2{d}, d1{d}, 'Tail','left'); % one-tailed paired non-parametric test
            end
            
            c = zeros(size(combos,1), 6);
            c(:,1) = combos(:,1);
            c(:,2) = combos(:,2);
            c(:,3) = cellfun(@(x) mean(x), d1, 'Uni',1);
            c(:,4) = cellfun(@(x) mean(x), d2, 'Uni',1);
            c(:,5) = p0_act; % raw activation p-values
            c(:,6) = p0_sup; % raw suppression p-values

            % Get p-values for testing activation vs. suppression (wall vs. baseline)
            p0_act = mafdr(c(c(:,1)==0, 5), 'BHFDR',true); % wall > baseline
            p0_sup = mafdr(c(c(:,1)==0, 6), 'BHFDR',true); % wall < baseline
            q_resp = cat(1, p0_act, p0_sup); % concat corrected p-values

            dataTable.ttest_p_activated{neuron} = p0_act;
            dataTable.ttest_p_suppressed{neuron} = p0_sup;

            % Activated: any significant wall > baseline
            if any(dataTable.ttest_p_activated{neuron} < alpha)
                dataTable.activated(neuron) = true;
            end
            
            % Suppressed: any significant wall < baseline
            if any(dataTable.ttest_p_suppressed{neuron} < alpha)
                dataTable.suppressed(neuron) = true;
            end
            
            if any(q_resp < alpha)
                dataTable.wallResp(neuron) = true;
            end

            % ----------------------------------------------
            %                 TUNING TEST 
            % ----------------------------------------------
            
            % Get index of preferred distance
            [~, pref_dist] = max(wall_mean_z);
            [~, pref_dist_sup] = min(wall_mean_z);

            % Separately, run one-way stats to determine distance tuning (any
            % distance condition pair significantly different)
            y = cell2mat(wpd);
            distance = categorical(wd);
            [~, ~, stats] = anova1(y, distance, 'off'); % parametric paired t-test
            [c, means] = multcompare(stats, 'Display','off');
            
            % Manually compute stats with independent t-tests
            N = 1:nd;
            combos = nchoosek(1:nd, 2);
            p1 = zeros(size(combos,1),1);
            d1 = cell(size(combos,1),1);
            d2 = cell(size(combos,1),1);
            for i = 1:numel(combos(:,1))
                d2{i} = wpd{combos(i,2)}; % wall
                if combos(i,1)==0
                    d1{i} = bd{combos(i,2)}; % baseline
                else
                    d1{i} = wpd{combos(i,1)}; % wall
                end
                
                % [~, p1(i)] = ttest(d2{i}, d1{i}); % two-tailed paired t-test
                p1(i) = signrank(d2{i}, d1{i}); % two-tailed paired non-parametric test
            end
            
            c = zeros(size(combos,1), 5);
            c(:,1) = combos(:,1);
            c(:,2) = combos(:,2);
            c(:,3) = cellfun(@(x) mean(x), d1, 'Uni',1);
            c(:,4) = cellfun(@(x) mean(x), d2, 'Uni',1);
            c(:,5) = p1(:,1); % raw two-tailed p-values

            % Get p-values for testing tuning (any difference between wall distance pairs)
            p_tuned = c(:, 5);
            p_tuned = mafdr(p_tuned(:), 'BHFDR',true);

            % ----------------------------------------------
            %          AMBIGUOUS MAP/PROXIMITY TEST 
            % ----------------------------------------------

            % Get p-values for testing map vs. ambiguous
            N = 2:nd;
            combos = cat(2, ones(numel(N),1), N');
            p2 = zeros(size(combos,1),1);
            d1 = cell(size(combos,1),1);
            d2 = cell(size(combos,1),1);
            for i = 1:numel(combos(:,1))
                d2{i} = wpd{combos(i,2)}; % wall
                if combos(i,1)==0
                    d1{i} = bd{combos(i,2)}; % baseline
                else
                    d1{i} = wpd{combos(i,1)}; % wall
                end
                
%                 [~, p2(i)] = ttest(d2{i}, d1{i}, 'Tail','right'); % one-tailed paired t-test
                p2(i) = signrank(d2{i}, d1{i}, 'Tail','right'); % one-tailed paired t-test
            end
            
            c = zeros(size(combos,1), 5);
            c(:,1) = combos(:,1);
            c(:,2) = combos(:,2);
            c(:,3) = cellfun(@(x) mean(x), d1, 'Uni',1);
            c(:,4) = cellfun(@(x) mean(x), d2, 'Uni',1);
            c(:,5) = p2(:,1); % raw one-tailed p-values

            if dataTable.activated(neuron) && any(p_tuned < alpha)
                p_map = mafdr(c(:, 5), 'BHFDR',true); % corrected
            else
                p_map = [];
            end

            % Store p-values
            dataTable.stats_p_tuned{neuron} = p_tuned;
            dataTable.stats_p_map{neuron} = min(p_map);
        
            % Test if significantly wall-tuned (any distance different from any other distance)
            if any(p_tuned < alpha) % stats multcomp already corrects for multiple comparisons across distances
                % Denote neurons with significant stats p-value as "wall-tuned"
                dataTable.wallTuned(neuron) = true;
    
                % Get error bounds around baseline firing rate estimation (mean +/- SEM)
                % Could also try 95% or 99% confidence intervals
                N = numel(dataTable.fr_rep{neuron}{1}); % num repetitions
                
                % Plot sanity test
                if show_suppressed
                    if dataTable.suppressed(neuron)
                        figure(neuron); clf
                        hold on
                        y = dataTable.fr_mean{neuron};
                        err_lo = dataTable.fr_mean{neuron} - dataTable.fr_ci{neuron}(:,1);
                        err_hi = dataTable.fr_ci{neuron}(:,2) - dataTable.fr_mean{neuron};
                        errorbar(1:numel(y), y, err_lo, err_hi, 'Color','k')
                        plot(dataTable.fr_mean{neuron}, 'k', 'LineWidth',2)
                        base_mean = cellfun(@(x) mean(x,'omitnan'), dataTable.fr_base{neuron}, 'Uni',1);
                        y = mean(base_mean);
                        err_lo = repmat(y - dataTable.fr_base_ci{neuron}(1), numel(dataTable.fr_mean{neuron}), 1);
                        err_hi = repmat(dataTable.fr_base_ci{neuron}(2) - y, numel(dataTable.fr_mean{neuron}), 1);
                        y = repmat(y,numel(dataTable.fr_mean{neuron}),1);
                        errorbar(1:numel(y), y, err_lo, err_hi, 'Color','m', 'LineWidth',2)
                        pause(0.001)
                    end
                end
            end

            % Determine map vs. monotonic cell (for activated units)
            if dataTable.wallTuned(neuron)
                if dataTable.activated(neuron) 
                    dataTable.pref_dist(neuron) = pref_dist;

                    if pref_dist==1
                        dataTable.map{neuron} = 'proximity';
                    elseif pref_dist > 1
                        dataTable.map{neuron} = 'map';
                    end
    
                    % Determine if map/monotonic is ambiguous (closest distance is not significantly different from preferred distance)
                    if pref_dist==1 % proximity
                        dataTable.map_ambiguous(neuron) = false;
                    elseif pref_dist > 1 && ~any(p_map < alpha) % any(~(p_map(sig_dist+1==idx) < alpha))
                        dataTable.map{neuron} = 'ambiguous';
                        dataTable.map_ambiguous(neuron) = true;
                    elseif pref_dist > 1 % otherwise unabmiguous map cell
                        dataTable.map_ambiguous(neuron) = false;
                    end

                elseif dataTable.suppressed(neuron) 
                    dataTable.map{neuron} = 'suppressed';
                    dataTable.pref_dist(neuron) = pref_dist_sup;
                end
            end
        
        catch e
            disp(e)
        end
    end
end

% Fill in map labels for untuned and ambiguous subpopulations
for neuron = 1:size(dataTable,1)
    if isempty(dataTable.map{neuron})
        dataTable.map{neuron} = 'untuned';
    end
end
  
acceptable = sst.acceptable & ~dataTable.low_firing & ~dataTable.stats_err;
incl = dataTable.wallResp & sst.acceptable & ~dataTable.low_firing & ~dataTable.stats_err;
wall_tuned = dataTable.wallTuned & dataTable.wallResp;

disp(newline)

% Unit inclusion stats
disp([num2str(sum(acceptable)) ' out of ' num2str(size(dataTable,1)) ' units are acceptable'])
disp([num2str(sum(dataTable.low_firing)) ' out of ' num2str(size(dataTable,1)) ' units low firing (< 1 Hz)'])
disp([num2str(sum(dataTable.stats_err)) ' out of ' num2str(sum(sst.acceptable)) ' acceptable units errored on stats'])
disp([num2str(sum(acceptable)) ' units included in stats (acceptable and not low firing)'])

% Wall-responsive vs. wall-tuned
disp([num2str(sum(incl)) ' out of ' num2str(sum(acceptable)) ' acceptable units are wall-responsive'])
disp([num2str(sum(wall_tuned & incl)) ' out of ' num2str(sum(incl)) ' acceptable units are wall-tuned'])
disp([num2str(sum((~dataTable.wallTuned | ~dataTable.wallResp) & incl)) ' out of ' num2str(sum(incl)) ' acceptable units are untuned'])

disp([num2str(sum(wall_tuned & dataTable.activated & incl)) ' out of ' num2str(sum(wall_tuned & incl)) ' wall-tuned units are activated'])
disp([num2str(sum(wall_tuned & dataTable.suppressed & incl)) ' out of ' num2str(sum(wall_tuned & incl)) ' wall-tuned units are suppressed'])

disp([num2str(sum(dataTable.activated & incl)) ' out of ' num2str(sum(incl)) ' wall-resp neurons are activated'])
disp([num2str(sum(dataTable.suppressed & incl)) ' out of ' num2str(sum(incl)) ' wall-resp neurons are suppressed'])
disp([num2str(sum(dataTable.activated & dataTable.suppressed & incl)) ' out of ' num2str(sum(incl)) ' wall-resp neurons are activated and suppressed'])

disp([num2str(sum(~dataTable.suppressed & dataTable.activated & incl)) ' out of ' ...
    num2str(sum(incl)) ' wall-resp neurons are activated only cells'])
disp([num2str(sum(dataTable.suppressed & ~dataTable.activated & incl)) ' out of ' ...
    num2str(sum(incl)) ' wall-resp neurons are suppressed only cells'])

% Classify included units into groups

% Proximity (prefer closest distance)
proximity = strcmpi(dataTable.map, 'proximity') & wall_tuned & dataTable.activated & incl;
disp([num2str(sum(proximity)) '/' num2str(sum(incl)) ' (' num2str(100*sum(proximity)/sum(incl), '%.1f') '%) neurons are proximity cells'])

% Ambiguous (prefer non-closest distance but not significantly different from closest distance)
ambiguous = strcmpi(dataTable.map, 'ambiguous') & dataTable.map_ambiguous & wall_tuned & dataTable.activated & incl;
disp([num2str(sum(ambiguous)) '/' num2str(sum(incl)) ' (' num2str(100*sum(ambiguous)/sum(incl), '%.1f') '%) neurons are ambiguous map/proximity cells'])

% Map (prefer non-closest distance)
map = strcmpi(dataTable.map, 'map') & wall_tuned & dataTable.activated & incl;
disp([num2str(sum(map)) '/' num2str(sum(incl)) ' (' num2str(100*sum(map)/sum(incl), '%.1f') '%) neurons are map cells'])

% Suppressed only (not activated)
suppressed = dataTable.suppressed & ~dataTable.activated & wall_tuned & incl;
disp([num2str(sum(suppressed)) '/' num2str(sum(incl)) ' (' num2str(100*sum(suppressed)/sum(incl), '%.1f') '%) neurons are suppressed only cells'])

% Untuned
untuned = incl & strcmpi(dataTable.map, 'untuned');
disp([num2str(sum(untuned)) '/' num2str(sum(incl)) ' (' num2str(100*sum(untuned)/sum(incl), '%.1f') '%) neurons are untuned cells'])

disp([num2str(sum(untuned | suppressed | map | ambiguous | proximity)) '/' num2str(sum(incl)) ' total cells accounted for'])

% Save tuning table with stats
disp('Saving stats into data table ...')
save(fullfile(settings.data_folder, 'tuningTable_PrV_stats.mat'), 'dataTable')

%% Plot combined (naive + whisker cutting intact) population data
% ** Also saves "monotonic" unit indices for population decoding analysis

% Plot Population Tuning for dataset
clc; clearvars -except settings
% settings.dataset = 'naive'; % 'telc'; % 'laser'; % 
if strcmp(settings.dataset, 'naive')
    result_dir = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData';
    wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_naive';
    sst = matfile(fullfile(wv_stats, 'PrV_spike_stats_table_naive.mat')).sst;
elseif strcmp(settings.dataset, 'laser')
    result_dir = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\GtACR';
    wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_Laser';
    sst = matfile(fullfile(wv_stats, 'PrV_spike_stats_table_laser.mat')).sst;
elseif strcmp(settings.dataset, 'telc')
    result_dir = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\TeLC';
    wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_TeLC';
    sst = matfile(fullfile(wv_stats, 'PrV_spike_stats_table_telc.mat')).sst;
elseif strcmp(settings.dataset, 'onewhisker')
    result_dir = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\One_Whisker_Cut';
    wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_OneWhisker';
    sst = matfile(fullfile(wv_stats, 'PrV_spike_stats_table_onewhisker.mat')).sst;
elseif strcmp(settings.dataset, 'tungsten')
    result_dir = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\Tungsten';
    wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_Tungsten';
    sst = matfile(fullfile(wv_stats, 'PrV_spike_stats_table_tungsten.mat')).sst;    
end

%% Panel D, E: Plot combined (naive + whisker cutting intact) population heatmap and histogram
% Settings
use_interp = 1;
us = 20; % interpolation n-fold upsample
save_plots = 1;
fig_size = [200 200 300 800];

if exist('output_path')==1
    settings.input_path = settings.output_path;
end

% Load Combined Tuning Table with Stats Tests
dataTable = matfile(fullfile(result_dir,'tuningTable_PrV_stats.mat')).dataTable; 

% Mask unit inclusion (low-firing (<1 Hz) units excluded during stats)
mask = sst.acceptable & ~dataTable.low_firing & ~dataTable.stats_err & dataTable.wallTuned; % & ~dataTable.map_ambiguous;
t = dataTable;

disp('Plotting population tuning heatmaps ...')

% Square off table so all units have same number of distance values
max_distances = cell2mat(cellfun(@(x) max(round(x)), t.distanceMM, 'Uni',0));
min_distances = cell2mat(cellfun(@(x) min(round(x)), t.distanceMM, 'Uni',0));
max_distance = max(max_distances);
min_distance = min(min_distances);
fields = {'distanceMM','fr_mean','fr_sem','fr_mean_shuf','fr_sem_shuf',...
    'zfr_mean','zfr_sem','zfr_mean_shuf','zfr_sem_shuf'};
num_recs = size(t,1);
dists = min_distance : max_distance;
num_dist = numel(dists);

for i = 1:num_recs
    % Fill in nans for distances not recorded
    [~,mn] = min(abs(dists-min_distances(i)));
    [~,mx] = min(abs(dists-max_distances(i)));
    for f = 1:numel(fields)
        tmp = NaN(1,num_dist);
        try
            tmp(:,mn:mx) = t.(fields{f}){i};
        catch e
            ds = 2;
            tmp(:,mn:mx) = interp1(linspace(mn,mx,numel(t.(fields{f}){i})), t.(fields{f}){i}, mn:mx, 'cubic');
        end
        t.(fields{f}){i} = tmp;
    end
end

t_mask = t(mask,:);
fr_mean = cell2mat(t_mask.fr_mean);
if strcmp(settings.dataset, 'tungsten')
    act = (t_mask.activated | t_mask.suppressed);
else
    act = t_mask.activated;
end
sup = t_mask.suppressed;
ambiguous = t_mask.map_ambiguous;

% Filter by mean
fr_mean_filt = fr_mean;

% Next normalize to min and max for each unit
% **Should try implementing Z-Score (subtract baseline)

% Subtract min and normalize
min_fr = min(fr_mean_filt,[],2,'omitnan');
max_fr = max(fr_mean_filt,[],2,'omitnan');
fr_mean_norm = (fr_mean_filt - min_fr) ./ (max_fr  - min_fr);

if use_interp
    num_bins = size(fr_mean_norm,2);
    interp_dists = linspace(dists(1),dists(end),num_dist*us);
    for i = size(fr_mean_norm,1):-1:1
        fr_mean_norm_interp(i,:) = interp1(1:num_bins,fr_mean_norm(i,:),linspace(1,num_dist,num_dist*us),'linear');
    end
    fr_mean_norm = fr_mean_norm_interp;
    dists = interp_dists;
end

% Remove firing rates for distances outside the binning range (by setting to NaN)
if strcmp(settings.dataset, 'tungsten')
    dist_edges = 3.5:1:20.5;
else
    dist_edges = [6 8 10 12 14 16 18 20 22 24];
end
dist_centers = mean([dist_edges(1:end-1); dist_edges(2:end)]);
fr_mean_norm(:,dists<dist_edges(1) | dists > dist_edges(end)) = NaN;

% Finally, sort by index of max for activated population
[~,idx] = max(fr_mean_norm,[],2,'omitnan');
[~,order] = sort(idx,'ascend');
fr_mean_norm_sorted = fr_mean_norm(order,:);

% Plot population tuning as image, with rows rounded to nearest mm
[num_units,num_dist] = size(fr_mean_norm_sorted(act(order),:));
disp(['Using ' num2str(num_units) ' units in activated population'])
rows = 1:num_units;
cols = dists;

% Plot figure for activated cells
figure(10); clf;
fig_size = [2 2 6 16];
set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',fig_size)
colormap('parula')
imAlpha = ones(size(fr_mean_norm_sorted(act(order),:)));
imAlpha(isnan(fr_mean_norm_sorted(act(order),:))) = 0;
imagesc(cols, rows, fr_mean_norm_sorted(act(order),:), 'AlphaData',imAlpha);
box off; grid off; axis tight

set(gca, 'TickDir','out', 'Box','off', 'Color',[1 1 1], ...
    'Ydir','reverse', 'YTick',[1 num_units], ...
    'XTick',dist_centers, 'XTickLabelRotation',0)
view(2)
title([settings.region ' Activated Population Wall Tuning'])
ylabel('Unit number (far \rightarrow close)')
xlabel('Wall distance (mm)')
ylim([0.5 num_units+0.5])
xlim([dist_edges(1) dist_edges(end)])
cb = colorbar;
set(cb, 'TickDir','out')
ylabel(cb, 'Norm. Spike Rate', 'Rotation',270)
cb.Label.Position(1) = 4;
pause(0.01)

if save_plots
    filename = fullfile(settings.fig_dir, [settings.region '_activated_population_tuning_' settings.dataset]);
    savefig(gcf,filename);
    saveas(gcf,filename,'png');
    saveas(gcf,filename,'pdf');
    pause(0.01)
end

% Sort suppressed population by min
[~,idx] = min(fr_mean_norm,[],2,'omitnan');
[~,order] = sort(idx,'ascend');
fr_mean_norm_sorted = fr_mean_norm(order,:);

% Plot population tuning as image, with rows rounded to nearest mm
[num_units, num_dist] = size(fr_mean_norm_sorted(sup(order),:));
if num_units > 0
    disp(['Using ' num2str(num_units) ' units in suppressed population'])
    rows = 1:num_units;
    cols = dists;
    
    % Plot figure for suppressed cells
    figure(20); clf
    fig_size = [2 2 6 16];
    set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',fig_size)
    colormap('parula')
    imAlpha = ones(size(fr_mean_norm_sorted(sup(order),:)));
    imAlpha(isnan(fr_mean_norm_sorted(sup(order),:))) = 0;
    imagesc(cols, rows, fr_mean_norm_sorted(sup(order),:), 'AlphaData',imAlpha);
    box off; grid off; axis tight
    dist_centers = mean([dist_edges(1:end-1); dist_edges(2:end)]);
    set(gca, 'TickDir','out', 'Box','off', 'Color',[1 1 1], ...
        'Ydir','reverse', 'YTick',[1 num_units], ...
        'XTick',dist_centers, 'XTickLabelRotation',0)
    view(2)
    title([settings.region ' Suppressed Population Wall Tuning'])
    ylabel('Unit number (far \rightarrow close)')
    xlabel('Wall distance (mm)')
    ylim([0.5 num_units+0.5])
    xlim([dist_edges(1) dist_edges(end)])
    cb = colorbar;
    set(cb, 'TickDir','out')
    ylabel(cb, 'Norm. Spike Rate', 'Rotation',270)
    cb.Label.Position(1) = 4;
    pause(0.01)
    
    if save_plots
        filename = fullfile(settings.fig_dir, [settings.region '_suppressed_population_tuning_' settings.dataset]);
        savefig(gcf,filename);
        saveas(gcf,filename,'png');
        saveas(gcf,filename,'pdf');
        pause(0.01)
    end
end

% Bin heatmap to 9 distances
if strcmp(settings.dataset, 'tungsten')
    dist_edges = 3.5:1:20.5;
else
    dist_edges = [6 8 10 12 14 16 18 20 22 24];
end
dist_centers = mean([dist_edges(1:end-1); dist_edges(2:end)]);
wd = cell2mat(t_mask.distanceMM);
fr_mean_9dist = NaN(size(fr_mean,1),numel(dist_centers));
for u = 1:size(fr_mean,1)
    [~,~,dbins] = histcounts(wd(u,:), dist_edges);
    for d = 1:numel(dist_centers)
        fr_mean_9dist(u,d) = mean(fr_mean(u,dbins==d),'omitnan');
    end
end

% Normalize to range [0 1]
fr_mean_9dist_norm = (fr_mean_9dist - min(fr_mean_9dist,[],2)) ./ ...
    (max(fr_mean_9dist,[],2) - min(fr_mean_9dist,[],2));

% Sort units by max idx
[~,distance_bins] = max(fr_mean_9dist_norm,[],2);
[~,sortIdx] = sort(distance_bins,'ascend');
fr_mean_9dist_sorted = fr_mean_9dist_norm(sortIdx,:);
actSort = act(sortIdx);

% Get min index for suppression units
[~,distance_mins] = min(fr_mean_9dist_norm,[],2);

% Plot binned heatmap for activated units, sorted by max idx
figure(333); clf
fig_size = [2 2 6 16];
set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',fig_size)
C = flipud(fr_mean_9dist_sorted(actSort,:));
imagesc('XData',dist_centers, 'YData',1:size(C,1), 'CData',C)
dist_centers = mean([dist_edges(1:end-1); dist_edges(2:end)]);
set(gca, 'TickDir','out', 'Box','off', 'Color',[1 1 1], ...
    'YTick',[1 size(C,1)], ...
    'XTick',dist_centers, 'XTickLabelRotation',0)
cb = colorbar;
set(cb, 'TickDir','out')
xlim([min(dist_centers)-0.5 max(dist_centers)+0.5])
ylim([0.5 size(C,1) + 0.5])
title('PrV Activated Population')
ylabel('Neurons')
xlabel('Wall Distance (mm)')
pause(0.01)

if save_plots
    filename = fullfile(settings.fig_dir, [settings.region '_activated_population_tuning_binned_' settings.dataset]);
    savefig(gcf,filename);
    saveas(gcf,filename,'png');
    saveas(gcf,filename,'pdf');
    pause(0.01)
end

try
    % Sort units by min idx
    [~,b] = min(fr_mean_9dist_norm,[],2);
    [~,sortIdx] = sort(b,'ascend');
    fr_mean_9dist_sorted = fr_mean_9dist_norm(sortIdx,:);
    
    % Plot binned heatmap for suppressed units, sorted by max idx
    figure(334); clf
    fig_size = [2 2 6 16];
    set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',fig_size)
    C = flipud(fr_mean_9dist_sorted(sup(sortIdx),:));
    imagesc('XData',dist_centers, 'YData',1:size(C,1), 'CData',C)
    dist_centers = mean([dist_edges(1:end-1); dist_edges(2:end)]);
    set(gca, 'TickDir','out', 'Box','off', 'Color',[1 1 1], ...
        'YTick',[1 size(C,1)], ...
        'XTick',dist_centers, 'XTickLabelRotation',0)
    cb = colorbar;
    set(cb, 'TickDir','out')
    xlim([min(dist_centers)-0.5 max(dist_centers)+0.5])
    ylim([0.5 size(C,1) + 0.5])
    title('PrV Suppressed Population')
    ylabel('Neurons')
    xlabel('Wall Distance (mm)')
    pause(0.01)
    
    if save_plots
        filename = fullfile(settings.fig_dir, [settings.region '_suppressed_population_tuning_binned_' settings.dataset]);
        savefig(gcf,filename);
        saveas(gcf,filename,'png');
        saveas(gcf,filename,'pdf');
        pause(0.01)
    end
end

% Get indices of units with tuning that is monotonic inversely proportional vs. not
% (out of all tuned & activated units with acceptable spike stats)
bin_num = 1; % 2;
distances = dist_centers(1:bin_num);
proximity = sst.acceptable & ~dataTable.low_firing & ~dataTable.stats_err & dataTable.wallTuned & strcmp(dataTable.map, 'proximity');
map = sst.acceptable & ~dataTable.low_firing & ~dataTable.stats_err & dataTable.wallTuned & strcmp(dataTable.map, 'map') & ~dataTable.map_ambiguous;

suppressed = sst.acceptable & ~dataTable.low_firing & ~dataTable.stats_err & dataTable.wallTuned & dataTable.suppressed;

disp([num2str(sum(proximity)) ' units in Proximity subpop pref <' num2str(dist_edges(bin_num+1)) ' mm '])
disp([num2str(sum(map)) ' units in Map subpop pref >' num2str(dist_edges(bin_num+1)) ' mm'])
disp([num2str(sum(suppressed)) ' units in Suppressed subpop'])

% Plot histogram of preferred distances for activated population
colors = flipud(turbo(numel(dist_centers)));
sz = 8; % size of figure dimensions, in cm
dist_dist = histcounts(dist_centers(distance_bins(act)), dist_edges);

figure(11); clf
figpos = [sz sz sz sz];
set(gcf, 'Color','w', 'Units','centimeters', 'Position',figpos)
b = bar(dist_centers, dist_dist, 'EdgeColor','none');
b.FaceColor = 'flat';
for d = 1:size(colors,1)
    b.CData(d,:) = colors(d,:);
end
set(gca, 'TickDir','out', 'Box','off')
title(['Pref Dist (Max)' newline 'Activated Population'])
ylabel('Neurons')
xlabel('Distance (mm)')
xlim([dist_centers(1)-1.5 dist_centers(end)+1])
pause(0.01)

if save_plots
    filename = fullfile(settings.fig_dir, [settings.region '_naive_act_pop_pref_dist_hist_' settings.dataset]);
    savefig(gcf,filename);
    saveas(gcf,filename,'png');
    saveas(gcf,filename,'pdf');
    pause(0.01)
end

% Plot histogram of preferred wall distances for suppressed
try
    dist_dist = histcounts(dist_centers(distance_mins(sup)), dist_edges);
    figure(21); clf
    figpos = [sz+sz sz sz sz];
    set(gcf, 'Color','w', 'Units','centimeters', 'Position',figpos)
    b = bar(dist_centers, dist_dist, 'EdgeColor','none');
    b.FaceColor = 'flat';
    for d = 1:size(colors,1)
        b.CData(d,:) = colors(d,:);
    end
    set(gca, 'TickDir','out', 'Box','off')
    title(['Pref Dist (Min)' newline 'Suppressed Population'])
    ylabel('Neurons')
    xlabel('Distance (mm)')
    xlim([dist_centers(1)-1.5 dist_centers(end)+1])
    pause(0.01)
    
    if save_plots
        filename = fullfile(settings.fig_dir, [settings.region '_naive_sup_pop_pref_dist_hist_' settings.dataset]);
        savefig(gcf,filename);
        saveas(gcf,filename,'png');
        saveas(gcf,filename,'pdf');
        pause(0.01)
    end
end

% Also plot histogram of preferred wall distances for all neurons
try
    dist_dist = histcounts(dist_centers(distance_bins), dist_edges);
    figure(31); clf
    figpos = [sz+2*sz sz sz sz];
    set(gcf, 'Color','w', 'Units','centimeters', 'Position',figpos)
    b = bar(dist_centers, dist_dist, 'EdgeColor','none');
    b.FaceColor = 'flat';
    for d = 1:size(colors,1)
        b.CData(d,:) = colors(d,:);
    end
    set(gca, 'TickDir','out', 'Box','off')
    title(['Pref Dist (Max)' newline 'All Population'])
    ylabel('Neurons')
    xlabel('Distance (mm)')
    xlim([dist_centers(1)-1.5 dist_centers(end)+1])
    pause(0.01)
    
    if save_plots
        filename = fullfile(settings.fig_dir, [settings.region '_naive_all_pop_pref_dist_hist_' settings.dataset]);
        savefig(gcf,filename);
        saveas(gcf,filename,'png');
        saveas(gcf,filename,'pdf');
        pause(0.01)
    end
end

%% Panel F: Define proximity (monotonic) vs. map cells for population decoding and plot pie chart

save_plots = 1;

% Also display stats, activated vs. suppressed Venn diagram, and Map/Proximity/Suppressed Only Pie chart
t_mask = t(sst.acceptable & ~t.low_firing,:);
tuned = t_mask.wallTuned | (t_mask.wallTuned & ~t_mask.wallResp);
proximity = strcmp(t_mask.map,'proximity');
map = strcmp(t_mask.map,'map');
ambiguous = strcmp(t_mask.map,'ambiguous');
unit_names = cellfun(@(x,y) [x '_unit_' num2str(y)], t_mask.recID, num2cell(t_mask.unit,2), 'Uni',0);

% Save logical indices of population for population decoding (and other metadata)
decoding_dir = 'E:\Code\Analysis\NeuralDecoder';
filename = [settings.dataset '_proximity_map_population'];
save(fullfile(decoding_dir, filename), 'tuned', 'proximity', 'map', 'ambiguous', 'distances', 'distance_bins', 'unit_names')

% Define total population (responsive units = 768)
incl = dataTable.wallResp & sst.acceptable & ~dataTable.low_firing & ~dataTable.stats_err;

% Plot pie chart (total = responsive units, slices = untuned, proximity, map, suppressed)
proximity_slice = sum(incl & dataTable.wallTuned & strcmp(dataTable.map,'proximity')); % monotonic
ambiguous_slice = sum(incl & dataTable.wallTuned & strcmp(dataTable.map,'ambiguous')); % ambiguous map/proximity cells
map_slice = sum(incl & dataTable.wallTuned & strcmp(dataTable.map,'map'));
suppressed_slice = sum(incl & dataTable.wallTuned & strcmp(dataTable.map,'suppressed'));
% untuned_slice = sum(incl) - suppressed_slice - map_slice - proximity_slice - ambiguous_slice;
% untuned_slice = sum(incl & ~dataTable.wallTuned | (dataTable.wallTuned & ~dataTable.wallResp)); % with unresponsive but wall-tuned
untuned_slice = sum(incl & ~dataTable.wallTuned); % with unresponsive but wall-tuned
total = sum(incl)

slices = [proximity_slice, ambiguous_slice, map_slice, suppressed_slice, untuned_slice];
disp(['Piechart slices account for ' num2str(sum(slices)) '/' num2str(total) ' units'])
names = {'Proximity', 'Ambiguous', 'Map', 'Suppressed', 'Untuned'};

% Save data to mat file
filename = fullfile(settings.fig_dir, ['piechart_source_data_' settings.dataset]);
save(filename,'slices','names')

for n = 1:numel(names)
    labels{n} = [names{n} ' (' num2str(round(slices(n)/total*100,1)) '%)'];
    disp([labels{n} ' ' num2str(round(slices(n))) '/' num2str(total)])
end

figure(314); clf
set(gcf, 'Color','w', 'Renderer','painters')
p = pie(slices, labels);
patchHandles = findobj(p, 'Type','Patch'); 
newColors = [1 0.6 0.6; 1 0.8 0.6; 1 1 0.6; 0.6 0.6 1; 0.6 0.6 0.6];
set(patchHandles, {'FaceColor'}, mat2cell(newColors, ones(size(newColors,1),1), 3))
set(gca, 'tickDir','out')
title([upper(settings.dataset(1)) lower(settings.dataset(2:end)) ' (N=' num2str(sum(slices)) ')'])

if save_plots
    filename = fullfile(settings.fig_dir, ['piechart_' settings.dataset]);
    saveas(gcf, filename, 'pdf')
    saveas(gcf, filename, 'png')
end
