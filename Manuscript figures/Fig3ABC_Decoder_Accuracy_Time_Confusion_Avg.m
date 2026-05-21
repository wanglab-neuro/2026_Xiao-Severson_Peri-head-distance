%% PrV naive data population decoding
clear; clc

% Initialize decoder
cd('E:\Code\Analysis\NeuralDecoder')

% add the path to the NDT so add_ndt_paths_and_init_rand_generator can be called
toolbox_basedir_name = 'ndt.1.0.4/';
addpath(toolbox_basedir_name);
 % add the NDT paths using add_ndt_paths_and_init_rand_generator
add_ndt_paths_and_init_rand_generator


%% Rasterize the data for each neuron
save_dir = 'E:\Code\Analysis\NeuralDecoder\ndt.1.0.4\Wenxi_PrV_9dist_raster_data';

sst = matfile('E:\PrV_Wall_Recordings\WaveformStats_naive\PrV_spike_stats_table_naive.mat').sst;
dataTable = matfile('E:\PrV_Wall_Recordings\Analysis\IntermediateData\tuningTable_PrV_anova.mat').dataTable;
psth = matfile('E:\PrV_Wall_Recordings\Analysis\IntermediateData\PSTH_naive.mat').psth;
psth.dataTable.wall_tuned = dataTable.wallTuned;
psth.dataTable.acceptable = sst.acceptable;
psth.dataTable.low_firing = dataTable.low_firing;

% Mask included units
included_units = ~dataTable.low_firing & sst.acceptable == 1;
psth.dataTable = psth.dataTable(included_units, :);

% shape data into raster format
start = 0; % ms
stop = 14000; % ms
binsize = 1; % ms

num_units = size(psth.dataTable,1);
for unit = 1:num_units
    unit_name = psth.dataTable.Properties.RowNames{unit};
    disp(['Unit ' num2str(round(unit)) '/' num2str(num_units)])
    raster_data = [];
    raster_labels.stimulus_ID = [];
    dist_centers = [7 9 11 13 15 17 19 21 23];
    for d = 1:length(dist_centers)
        distance = ['Dist_' num2str(dist_centers(d)) 'mm'];
        distances = unique(psth.dataTable.wall_distance{unit,1});
        [min_dist, nearest_neighbor] = min(abs(distances - dist_centers(d)));
        nearest_distance = distances(nearest_neighbor);
        idx = find(psth.dataTable.wall_distance{unit,1} == nearest_distance);
    
        psth_bin1ms.(distance) = spike2eventRasteandPSTH_NP( ...
                    psth.dataTable.spike_times{unit,1}, ...
                    psth.dataTable.pass_starts{unit,1}(idx), ...
                    binsize, start, stop);
    
        raster_data = cat(1, raster_data, psth_bin1ms.(distance).scmatrix);
        raster_labels.stimulus_ID = cat(2, raster_labels.stimulus_ID, repmat({distance},1,numel(idx)));
    end
    disp([num2str(length(raster_labels.stimulus_ID)) ' trials'])

    % Save data
    save(fullfile(save_dir, unit_name), 'raster_data', '-v7.3');
    save(fullfile(save_dir, unit_name), "raster_labels", '-append');
    
    % Append raster size info
    unit_num = strcmp(sst.Properties.RowNames,unit_name);
    parts = split(unit_name,'_');
    raster_site_info.session_ID = strjoin(parts(1:4),'_');
    raster_site_info.recording_channel = sst.ch_id(unit_num);
    raster_site_info.unit = str2double(parts{end});
    raster_site_info.alignment_event_time = 1;
    
    save(fullfile(save_dir, unit_name), "raster_site_info", '-append');
end

%% Bin the data for each neuron
cd('E:\Code\Analysis\NeuralDecoder\ndt.1.0.4')
raster_file_directory_name = 'Wenxi_PrV_9dist_raster_data/';
save_prefix_name = 'Binned_Wenxi_PrV_9dist_raster_data';
bin_width = 150;
step_size = 50;

create_binned_data_from_raster_data(raster_file_directory_name, save_prefix_name, bin_width, step_size);

%% Subsample 'activated' population
% For equalizing subpopulation size (monotonic vs. non-monotonic)
% And showing performance vs. population size
clear; clc

% Fix random seed
rng('Default')
rng(42)

% Settings
decoding_dir = 'E:\Code\Analysis\NeuralDecoder';
cd('E:\Code\Analysis\NeuralDecoder\ndt.1.0.4')
raster_file_directory_name = 'Wenxi_PrV_9dist_raster_data/';
save_prefix_name = 'Binned_Wenxi_PrV_9dist_raster_data';

filename = 'naive_proximity_map_population_mixed.mat';
datetoday = '20251207';

% Load tables for indexing units
sst = matfile('E:\PrV_Wall_Recordings\WaveformStats_naive\PrV_spike_stats_table_naive.mat').sst;
t = matfile('E:\PrV_Wall_Recordings\Analysis\IntermediateData\tuningTable_PrV_anova.mat').dataTable;
t_mask = t(t.wallResp & sst.acceptable & ~t.low_firing,:); % 1141 with acceptable spike stats

activated = shiftdim(t_mask.activated & t_mask.wallTuned); % 383/1140 tuned and activated
suppressed = shiftdim(t_mask.suppressed & t_mask.wallTuned);
proximity = strcmp(t_mask.map, 'proximity');
map = strcmp(t_mask.map, 'map');
ambiguous = strcmp(t_mask.map, 'ambiguous');
untuned = strcmp(t_mask.map, 'untuned');

untuned_pop = untuned;
proximity_pop = activated & proximity;
ambiguous_pop = activated & ambiguous;
map_pop = activated & map & ~ambiguous;
shuffled_pop = proximity_pop | map_pop;
suppressed_pop = strcmp(t_mask.map, 'suppressed'); % suppressed only

disp(['Untuned: ' num2str(sum(untuned_pop))])
disp(['Proximity: ' num2str(sum(proximity_pop))])
disp(['Ambiguous: ' num2str(sum(ambiguous_pop))])
disp(['Map: ' num2str(sum(map_pop))])
disp(['Shuffled: ' num2str(sum(shuffled_pop))])
disp(['Suppressed: ' num2str(sum(suppressed_pop))])
% 592 untuned, 263 proximity, 139 map, 102 suppressed, 402 shuffled

% Create new folders for each subpopulation
population_names = {'untuned', 'proximity', 'map', 'mixed', 'suppressed', 'shuffled'};
for p = 1:numel(population_names)
    subpop_folder = fullfile(decoding_dir,'ndt.1.0.4','subpop',population_names{p});
    if ~isfolder(subpop_folder); mkdir(subpop_folder); end
end

population_unit_inds = {
    find(untuned_pop), ...
    find(proximity_pop), ...
    find(map_pop), ...
    find(map_pop | proximity_pop), ...
    find(suppressed_pop), ...
    find(shuffled_pop), ...
    };
population_size = cellfun(@(x) numel(x), population_unit_inds, 'Uni',1)
[min_population, min_idx] = min(population_size);
[pop_sizes] = [min_population]; % [133 80 40 20 10 5]; % 
smallest_population = population_names{min_idx}
disp(min_population)

num_subsamples = 10;
bin_width = 150;
step_size = 50;

%% Load subsampled data for each population
datafile = 'Binned_Wenxi_PrV_9dist_raster_data_150ms_bins_50ms_sampled.mat';
data = load(fullfile(decoding_dir,'ndt.1.0.4',datafile));

save_subsamples = true; % false;
mixed_map_size = 36;
mixed_prox_size = 82;

% Get unit names from binned data (order is not in "natsort" order)
psth_unit_names = t_mask.Properties.RowNames;
psth_class = t_mask.map;
binned_unit_names = cellfun(@(x,u) [x '_' num2str(u)], data.binned_site_info.session_ID', num2cell(data.binned_site_info.unit), 'Uni',0);

% Subsample data n times and save subsampled data for decoding
% population_names = {'untuned', 'proximity', 'map', 'mixed', 'suppressed', 'shuffled'};
for p = 1:numel(population_names)
    disp(['Sub-sampling binned data from population: ' population_names{p}])
    p_inds = population_unit_inds{p};
    subpop_folder = fullfile(decoding_dir,'ndt.1.0.4','subpop',population_names{p});
    for k = 1:numel(pop_sizes)
        for n = 1:num_subsamples

            if strcmp(population_names{p}, 'mixed')
                % Take random sample of with aproximate ration observed in
                % population (~2:1 prox:map)
                map_inds = p_inds(strcmp(psth_class(p_inds),'map'));
                map_units = datasample(map_inds, mixed_map_size, 'Replace',false);
                prox_inds = p_inds(strcmp(psth_class(p_inds),'proximity'));
                prox_units = datasample(prox_inds, mixed_prox_size, 'Replace',false);
                rand_units = [map_inds; prox_inds];
            else
                % Take random sample of units equal to population size (without replacement)
                rand_units = datasample(p_inds, pop_sizes(k), 'Replace',false);
                rand_unit_names = psth_unit_names(rand_units);
                rand_unit_class = psth_class(rand_units);
            end
        
            % Subsample cells' "binned_data" (PSTH)
            rand_binned_units = ismember(binned_unit_names, rand_unit_names);
            binned_data = data.binned_data(rand_binned_units);

            % If shuffled control, shuffle the trial data with respect to stimulus ID
            if strcmp(population_names{p}, 'shuffled')
                for c = 1:numel(binned_data)
                    rand_trials = datasample(1:size(binned_data{c},1), size(binned_data{c},1), 'Replace',false);
                    binned_data{c} = binned_data{c}(rand_trials,:);
                end
            end
            
            % Subsample cells' "binned_labels.stimulus_ID" to match stims with responses
            binned_labels.stimulus_ID = data.binned_labels.stimulus_ID(rand_binned_units);
            
            % Subsample cells' "binned_site_info"
            binned_site_info.session_ID = data.binned_site_info.session_ID(rand_binned_units);
            binned_site_info.recording_channel = data.binned_site_info.recording_channel(rand_binned_units);
            binned_site_info.unit = data.binned_site_info.unit(rand_binned_units);
            binned_site_info.alignment_event_time = data.binned_site_info.alignment_event_time(rand_binned_units);
            binned_site_info.binning_parameters = data.binned_site_info.binning_parameters;
        
            % Populate the output binned data filename
            filename = [datetoday '_' save_prefix_name '_' num2str(bin_width) 'ms_bins_' ...
                num2str(step_size) 'ms_sampled_' ...
                population_names{p} '_' ...
                num2str(pop_sizes(k)) ...
                'units_' num2str(n)];

            % Save binned data
            if save_subsamples
                disp(['Saving ' filename '...'])
                save(fullfile(subpop_folder,filename), 'binned_data', 'binned_labels', 'binned_site_info')
            elseif n == 1
                disp('Not saving subsamples')
            end
        end
    end
end

% Test - plot rasters for subsampled data
% clear; clc; close all
% base_dir = 'E:\Code\Analysis\NeuralDecoder\ndt.1.0.4';
% p = 'untuned'; % 'mixed'; % 'suppressed'; % 'map'; % 'proximity'; % 'shuffled'; % 
% datetoday = '20251207';
% filenames = dir(fullfile(base_dir, 'subpop', p, '*Binned_Wenxi*'));
% f = 1;
% data = load(fullfile(base_dir, 'subpop', p, filenames(f).name));
% for u = 1:100
%     figure(u); clf
%     imagesc(data.binned_data{u})
%     pause(0.001)
% end

%% Run population decoding for subsampled subpopulation data
clear; clc

population_names = {'untuned', 'proximity', 'map', 'mixed', 'suppressed', 'shuffled'};
datetoday = '20251207';
num_subsample = 10;

for p = 1:numel(population_names)
    clearvars -except p population_names datetoday num_subsample
    pop = population_names{p}

    toolbox_basedir_name = 'E:\Code\Analysis\NeuralDecoder\ndt.1.0.4';
    cd(toolbox_basedir_name)
    addpath(toolbox_basedir_name);
    add_ndt_paths_and_init_rand_generator
    
    subpop_folder = 'E:\Code\Analysis\NeuralDecoder\ndt.1.0.4\subpop';
    save_prefix_name = [datetoday '_Binned_Wenxi_PrV_9dist_raster_data'];

    subsample_files = dir(fullfile(subpop_folder, pop, [save_prefix_name '*']));
    subpop_folder = fullfile(subpop_folder,pop);
    for f = 1:min([num_subsample numel(subsample_files)])
        tic
        
        % Define binned data filename
        fn = subsample_files(f).name
        disp(['Decoding ' fn])
        bd = fullfile(subpop_folder, fn);
        % load(binned_data);

        % Copy binned data to toolbox folder
        copyfile(bd, toolbox_basedir_name)
        pause(0.1)
         
        % Number of CV splits is the fewest number of repeats for trial condition
        num_cv_splits = 15; % 15*9 = 135
        
        % the name of the file that has the data in binned-format
        binned_format_file_name = fn;
         
        % will decode the identity of which object was shown (regardless of its position)
        specific_label_name_to_use = 'stimulus_ID';
        ds = basic_DS(binned_format_file_name, specific_label_name_to_use, num_cv_splits);
        
        % Create feature preprocessor to z-score normalize each neuron
        % note that the FP objects are stored in a cell array 
        % which allows multiple FP objects to be used in one analysis
        the_feature_preprocessors{1} = zscore_normalize_FP;
        
        % Create the classifier (CL) object
        the_classifier = max_correlation_coefficient_CL;
        
        % create the cross-validator (CV) object
        the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
         
        % set how many times the outer 'resample' loop is run
        % generally we use more than 2 resample runs which will give more accurate results
        % but to save time in this tutorial we are using a small number.
        the_cross_validator.num_resample_runs = 5;
        
        % Run population decoding analysis
        DECODING_RESULTS = the_cross_validator.run_cv_decoding;
        save_file_name = ['DECODING_RESULTS_' fn(1:end-4)];
        save(save_file_name, 'DECODING_RESULTS');
        pause(0.001)

        delete(fullfile(toolbox_basedir_name,fn))
        pause(0.001)
        toc
    end
end

%% Panel B: Subpopulation Confusion Matrices + Accuracy vs. Time

% Load decoding results across distances
% clear; clc
decode_dir = 'E:\Code\Analysis\NeuralDecoder\ndt.1.0.4';
cd(decode_dir)

population_names = {'Untuned', 'Proximity', 'Map', 'Mixed', 'Suppressed', 'Shuffled'}; % (mixed = map + ambiguous?)
rundate = '20250827';
rundate2 = '20251207';
n = '118'; % number of units in smallest population to subsample

% save_file_name = 'DECODING_RESULTS_Binned_Wenxi_PrV_9dist_raster_data_150ms_bins_50ms_sampled_activated-non-mono_169units_1.mat';
save_file_names = {
    ['DECODING_RESULTS_' rundate2 '_Binned_Wenxi_PrV_9dist_raster_data_150ms_bins_50ms_sampled_untuned_' n 'units'];
    ['DECODING_RESULTS_' rundate '_Binned_Wenxi_PrV_9dist_raster_data_150ms_bins_50ms_sampled_proximity_' n 'units'];
    ['DECODING_RESULTS_' rundate '_Binned_Wenxi_PrV_9dist_raster_data_150ms_bins_50ms_sampled_map_' n 'units'];
    ['DECODING_RESULTS_' rundate2 '_Binned_Wenxi_PrV_9dist_raster_data_150ms_bins_50ms_sampled_mixed_' n 'units'];
    ['DECODING_RESULTS_' rundate2 '_Binned_Wenxi_PrV_9dist_raster_data_150ms_bins_50ms_sampled_suppressed_' n 'units'];
    ['DECODING_RESULTS_' rundate '_Binned_Wenxi_PrV_9dist_raster_data_150ms_bins_50ms_sampled_shuffled_' n 'units']};

accuracy_by_distance = cell(1, numel(save_file_names));
error_by_distance = cell(1, numel(save_file_names));
for f = numel(save_file_names):-1:1
    % Get all subsamples with the same filename pattern
    file_names = dir(fullfile(decode_dir,[save_file_names{f} '*']));
    for ff = numel(file_names):-1:1
        % Load decoding results
        save_file_name = fullfile(decode_dir,file_names(ff).name);
        DECODING_RESULTS = matfile(save_file_name).DECODING_RESULTS;
        
        % Sort the labels by distance (sorted alphabetically by default in NDT)
        labels = DECODING_RESULTS.DS_PARAMETERS.label_names_to_use;
        parts = squeeze(split(labels,'_'));
        label_dists = cellfun(@(x) str2double(x(1:end-2)), parts(:,end), 'Uni',1);
        [~,sortOrder] = sort(label_dists,'ascend'); % sort short to far

        % Compute error matrix (error between decoded and actual distance)
        [distance_matrix_actual, distance_matrix_predicted] = meshgrid(label_dists(sortOrder));
        distance_matrix_error = distance_matrix_actual - distance_matrix_predicted;
        
        % Read data from the confusion matrix
        confusion_matrix = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix; % 9x9xTimepoints
        clear acc_by_distance err_by_distance
        for t = size(confusion_matrix,3):-1:1
            % Get confusion matrix and sort by increasing distance labels
            cm = squeeze(confusion_matrix(sortOrder,sortOrder,t));

            % Compute False Positive (FP), False Negative (FN), True Positive (TP), and True Negative (TN) Rates
            FP = shiftdim(sum(cm, 1, 'omitnan')' - diag(cm));
            FN = shiftdim(sum(cm, 2, 'omitnan') - diag(cm));
            TP = shiftdim(diag(cm));
            TN = shiftdim(sum(sum(cm,'omitnan'),'omitnan') - FP - FN - TP);
            
            % Compute True Positive Rate (TPR) and False Positive Rate (FPR)
            TPR = TP ./ (TP + FN); % True Positive Rate ** Best Accuracy Metric **
            FPR = FP ./ (FP + TN);
            
            % Compute classification accuracy
            acc_by_distance(:,t) = TPR;

            % Compute overall distance error
            err_by_distance(:,t) = sum(abs(cm .* distance_matrix_error),1) ./ sum(cm,1);
        end
        % Convert accuracy to percent
        accuracy_by_distance{f}(:,:,ff) = acc_by_distance;
        accuracy_by_distance{f}(:,:,ff) = 100*accuracy_by_distance{f}(:,:,ff);
        
        % Also get total accuracy from data struct
        overall_accuracy(:,ff) = 100*diag(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.mean_decoding_results);

        % Get average error (in mm)
        error_by_distance{f}(:,:,ff) = err_by_distance;

        % Plot confusion matrix
        % Plot Population Decoding Summary
        x = DECODING_RESULTS.DS_PARAMETERS.binned_site_info.binning_parameters.the_bin_start_times';
        summary_period = [3000 11000]; % 3-11 second wall-pass period
        summary_bins = x >= summary_period(1) & x <= summary_period(2);
        cm = confusion_matrix(sortOrder,sortOrder,summary_bins);
        cm = cm ./ sum(cm,1);
        CM(:,:,ff) = mean(cm,3);
    end
        
    num_dist = size(accuracy_by_distance{1},1);
    colors = flipud(turbo(num_dist));
    accuracy{f} = overall_accuracy;

    figure(7110+f); clf
    figpos = [2 2 5 5];
    set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',figpos)
    X = unique(distance_matrix_actual);
    Y = unique(distance_matrix_predicted);
    C = mean(CM,3);
    hold on
    imagesc('XData',X, 'YData',Y, 'CData',C)
    for row = 1:numel(Y)
        for col = 1:numel(X)
            if row == col
                txt = strip(sprintf('%.2f', round(C(row,col),2)), 'left', '0');
                text(Y(col), X(row), txt, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', ...
                    'Color','w', ...
                    'FontSize',6); 
            end
        end
    end
    axis square
    set(gca, 'TickDir','out', 'Box','off', 'XTick',X, 'YTick',Y)
    xlabel('Actual Distance')
    ylabel('Predicted Distance')
    colormap('copper')
    caxis([0 1])
    cb = colorbar;
    set(cb, 'TickDir','out')
    ylabel(cb, 'Proportion')
    xlim([min(X)-1 max(X)+1])
    ylim([min(Y)-1 max(Y)+1])
    title(population_names{f})
    pause(0.001)
    
    fig_dir = 'E:\PrV_Wall_Recordings\Figures';
    filename = fullfile(fig_dir, ['decoder_confusion_matrix_' population_names{f}]);
    saveas(gcf, filename, 'pdf')
    exportgraphics(gcf, [filename '.png'], 'Resolution',600)
    
    % figure(10+f); clf
    % figpos = [5 5 16 8];
    % set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',figpos)
    % hold on
    % x = DECODING_RESULTS.DS_PARAMETERS.binned_site_info.binning_parameters.the_bin_start_times;
    % chance = 100/num_dist;
    % plot(x, repmat(chance,1,numel(x)), '--', 'Color',[0.5 0.5 0.5])
    % for d = 1:num_dist
    %     plot(x, mean(squeeze(accuracy_by_distance{f}(d,:,:)),2,'omitnan')', 'Color',colors(d,:))
    % end
    % set(gca, 'TickDir','out', 'Box','off')
    % ylim([0 100])
    % xlabel('Time (ms)')
    % ylabel('Classifier Accuracy (%)')
    % pause(0.001)
end

% Panel S5A: Distance Error Matrix
figure(16453); clf
figpos = [2 2 5 5];
set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',figpos)
X = unique(distance_matrix_actual);
Y = unique(distance_matrix_predicted);
C = distance_matrix_error';
imagesc('XData',X, 'YData',Y, 'CData',C)
axis square
set(gca, 'TickDir','out', 'Box','off', 'XTick',X, 'YTick',Y)
xlabel('Actual Distance')
ylabel('Predicted Distance')
colormap(flipud(cbrewer2('RdBu',256)))
cb = colorbar;
set(cb, 'YTick',min(min(C)):4:max(max(C)), 'TickDir','out')
ylabel(cb, 'Distance Error')
xlim([min(X)-1 max(X)+1])
ylim([min(Y)-1 max(Y)+1])

fig_dir = 'E:\PrV_Wall_Recordings\Figures';
filename = fullfile(fig_dir, 'decoder_error_matrix');
saveas(gcf, filename, 'pdf')
exportgraphics(gcf, [filename '.png'], 'Resolution',600)

%% Panel A: Plot Population Decoding Summary vs. Time

% Overall decoder accuracy vs time, averaged across all distances
cmap = flipud(turbo(256));
cmap(:,end+1) = 1; % alpha
colors = zeros(size(accuracy,2),4);
colors(strcmp(population_names,'Proximity'),:) = cmap(1,:);
colors(strcmp(population_names,'Untuned'),:) = [0 0 0 1];
colors(strcmp(population_names,'Map'),:) = [0 0 0 1];
colors(strcmp(population_names,'Mixed'),:) = [1 0 1 1];
colors(strcmp(population_names,'Suppressed'),:) = [0 0 1 1];
colors(strcmp(population_names,'Shuffled'),:) = [0.5 0.5 0.5 1];

figure(44); clf
figpos = [5 5 16 8];
set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',figpos)
hold on
x = DECODING_RESULTS.DS_PARAMETERS.binned_site_info.binning_parameters.the_bin_start_times';
clear h
for f = numel(accuracy):-1:1
    y = mean(accuracy{f},2); % (time bins, n_models)
    err = bootci(1000, @(x) mean(x), squeeze(accuracy{f}')); % 95% CI
    err = mean(abs([y'; y']-err)); % make one-sided
    MPlot.ErrorShade(y', err, x'/1000, 'Color',colors(f,1:3)); % shade SEM/CI
    h(f) = plot(x/1000, y, 'Color',colors(f,:), 'LineWidth',0.75);
end
set(gca, 'TickDir','out', 'Box','off')
ylim([0 100])
xlabel('Time (s)')
ylabel('Classifier Accuracy (%)')
pause(0.001)
legend(h, population_names, 'Location','northeast', 'Box','off')

fig_dir = 'E:\PrV_Wall_Recordings\Figures';
filename = fullfile(fig_dir, 'decoding_accuracy_proximity_map_vs_time');
saveas(gcf, filename, 'pdf')
exportgraphics(gcf, [filename '.png'], 'Resolution',600)

% Plot Population Decoding Summary
summary_period = [3000 11000]; % 3-11 second wall-pass period

% Get logical index for summary period (bins during pass stimulus)
summary_bins = x >= summary_period(1) & x <= summary_period(2);

% Get distances
num_dist = size(accuracy_by_distance{1},1);
chance = 100/num_dist;

%% Panel C (left): Average Decoding Accuracy vs. Distance for each subpopulation 
% Compute accuracy-by-distance for summary period
y = cell2mat(cellfun(@(x) mean(x(:,summary_bins,:),2,'omitnan'), accuracy_by_distance, 'Uni',0));

% Get colormap
cmap = flipud(turbo(256));
colors = zeros(size(accuracy,2),3);
colors(strcmp(population_names,'Proximity'),:) = cmap(1,:); % proximity = red
colors(strcmp(population_names,'Untuned'),:) = [0 0 0]; % untuned = dark gray
colors(strcmp(population_names,'Map'),:) = [1 1 1]; % map = white (yellow)
colors(strcmp(population_names,'Mixed'),:) = [1 0 1]; % mixed = magenta
colors(strcmp(population_names,'Suppressed'),:) = [0 0 1]; % suppressed = blue
colors(strcmp(population_names,'Shuffled'),:) = [0.5 0.5 0.5]; % shuffled = light gray

clear h1 h2 h3 h4 h5 h6
figure(2165); clf
figpos = [2 2 10 10];
set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',figpos)
hold on

% Conform color scheme with pie chart (red=prox, yellow=map, gray=untuned, white=shuffled)

% Get mean +/- std across subsample repetitions
err = std(y,[],3) / sqrt(size(y,3));
y = mean(y,3);
bar_width = 0.15;

% Plot 'proximity' population decoding
p = find(contains(population_names,'proximity', 'IgnoreCase',true));
for d = 1:num_dist
    x = d-bar_width*2.5;
    h1(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5);
    h1(d).FaceColor = 'flat';
    h1(d).FaceAlpha = 0.6;
    h1(d).CData = [colors(p,:)];
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'map' population decoding
p = find(contains(population_names,'map', 'IgnoreCase',true));
for d = 1:num_dist
    x = d-bar_width*0.5;
    h2(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5);
    h2(d).FaceColor = 'flat';
    h2(d).CData = [colors(p,:)];
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'suppressed' population decoding
p = find(contains(population_names,'suppressed', 'IgnoreCase',true));
for d = 1:num_dist
    x = d+bar_width*0.5;
    h3(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5);
    h3(d).FaceColor = 'flat';
    h3(d).FaceAlpha = 0.6;
    h3(d).CData = [colors(p,:)];
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'untuned' population decoding
p = find(contains(population_names,'untuned', 'IgnoreCase',true));
for d = 1:num_dist
    x = d+bar_width*1.5;
    h4(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5, 'FaceColor',[0.5 0.5 0.5]);
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'shuffled' population decoding
p = find(contains(population_names,'shuffle', 'IgnoreCase',true));
for d = 1:num_dist
    x = d+bar_width*2.5;
    h5(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5, 'FaceColor',[0.75 0.75 0.75]);
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'mixed' population decoding
p = find(contains(population_names,'mixed', 'IgnoreCase',true));
for d = 1:num_dist
    x = d-bar_width*1.5;
    h6(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5);
    h6(d).FaceColor = 'flat';
    h6(d).FaceAlpha = 0.6;
    h6(d).CData = [colors(p,:)];
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

legend([h1(1), h2(1), h6(1), h3(1), h4(1), h5(1)], population_names([2 3 4 5 1 6]), 'Location','northeast', 'Box','off')
set(gca, 'TickDir','out', 'Box','off', ...
    'XTick',1:numel(label_dists), 'XTickLabels',label_dists(sortOrder), 'XTickLabelRotation',0)
ylim([0 100])
ylabel('Accuracy (%)')
xlabel('Wall Distance (mm)')

fig_dir = 'E:\PrV_Wall_Recordings\Figures';
filename = fullfile(fig_dir, 'decoding_accuracy_proximity_map_bar');
saveas(gcf, filename, 'pdf')
exportgraphics(gcf, [filename '.png'], 'Resolution',600)

%% Panel C (right): Compute error-by-distance for summary period
y = cell2mat(cellfun(@(x) mean(x(:,summary_bins,:),2,'omitnan'), error_by_distance, 'Uni',0));

clear h1 h2 h3 h4 h5
figure(3498); clf
figpos = [2 2 10 10];
set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',figpos)
hold on

% Get mean +/- std across subsample repetitions
err = std(y,[],3) / sqrt(size(y,3));
y = mean(y,3);
bar_width = 0.15; % 0.1667;

% Plot 'proximity' population decoding
p = find(contains(population_names,'proximity', 'IgnoreCase',true));
for d = 1:num_dist
    x = d-bar_width*2.5; % d-bar_width*2;
    h1(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5);
    h1(d).FaceColor = 'flat';
    h1(d).FaceAlpha = 0.6;
    h1(d).CData = [colors(p,:)];
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'map' population decoding
p = find(contains(population_names,'map', 'IgnoreCase',true));
for d = 1:num_dist
    x = d-bar_width*0.5; % d-bar_width;
    h2(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5);
    h2(d).FaceColor = 'flat';
    h2(d).CData = [colors(p,:)];
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'suppressed' population decoding
p = find(contains(population_names,'suppressed', 'IgnoreCase',true));
for d = 1:num_dist
    x = d+bar_width*0.5;
    h3(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5);
    h3(d).FaceColor = 'flat';
    h3(d).FaceAlpha = 0.6;
    h3(d).CData = [colors(p,:)];
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'untuned' population decoding
p = find(contains(population_names,'untuned', 'IgnoreCase',true));
for d = 1:num_dist
    x = d+bar_width*1.5;
    h4(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5, 'FaceColor',[0.5 0.5 0.5]);
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'shuffled' population decoding
p = find(contains(population_names,'shuffle', 'IgnoreCase',true));
for d = 1:num_dist
    x = d+bar_width*2.5;
    h5(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5, 'FaceColor',[0.75 0.75 0.75]);
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

% Plot 'proximity' population decoding
p = find(contains(population_names,'mixed', 'IgnoreCase',true));
for d = 1:num_dist
    x = d-bar_width*1.5; 
    h6(d) = bar(x, y(d,p), bar_width, 'LineWidth',0.5);
    h6(d).FaceColor = 'flat';
    h6(d).FaceAlpha = 0.6;
    h6(d).CData = [colors(p,:)];
    errorbar(x, y(d,p), err(d,p), 'Color','k', 'LineWidth',0.5, 'CapSize',2)
end

legend([h1(1), h2(1), h6(1), h3(1), h4(1), h5(1)], population_names([2 3 4 5 1 6]), 'Location','northeast', 'Box','off')
set(gca, 'TickDir','out', 'Box','off', ...
    'XTick',1:numel(label_dists), 'XTickLabels',label_dists(sortOrder), 'XTickLabelRotation',0)
ylim([0 8.5])
ylabel('Classification Error (mm)')
xlabel('Wall Distance (mm)')

fig_dir = 'E:\PrV_Wall_Recordings\Figures';
filename = fullfile(fig_dir, 'decoding_error_proximity_map_bar');
saveas(gcf, filename, 'pdf')
exportgraphics(gcf, [filename '.png'], 'Resolution',600)

% Save accuracy source data (for stats and table)
distances = label_dists(sortOrder);
mean_accuracy = cellfun(@(x) squeeze(mean(x(:,summary_bins,:),2)), accuracy_by_distance, 'Uni',0);
save(fullfile(fig_dir, [rundate '_decoder_accuracy_source_data']), 'mean_accuracy', 'population_names', 'distances')

mean_error = cellfun(@(x) squeeze(mean(x(:,summary_bins,:),2)), error_by_distance, 'Uni',0);
save(fullfile(fig_dir, [rundate '_decoder_error_source_data']), 'mean_error', 'population_names', 'distances')

%% Run statistics on summary results (excluding suppression)
% clear; clc; close all

% Load source data ('mean_accuracy', 'population_names', 'distances')
fig_dir = 'E:\PrV_Wall_Recordings\Figures';
rundate = '20250827';

for load_accuracy = [true false]
    if load_accuracy
        filename = fullfile(fig_dir, [rundate '_decoder_accuracy_source_data']);
        data = matfile(filename).mean_accuracy;
        metrics = 'Accuracy';
        units = '(%)';
    else
        filename = fullfile(fig_dir, [rundate '_decoder_error_source_data']);
        data = matfile(filename).mean_error;
        metrics = 'Error';
        units = '(mm)';
    end

    disp(metrics)

    load(filename)
    all_population_names = population_names;
    
    % Compile 'y' and 'groups' for ANOVA
    [y,g1,g2] = deal([]);
    for p = 1:numel(population_names)
        % Get data
        yy = data{p}';
    
        % First level groups is the distance labels
        dd = repelem(distances,size(yy,1),1);
        g1 = [g1; dd(:)];
    
        % Second level groups is the population name (e.g.,'Proximity','Map','Shuffled')
        pp = repmat({population_names{p}},size(yy,1),size(yy,2));
        g2 = [g2; pp(:)];
    
        % Concatenate data
        y = [y; yy(:)];
    end
    
    % Run 2-way anova (X1 = distance, X2 = subpop)
    [pvalues, tbl, stats] = anovan(y, {g1, g2}, 'Display','off', 'Model','full');  

    % Run multiple comparisons correction
    [c12, m, h, gnames] = multcompare(stats, "Dimension",[1 2], "Display","off"); % default = Tukey's honestly significant difference procedure
    tbl12 = array2table(c12, "VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl12.("Group A") = gnames(tbl12.("Group A"));
    tbl12.("Group B") = gnames(tbl12.("Group B"));
        
    % Compare p-values for same distances across populations
    combos = nchoosek(1:numel(population_names),2);
    [p,md] = deal(zeros(numel(distances),size(combos,1)));
    for d = 1:numel(distances)
        gidx = find(contains(gnames,['X1=' num2str(distances(d)) ',']));
        for c = 1:size(combos,1)
            GroupA = gnames(strcmp(gnames,['X1=' num2str(distances(d)) ',X2=' num2str(population_names{combos(c,1)})]));
            GroupB = gnames(strcmp(gnames,['X1=' num2str(distances(d)) ',X2=' num2str(population_names{combos(c,2)})]));
            idx = find(strcmp(tbl12.("Group A"), GroupA) & strcmp(tbl12.("Group B"), GroupB));
            p(d,c) = tbl12.("P-value")(idx);
            md(d,c) = tbl12.("A-B")(idx);
        end
    end

    % Compare 'proximity' vs 'map' at each distance using Wilcoxon rank-sum test
    distances = sort(unique(g1),'ascend');
    num_dist = numel(distances);
    pvalues = nan(num_dist,1);
    prox = strcmpi(g2,'proximity');
    map = strcmpi(g2,'map');
    for d = 1:num_dist
        dd = g1 == distances(d);
        [pvalues(d), ~] = ranksum(y(dd & prox), y(dd & map));
    end

    % Compute FDR correction
    % [~, ~, ~, adj_p_values] = fdr_bh(p_values);
    pvalues_adj = mafdr(pvalues, 'BHFDR',true);
    disp('Adjusted p-values for proximity vs map at each distance:');
    disp(pvalues_adj);

    % Compare 'proximity' vs 'mixed' at each distance using Wilcoxon rank-sum test
    distances = sort(unique(g1),'ascend');
    num_dist = numel(distances);
    pvalues = nan(num_dist,1);
    prox = strcmpi(g2,'proximity');
    map = strcmpi(g2,'mixed');
    for d = 1:num_dist
        dd = g1 == distances(d);
        [pvalues(d), ~] = ranksum(y(dd & prox), y(dd & map));
    end

    % Compute FDR correction
    % [~, ~, ~, adj_p_values] = fdr_bh(p_values);
    pvalues_adj = mafdr(pvalues, 'BHFDR',true);
    disp('Adjusted p-values for proximity vs mixed at each distance:');
    disp(pvalues_adj);

    % Compare 'map' vs 'mixed' at each distance using Wilcoxon rank-sum test
    distances = sort(unique(g1),'ascend');
    num_dist = numel(distances);
    pvalues = nan(num_dist,1);
    prox = strcmpi(g2,'map');
    map = strcmpi(g2,'mixed');
    for d = 1:num_dist
        dd = g1 == distances(d);
        [pvalues(d), ~] = ranksum(y(dd & prox), y(dd & map));
    end

    % Compute FDR correction
    % [~, ~, ~, adj_p_values] = fdr_bh(p_values);
    pvalues_adj = mafdr(pvalues, 'BHFDR',true);
    disp('Adjusted p-values for map vs mixed at each distance:');
    disp(pvalues_adj);

    % Compare 'suppressed' vs 'shuffled' at each distance using Wilcoxon rank-sum test
    distances = sort(unique(g1),'ascend');
    num_dist = numel(distances);
    pvalues = nan(num_dist,1);
    prox = strcmpi(g2,'suppressed');
    map = strcmpi(g2,'shuffled');
    for d = 1:num_dist
        dd = g1 == distances(d);
        [pvalues(d), ~] = ranksum(y(dd & prox), y(dd & map));
    end

    % Compute FDR correction
    % [~, ~, ~, adj_p_values] = fdr_bh(p_values);
    pvalues_adj = mafdr(pvalues, 'BHFDR',true);
    disp('Adjusted p-values for suppressed vs. shuffled at each distance:');
    disp(pvalues_adj);
    
    % Get p-values for every possible comparison
    mean_diff = NaN(size(data{combos(c,1)},1),size(combos,1));
    for c = 1:size(combos,1)
        p1 = combos(c,1);
        p2 = combos(c,2);
        mean_diff(:,c) = mean(data{p1},2) - mean(data{p2},2);
    end
    
    % Compile table
    t = table;
    t.('Distance (mm)') = distances;
    for pp = 1:numel(population_names)
        pop = find(strcmp(all_population_names, population_names{pp}));
        t.([population_names{pp} ' ' metrics ' ' units]) = mean(data{pop},2); % mean
        t.([population_names{pp} ' 95% CI']) = 1.96 * std(data{pop},[],2) / sqrt(size(data{pop},2)); % CI
    end
    for c = 1:size(p,2)
        combo_name = [population_names{combos(c,1)} ' vs ' population_names{combos(c,2)} ' p-value'];
        t.(combo_name) = p(:,c);
    end
    
    disp(t)

    % Save table to CSV
    filename = fullfile(fig_dir, ['Decoder_Table_All_Combo_Stats_' metrics '.csv']);
    writetable(t, filename)
end


