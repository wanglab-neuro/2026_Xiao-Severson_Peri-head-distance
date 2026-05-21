clearvars; clc; close all;

compute_sdfs = true;
plot_sdfs = true;
cluster_cells = false;

%% Load data from Naive dataset 
data_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData';
psth{1} = matfile(fullfile(data_path, 'PSTH_naive.mat')).psth;
condition_names = {'naive'};
unit_ids = [];

%% Load data from Cut dataset 
data_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\One_Whisker_cut';
psth{1} = matfile(fullfile(data_path, 'PSTH_cut.mat')).psth;
psth{2} = matfile(fullfile(data_path, 'PSTH_cut_condition.mat')).psth;
condition_names = {'pre-cut','post-cut'};
unit_ids = {
    %'20240522_sp_session_1_unit_8'
    '20250327_sp_session_2_unit_21'

    % figure 2.2:
    %'20240626_sp_session_1_unit_13'
    %'20240624_sp_session_1_unit_2'
    % figure 2.3:
    %'20240926_sp_session_1_unit_11'
    %'20240626_sp_session_1_unit_6'
    %'20240926_sp_session_1_unit_9'
    % figure 2.4:
    %'20240703_sp_session_1_unit_12'
    %'20240703_sp_session_1_unit_5'
    %'20240815_sp_session_1_unit_16'
    % figure 2.5:
    %'20240905_sp_session_1_unit_4'
    %'20240905_sp_session_1_unit_9'
    %'20240924_sp_session_1_unit_8'

};

%% Load data from GtACR Laser dataset 

% %% Load data from TeLC dataset (Chapter 4)
% data_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\TeLC';
% psth = matfile(fullfile(data_path, 'PSTH_telc.mat')).psth;

data_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\GtACR';
psth{1} = matfile(fullfile(data_path, 'PSTH_laser.mat')).psth;
psth{2} = matfile(fullfile(data_path, 'PSTH_laser_condition.mat')).psth;
condition_names = {'no-laser','laser'};
unit_ids = {
    % Fig.S6E, top
    '20240307_sp_session_1_unit_22'

    % figure 3.2:
   % '20240312_sp_session_1_unit_5'
    %'20240311_sp_session_1_unit_18'
   % '20240312_sp_session_1_unit_3'
    %'20240311_sp_session_1_unit_21'
    %'20240308_sp_session_2_unit_1'
    %'20240307_sp_session_1_unit_9'
    % figure 3.3:
   % '20240307_sp_session_1_unit_2'
   % '20240311_sp_session_1_unit_14'
};
%% Load data from Tungsten dataset 
data_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\Tungsten';
psth{1} = matfile(fullfile(data_path, 'PSTH_tungsten.mat')).psth;
condition_names = {'Tungsten'};
unit_ids = {
    '20211015_t_br_session_2_unit_1'
    '20220706_t_br_session_2_unit_1'
};

%% Save Spike Density Function
cd(data_path)

[all_sdf, all_sem] = deal(cell(length(psth)));
for c = 1:length(psth)

    bin_width = psth{c}.bin_width;             % 0.2
    time_baseline = psth{c}.time_baseline;     % -4.4:bin_width:1.4
    time_pass = psth{c}.time_pass;             % 2.1:bin_width:11.9
    
    % Compute spike density function for each cell
    % * Use a gaussian kernel with 1ms width
    % * Use a bin_width time bin
    
    % Get spike times
    spike_times = psth{c}.dataTable.spike_times;
    
    % Get pass start times
    pass_starts = psth{c}.dataTable.pass_starts;
    
    % Get wall distance conditions
    conditions = psth{c}.dataTable.wall_distance;
    
    % Get number of cells
    num_cells = length(spike_times);
    
    trial_duration = time_pass(end) - time_baseline(1);
    
    keep_cells = false(size(psth{c}.dataTable,1),1);
    keep_cells(ismember(psth{c}.dataTable.Properties.RowNames, unit_ids)) = true;
    % keep_cells([1009]) = true;
    
    % Compute rasters
    % too big - requires 80GB memory
    % unitID = repelem(1:num_cells, cellfun(@length, spike_times));
    % all_rasters = EphysFun.MakeRasters(vertcat(spike_times{:}), unitID, 1);
    
    % Initialize variables
    timeUnit = 1; % 1ms
    time_res = timeUnit/1000; % 1 ms
    sigma = 200; % 200 = 5 ms
    causal = 0;
    format = 'spiketimes';
    
    % Define timeline
    time_post = 2.5; % time after end of wall pass
    % timeline = time_baseline(1):time_res:time_pass(end)-time_res+time_post;
    
    % Compute SDFs
    if compute_sdfs
        [all_sdf, all_sem, timeline] = wall_pass_compute_sdf( ...
            condition_names{c}, ...
            spike_times(keep_cells), ...
            pass_starts(keep_cells), ...
            conditions(keep_cells),...
            time_baseline, ...
            time_pass, ...
            time_res, ...
            sigma, ...
            causal, ...
            format, ...
            timeUnit, ...
            data_path);
    end
end

%% Plot SDFs

if strcmpi(condition_names{1}, 'tungsten')
    dist_centers = [];
else
    dist_centers = matfile('E:\Code\Analysis\NeuralDecoder\dist_centers.mat').dist_centers;
end

if plot_sdfs
    ymax = zeros(1,sum(keep_cells));
    for c = 1:length(psth)
        load(fullfile(data_path, 'processed_data', ['wp_resps_' condition_names{c} '.mat']));
        [all_sdf, all_sem] = deal(cell(1, length(wp_resps)));
        for cellNum = 1 : length(wp_resps)
            all_sdf{cellNum} = wp_resps{cellNum}.trials_sdf;
            all_sem{cellNum} = wp_resps{cellNum}.trials_sems;
        end
        timeline = wp_resps{1}.timeline;

        ymax = max([ymax; cellfun(@(sdf,sem) max(max(sdf+sem)), all_sdf,all_sem, 'Uni',1)]);
    end

    % Make y-limits the same across conditions
    ylimits = cat(1, zeros(size(ymax)), ymax)';

    % Call SDF plot function
    for c = 1:length(psth)
        load(fullfile(data_path, 'processed_data', ['wp_resps_' condition_names{c} '.mat']));
        [all_sdf, all_sem] = deal(cell(1, length(wp_resps)));
        for cellNum = 1 : length(wp_resps)
            all_sdf{cellNum} = wp_resps{cellNum}.trials_sdf;
            all_sem{cellNum} = wp_resps{cellNum}.trials_sems;
        end
        timeline = wp_resps{1}.timeline;

        wall_pass_plot_sdf( ...
            condition_names{c}, ...
            all_sdf, ...
            all_sem, ...
            timeline, ...
            conditions(keep_cells), ...
            dist_centers, ...
            data_path, ...
            psth{c}.dataTable(keep_cells,:).Properties.RowNames, ...
            keep_cells, ...
            ylimits);
    end
end


