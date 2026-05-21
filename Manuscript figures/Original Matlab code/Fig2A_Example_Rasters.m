% %% RASTER examplesfor Si Probe PrV Data
% Figure 2, panel A

%% Digitized raster by trials
clear; clc

% Manually select sessions and units
base_dir = 'F:\PrV_Wall_Recordings'; 
rec_dir = '20240905_sp\session_1';
session_path = fullfile(base_dir, rec_dir, 'session.mat');

% Raster settings
bin_width = 0.05; % time bin width in seconds
use_color = true; % false; %
use_shade = true; % false; %
save_plots = true;

% ephys_path = fullfile(base_dir, rec_dir, 'ephys');
% spike_sort_path = 'SpikeSorting\20230323_sp_ephys27_session2';
% mouseID = 'Ephys27';
% TTL_file = '20230323_sp_ephys27_session2_trialTTLs.csv';
% wall_dist_file = 'trialInfo_ephys27_3232023.csv';

% Load session master
master = SessionLoader2.SessionMaster();
for i = numel(master):-1:1
    master_id(i) = strcmp(master(i).folder, rec_dir);
end
master_id = find(master_id, 1, 'first');

% Load session
session = matfile(session_path).session;


% cd(ephys_path)

% Load uniform distances from decoder analysis
decoder_folder = 'E:\Code\Analysis\NeuralDecoder';
dist_centers = matfile(fullfile(decoder_folder, 'dist_centers.mat')).dist_centers;
dist_edges = matfile(fullfile(decoder_folder, 'dist_edges.mat')).dist_edges;

% Load data from session
spike_times = session.sortedSpikeTimes;
units = 1:numel(spike_times);
mouseID = session.mouseID;

for unit_num = 5 %:numel(units)

    u = units(unit_num);
    st = spike_times{u};
    
    % Get unique wall distances
    wall_distance = session.dataTable.wallDistanceMM;
    wall_distance_table = array2table(wall_distance, 'VariableNames',{'wall_distance'});
    
    TTLtable = table(session.pass_start_times, session.pass_end_times);
    
    % Only use pre-cut condition
    if any(strcmp(session.dataTable.Properties.VariableNames,'trialConditions'))
        trial_conditions = cellfun(@(x) isempty(x), session.dataTable.trialConditions, 'Uni',1);
    else
        trial_conditions = true(size(TTLtable,1));
    end
    
    % Remove excluded Zaber loop trials (in SessionMaster)
    zloop = master(master_id).zaber_loops;
    trials_included = trial_conditions;
    distances = unique(wall_distance);
    n = numel(distances);

    % Include trials with distance in uniform set using nearest neighbors
    idx = false(size(trials_included));
    for d = 1:numel(dist_centers)
        [~, nearest_neighbor] = min(abs(distances - dist_centers(d)));
        nearest_distance = distances(nearest_neighbor);
        idx(wall_distance == nearest_distance) = true;
    end
    trials_included(~idx) = false;

    % Hack for this session
    if strcmp(session.mouseID,'ephys106') && strcmp(session.recID,'20240924_sp_session_1')
        trials_included(181:end) = false;
    end

    % Include trials in TTL and wall distance tables
    TTLtable = TTLtable(trials_included,:);
    wall_distance_table = wall_distance_table(trials_included,:);
    
    % Chunk the time series data
    post_time = [0 4]; % seconds after "endTime"
    numTrials = size(TTLtable, 1);
    st_by_trials = cell(numTrials, 1);
    for i = 1:numTrials
        startTime = TTLtable.Var1(i) + post_time(1);
        endTime = TTLtable.Var2(i) + post_time(2);
        st_by_trials{i} = st(st >= startTime & st <= endTime) - startTime;
    end
    
    % Compile spike times into table
    st_by_trials_table = table(st_by_trials, 'VariableNames', {'st_by_trials'});
    
    % Assert TTLtable and wall_distance_table are same length
    num_TTLs = size(TTLtable,1);
    num_wall_pass = size(wall_distance_table,1);
    if num_TTLs > num_wall_pass
        TTLtable = TTLtable(2:num_TTLs,:);
        st_by_trials_table = st_by_trials_table(2:num_TTLs,:);
    elseif num_wall_pass > num_TTLs
        wall_distance_table = wall_distance_table(1:num_TTLs,:);
    end
    
    % Concatenate "TTLtable" and "wall_distance_table" horizontally
    TTLtable = [TTLtable, wall_distance_table];
    
    % Add spike times into the table
    TTLtable = [TTLtable, st_by_trials_table];
    
    % Sort the table based on the values in the "wall_distance" column
    TTLtable = sortrows(TTLtable, 'wall_distance');
    
    % Convert spike times into spike counts matrix
    pass_time = (TTLtable.Var2 + post_time(2)) - (TTLtable.Var1 + post_time(1));
    pass_time_end = floor(median(pass_time));
    bin_edges = 0 : bin_width : pass_time_end;
    bin_centers = mean([bin_edges(1:end-1); bin_edges(2:end)]);
    spike_counts = cell2mat(shiftdim(cellfun(@(x) histcounts(x, bin_edges), TTLtable.st_by_trials, 'Uni',0)));
    
    % Get unique wall distances
    unique_wall_distances = unique(TTLtable.wall_distance);
    unique_spike_counts = unique(spike_counts);
    
    % Compile color map for raster (wall distance color if bin contains spike, otherwise white)
    cmap = flipud(turbo(numel(unique_wall_distances)));
    shademap = flipud(ones(numel(unique_spike_counts),3) .* linspace(0,1,numel(unique_spike_counts))');
    colors = ones(size(spike_counts,1),size(spike_counts,2),3);
    
    for i = 1:size(colors,1)
        c = cmap(unique_wall_distances==TTLtable.wall_distance(i),:);
        inds = spike_counts(i,:) > 0;
        for ind = 1:numel(inds)
            if spike_counts(i,ind) > 0
                shade = shademap(unique_spike_counts == spike_counts(i,ind),:);
                if use_color
                    if use_shade
                        colors(i,ind,:) = c .* shade;
                    else
                        colors(i,ind,:) = c;
                    end
                else
                    if use_shade
                        colors(i,ind,:) = shade;
                    else
                        colors(i,ind,:) = [0,0,0];
                    end
                end
            end
        end
    end
    
    X = (0:size(colors,2)-1) / (1/bin_width) + post_time(1);
    Y = 1:size(colors,1);
    C = colors;
    
    clear dist_trials
    for d = numel(dist_centers):-1:1
        dist_trials(d) = sum(wall_distance_table.wall_distance == distances(d));
    end
    dist_trials = cumsum(dist_trials);
    
    % Plot raster as image
    figpos = [50 50 800 800];
    h1 = figure(1); clf
    set(gcf, ...
        'Color','w', ...
        'renderer','painters', ...
        'Units','pixels', ...
        'Position',figpos)
    imagesc(X,Y,C)
    flipped_cmap = flipud(cmap);
    colormap(flipped_cmap);
    cb = colorbar;
    set(cb, 'Ticks',[], 'TickLabels',[], 'TickDir','out')
    caxis([min(unique_wall_distances) max(unique_wall_distances)]);
    num_ticks = pass_time_end / 2; % Ticks at 2-second intervals
    % xticklabel = round(bin_centers(xt));
    axis square
    
    set(gca, 'TickDir','out', 'box','off', 'YTick',dist_trials)
    ylabel(cb, 'Trial number')
    xlabel('Time (s)')
    title(gca,strjoin({mouseID, rec_dir, 'Unit', num2str(unit_num)}, ' '))

    if save_plots
        fig_dir = 'E:\PrV_Wall_Recordings\Figures\Rasters';
        if ~isfolder(fig_dir); mkdir(fig_dir); end
        filename = fullfile(fig_dir, strjoin({mouseID, strrep(rec_dir,filesep,'_'), 'Unit', num2str(unit_num)}, ' '));
        savefig(gcf,filename)
        saveas(gcf,filename,'png')
        saveas(gcf,filename,'pdf')
    end

end

%close all