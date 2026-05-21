%% Figure S1: Whisker retraction
% Whisker-Wall Tracking
clear; clc
settings = struct;
settings.fig_dir = 'C:\Users\kssev\Desktop\Manuscript\Figure S1 - whisker retraction';
settings.data_dir = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\Whisker-Wall Tracking';

% Binning settings
settings.dist_centers = 7:2:23;
settings.dist_edges = 6:2:24;

% Other settings
settings.frame_rate = 500; % fps
settings.baseline_period = [-4.5 0];
settings.wall_period = [3 11];
settings.psth_period = [0 12]; 

%% Panel A: Example video frame, angle measurement schematic
% beh1

fn = 'E:\PrV_Wall_Recordings\20231212_sp\session_1\video\280000_282499.mp4';
Martiny('tifPath',fn)


%% Panel B: Example angle trace, midpoint measurement

% Get data from example mouse (beh 103)
example_dir = 'E:\PrV_Wall_Recordings\20231212_sp\session_1';
input_file = 'whiskerTracking_session.mat';

session = matfile(fullfile(example_dir,input_file)).session;

frames = 17801:18400;
w = 2; % C1 whisker of interest;
c = 1; % condition (1=awake, 2=anesth)

% Get time (x) vs. whisker angle (y1) and whisker midpoint (y2)
x = (frames - 1) / settings.frame_rate; % time in sec
theta_raw = session.whiskerData.theta(frames,w);

% Get midpoint
midpoint = get_whisker_data(session, 'midpoint');
midpoint = midpoint(frames,w);
                
% Plot trace with whisking during wall pass zoomed in
figure(1000); clf
figpos = [5 5 15 5];
set(gcf, 'Color','w', 'Units','centimeters', 'Renderer','painters', 'Position',figpos)
hold on
plot(x, midpoint, '-', 'LineWidth',2, 'Color',[0.5 0.5 0.5]) % midpoint angle
plot(x, theta_raw, '-', 'LineWidth',2, 'Color','k') % theta angle

% Format axes
set(gca, 'TickDir','out', 'Box','off')
xlabel('Time (s)')
ylabel('Angle (deg)')
pause(0.001)

% Save figure
filename = 'FigS1_PanelB_whisker_theta_midpoint_trace';
saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf')
saveas(gcf, fullfile(settings.fig_dir, filename), 'png')

%% Panel C: Example traces aligned to wall pass

clearvars -except settings

% Get data from example mouse (beh 104)
% example_dir = 'E:\PrV_Wall_Recordings\20231204_sp\session_1';
% input_file = 'whiskerTracking_session.mat';
% example_mouse = 'beh104';

example_dir = 'E:\PrV_Wall_Recordings\20231213_sp\session_1';
input_file = 'whiskerTracking_session.mat';
example_mouse = 'beh102';

% Load session data
load(fullfile(example_dir, input_file));

% Fix wrong condition label
session.dataTable.trialConditions{201} = 'anesth';

% Settings
example_distance = 9;
whisker_id = 2; % C1
baseline_period = settings.baseline_period;
wall_period = settings.wall_period;
psth_period = settings.psth_period;

% Get pass indices for this distance
wd = unique(session.dataTable.wallDistanceMM);
[~,wdex] = min(abs(wd-example_distance));
wdex = wd(wdex);

% Process whisker data
var = get_whisker_data(session, 'kappa');
var = smoothdata(var, 1, 'sgolay',51, 'degree',3);

% var = get_whisker_data(session, 'midpoint');
% var = smoothdata(var, 1, 'sgolay',51, 'degree',3);

% Synchronize ephys with video frame times
frame_times = session.frame_times;

% Get passes from awake
% passIdx = session.dataTable.wallDistanceMM==wdex & cellfun(@(x) strcmp(x,'anesth'), session.dataTable.trialConditions, 'Uni',1);
passIdx = session.dataTable.wallDistanceMM==wdex & cellfun(@(x) isempty(x), session.dataTable.trialConditions, 'Uni',1);
pass_start_times = session.pass_start_times(passIdx);

[baseline, wall, psth] = deal(cell(numel(pass_start_times),1));
for i = 1:numel(pass_start_times)
    % Get baseline period
    idx = frame_times > pass_start_times(i)+baseline_period(1) & frame_times < pass_start_times(i)+baseline_period(2);
    baseline{i} = var(idx, whisker_id);

    % Get wall period data
    idx = frame_times > pass_start_times(i)+wall_period(1) & frame_times < pass_start_times(i)+wall_period(2);
    wall{i} = var(idx, whisker_id);

    % Get full psth period
    idx = find(frame_times > (pass_start_times(i)+psth_period(1)) & frame_times < (pass_start_times(i)+psth_period(2)));
    psth{i} = var(idx, whisker_id);
end

% Assemble table
T = table(pass_start_times, baseline, wall, psth);
mean_baseline = mean(cellfun(@(x) mean(x, 'omitnan'), T.baseline, 'Uni',1));

% Subtract mean baseline to get delta-kappa
T.baseline = cellfun(@(x) x - mean_baseline, T.baseline, 'Uni',0);
T.psth = cellfun(@(x) x - mean_baseline, T.psth, 'Uni',0);

% Plot
figure(2); clf; hold on
set(gcf, 'Color','w', 'Units','centimeters', 'Renderer','painters', 'Position',[2 2 8 4])
for i = 1:size(T,1)
    x = (0:(numel(T.psth{i})-1)) / settings.frame_rate;
    y = T.psth{i};
    plot(x, y, 'Color',[0 0 0 0.5], 'LineWidth',1)
end

% Format axes
title([example_mouse ' C1 whisker, distance = ' num2str(example_distance) ' mm'], 'FontWeight','normal', 'FontSize',6)
xlabel('Time (s)')
ylabel('\Delta\kappa (mm^{-1})')
set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)

% Save figure
filename = 'FigS1_PanelC_wallpass_curvature_traces';
saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf')
saveas(gcf, fullfile(settings.fig_dir, filename), 'png')

%% Panel D: Awake vs. anesthetized average PSTH for each distance

% Get data from example mouse (beh 102)
example_dir = 'E:\PrV_Wall_Recordings\20231213_sp\session_1';
input_file = 'whiskerTracking_session.mat';
example_mouse = 'beh102';
whisker_id = 2; % C1
baseline_period = settings.baseline_period;
wall_period = settings.wall_period;
psth_period = settings.psth_period;
settings.varname = 'midpoint';

colors = flipud(turbo(numel(settings.dist_centers)));
distance_names = arrayfun(@(x) [num2str(x) ' mm'], settings.dist_centers, 'Uni',0);

% Load session data
load(fullfile(example_dir, input_file));

% Fix wrong condition label
session.dataTable.trialConditions{201} = 'anesth';

% Process whisker data
var = get_whisker_data(session, settings.varname);

% Synchronize ephys with video frame times
frame_times = session.frame_times;
median_frame_duration = mean(diff(frame_times));
frame_times = frame_times(1) : median_frame_duration : (frame_times(1) + (size(var,1)-1)*median_frame_duration);

% Get pass indices for each distance, awake (1) vs. anesth (2), Cut to even size
baseline_length = diff(baseline_period) * settings.frame_rate;
wall_length = diff(wall_period) * settings.frame_rate;
psth_length = diff(psth_period) * settings.frame_rate;

% Downsample psth
ds = 10; % Set downsample factor (binning)
edges = 0 : ds : psth_length; % Define bin edges for equally spaced intervals
time_edges = edges / settings.frame_rate;
[~, ~, bin] = histcounts(1:psth_length, edges); % Get the bin index for each element
bin_centers = mean(cat(1, time_edges(1:end-1), time_edges(2:end)));

% Save psths as matrix (time, distance, condition)
clear psth_mean psth_sem
[psth_mean, psth_sem] = deal(NaN(psth_length / ds, numel(settings.dist_centers), 2));
baseline_mean = NaN(baseline_length, numel(settings.dist_centers), 2);
for c = 1:2
    for d = 1:numel(settings.dist_centers)
        wd = unique(session.dataTable.wallDistanceMM);
        [~,wdex] = min(abs(wd-settings.dist_centers(d)));
        wdex = wd(wdex);

        % Get passes from awake
        if c==1
            passIdx = session.dataTable.wallDistanceMM==wdex & cellfun(@(x) isempty(x), session.dataTable.trialConditions, 'Uni',1);
        else
            passIdx = session.dataTable.wallDistanceMM==wdex & cellfun(@(x) strcmp(x,'anesth'), session.dataTable.trialConditions, 'Uni',1);
        end
        pass_start_times = session.pass_start_times(passIdx);

        [baseline, wall, psth] = deal(cell(numel(pass_start_times),1));
        for i = 1:numel(pass_start_times)
            % Get baseline period
            idx = frame_times > pass_start_times(i)+baseline_period(1) & frame_times < pass_start_times(i)+baseline_period(2);
            baseline{i} = var(idx, whisker_id);
        
            % Get wall period data
            idx = frame_times > pass_start_times(i)+wall_period(1) & frame_times < pass_start_times(i)+wall_period(2);
            wall{i} = var(idx, whisker_id);
        
            % Get full psth period
            idx = find(frame_times > (pass_start_times(i)+psth_period(1)) & frame_times < (pass_start_times(i)+psth_period(2)));
            psth{i} = var(idx, whisker_id);
        end

        % Interpolate if necessary
        for i = 1:numel(pass_start_times)
            baseline{i} = interp1(1:numel(baseline{i}), baseline{i}, 1:baseline_length, 'linear', NaN);
            wall{i} = interp1(1:numel(wall{i}), wall{i}, 1:wall_length, 'linear', NaN);
            psth{i} = interp1(1:numel(psth{i}), psth{i}, 1:psth_length, 'linear', NaN);
        end
        
        % Downsample temporal dimension
        tmp_wall = cell2mat(psth);
        tmp_wall = reshape(tmp_wall, size(tmp_wall,1), ds, (size(tmp_wall,2)/ds));
        tmp_wall = squeeze(mean(tmp_wall,2, 'omitnan'));
        
        % Insert mean into matrix
        psth_mean(:,d,c) = mean(tmp_wall,1,'omitnan');
        psth_sem(:,d,c) = std(tmp_wall,1,'omitnan') / sqrt(numel(pass_start_times));

        % Insert mean into matrix
        baseline_mean(:,d,c) = mean(cell2mat(baseline),1,'omitnan');
    end
end

% Get mean (baseline-subtracted) psth for each distance
mean_baseline = squeeze(mean(mean(baseline_mean,2,'omitnan'),1,'omitnan'));
mean_psth_awake = squeeze(psth_mean(:,:,1)) - mean_baseline(1);
mean_psth_anesth = squeeze(psth_mean(:,:,2)) - mean_baseline(2);

%% Panel D (left): Plot awake mean PSTH +/- SEM

warning('off')
yl = [-22 17];
xl = psth_period;

figure(1); clf; hold on
set(gcf, 'Color','w', 'Units','centimeters', 'Renderer','painters', 'Position',[2 2 9 6])
clear h
x = shiftdim(psth_period(1) : (1/settings.frame_rate * ds) : (psth_period(2) - 1/settings.frame_rate * ds));
for d = 1:numel(settings.dist_centers)
    y = squeeze(mean_psth_awake(:,d));
    err = squeeze(psth_sem(:,d,1));
    error_shade(y, err, x, 'Color',colors(d,:))
    h(d) = plot(x, y, 'LineWidth',1, 'Color',colors(d,:));
end

% Insert legend
legend(h, distance_names, 'Location','northeastoutside', 'Box','off')

% Format axes
axis square
xlim(xl); ylim(yl)
title([example_mouse ' awake, C1 whisker, distance = ' ' mm'], 'FontWeight','normal', 'FontSize',6)
xlabel('Time (s)')
ylabel('Midpoint (deg)')
set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)

% Save figure
pause(0.001)
filename = ['FigS1_PanelD_awake_wallpass_average_psth_by_distance_' example_mouse];
saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf')
saveas(gcf, fullfile(settings.fig_dir, filename), 'png')
pause(0.001)

% Panel D (right): Plot anesthetized mean PSTH +/- SEM
figure(2); clf; hold on
set(gcf, 'Color','w', 'Units','centimeters', 'Renderer','painters', 'Position',[2 2 9 6])
clear h
x = shiftdim(psth_period(1) : (1/settings.frame_rate * ds) : (psth_period(2) - 1/settings.frame_rate * ds));
for d = 1:numel(settings.dist_centers)
    y = squeeze(mean_psth_anesth(:,d));
    err = squeeze(psth_sem(:,d,2));
    error_shade(y, err, x, 'Color',colors(d,:))
    h(d) = plot(x, y, 'LineWidth',1, 'Color',colors(d,:));
end

% Insert legend
legend(h, distance_names, 'Location','northeastoutside', 'Box','off')

% Format axes
axis square
xlim(xl); ylim(yl)
title([example_mouse ' anesth, C1 whisker, distance = ' ' mm'], 'FontWeight','normal', 'FontSize',6)
xlabel('Time (s)')
ylabel(['\Delta' '\theta_{mid} (deg)'])
set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)

% Save figure
pause(0.001)
filename = ['FigS1_PanelE_anesth_wallpass_average_psth_by_distance_' example_mouse];
saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf')
saveas(gcf, fullfile(settings.fig_dir, filename), 'png')
pause(0.001)


%% Panel E: Display mean displacement by distance for one example mouse
% Take mean over time first, show error and do stats?

warning('off')
settings.varname = 'midpoint';

% Settings
settings.varname = 'midpoint';
baseline_period = settings.baseline_period;
wall_period = settings.wall_period;
psth_period = settings.psth_period;
colors = flipud(turbo(numel(settings.dist_centers)));
distance_names = arrayfun(@(x) [num2str(x) ' mm'], settings.dist_centers, 'Uni',0);

base_dir = 'E:\PrV_Wall_Recordings';
example_dirs = {'20231213_sp\session_1';};
subjects = {'beh102'};
input_file = 'whiskerTracking_session.mat';

% Get data from all subjects
[mean_var_awake, mean_var_anesth, mean_baseline] = deal(cell(numel(subjects),1));
for i = 1:numel(subjects)
    % Load session data
    example_dir = fullfile(base_dir, example_dirs{i});
    filename = fullfile(example_dir, input_file);
    disp(['(' num2str(i) '/' num2str(numel(subjects)) '): Loading whisker data from ' subjects{i}])
    load(filename); % 'session'
    [mean_var_awake{i}, mean_var_anesth{i}, mean_baseline{i}] = get_subject_data(session, settings);
end

% TODO: Do stats

%% Panel E: Plot for each whisker Greek, C1, C2
% Get data as mean +/- SEM change in midpoint from baseline across repetitions
x = shiftdim(settings.dist_centers);
y_awake = shiftdim(mean(mean_var_awake{1},1));
err_awake = shiftdim(std(mean_var_awake{1},[],1) / sqrt(size(mean_var_awake{1},1)));
y_anesth = shiftdim(mean(mean_var_anesth{1},1));
err_anesth = shiftdim(std(mean_var_anesth{1},[],1) / sqrt(size(mean_var_anesth{1},1)));

% Plot figure for whisker
whisker_ids = {'Greek','C1','C2'};
for wid = 2
    figure(2+wid); clf; hold on
    set(gcf, 'Color','w', 'Units','centimeters', 'Position',[2 2 8 4])
    error_shade(y_awake(:,wid), err_awake(:,wid), x, 'Color','g');
    error_shade(y_anesth(:,wid), err_anesth(:,wid), x, 'Color','b');
    h(1) = plot(x, y_awake(:,wid), 'g'); % awake
    h(2) = plot(x, y_anesth(:,wid), 'b'); % anesth
    legend(h, {'Awake','Anesthetized'}, 'Location','northeastoutside', 'Box','off')
    
    % Format axes
    xlim([settings.dist_centers(1) settings.dist_centers(end)])
    ylim([-25 5])
    axis square
    title([subjects{1} ' ' whisker_ids{wid} ' mean midpoint change'], 'FontWeight','normal', 'FontSize',6)
    xlabel('Wall Distance (mm)')
    ylabel(['\Delta' '\theta_{mid} (deg)'])
    set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)
    
    % Save figure
    pause(0.001)
    filename = ['FigS1_PanelF_mean_midpoint_change_by_distance_' whisker_ids{wid} '_' subjects{1}];
    saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf');
    saveas(gcf, fullfile(settings.fig_dir, filename), 'png');
    pause(0.001)
end

%% Panel F: Plot average across mice for each whisker

clearvars -except settings
warning('off')

% Settings
settings.varname = 'midpoint';
baseline_period = settings.baseline_period;
wall_period = settings.wall_period;
psth_period = settings.psth_period;
colors = flipud(turbo(numel(settings.dist_centers)));
distance_names = arrayfun(@(x) [num2str(x) ' mm'], settings.dist_centers, 'Uni',0);

base_dir = 'E:\PrV_Wall_Recordings';
example_dirs = {
    '20231204_sp\session_1';
    '20231212_sp\session_1';
    '20231212_sp\session_3';
    '20231213_sp\session_1';
    };
subjects = {'beh104', 'beh103', 'beh101', 'beh102'};
input_file = 'whiskerTracking_session.mat';

% Get data from all subjects
[mean_var_awake, mean_var_anesth, mean_baseline_awake, mean_baseline_anesth] = deal(cell(numel(subjects),1));
for i = 1:numel(subjects)
    % Load session data
    example_dir = fullfile(base_dir, example_dirs{i});
    filename = fullfile(example_dir, input_file);
    disp(['(' num2str(i) '/' num2str(numel(subjects)) '): Loading whisker data from ' subjects{i}])
    load(filename); % 'session'
    [mean_var_awake{i}, mean_var_anesth{i}, mean_baseline_awake{i}, mean_baseline_anesth{i}] = get_subject_data(session, settings);
end

%% Assemble source data table
T = build_whisker_source_table(subjects, mean_var_awake, mean_var_anesth, distance_names, settings.data_dir, settings.varname);

% Save source data as CSV
unique_subjects = categories(T.subject);
for i = 1:numel(unique_subjects)
    subj = char(unique_subjects(i));
    Tsubj = T(T.subject == subj, :);
    out_csv = fullfile(settings.data_dir, ...
        sprintf('whisker_retraction_source_%s_%s.csv', settings.varname, subj));
    writetable(Tsubj, out_csv);
end

%% *** (TODO) Do stats on each subject, compare effects across all subjects
T = [];
for i = 1:numel(subjects)
    t = readtable(fullfile(settings.data_dir, sprintf('whisker_retraction_source_midpoint_%s.csv', subjects{i})));
    T = cat(T, t);
    clear t
end

%% Panel XX: Plot mean whisker midpoint change, averaged across subjects

% Load source data

% Format data from table

% Concatenate data across mice
num_whiskers = 3;
[y_awake_all, y_anesth_all] = deal(zeros(numel(subjects), numel(settings.dist_centers), num_whiskers));
for i = 1:numel(subjects)
    % Get data as mean +/- SEM change in midpoint from baseline across repetitions
    y_awake_all(i,:,:) = squeeze(mean(mean_var_awake{i},1));
    y_anesth_all(i,:,:) = squeeze(mean(mean_var_anesth{i},1));
end

% Get mean +/- sem across subjects
y_awake = squeeze(mean(y_awake_all,1));
err_awake = squeeze(std(y_awake_all,[],1) / sqrt(size(y_awake_all,1)));
y_anesth = squeeze(mean(y_anesth_all,1));
err_anesth = squeeze(std(y_anesth_all,[],1) / sqrt(size(y_anesth_all,1)));
x = shiftdim(settings.dist_centers);
   
whisker_names = {'Greek Arc','1st Arc','2nd Arc'};
for wid = 1:num_whiskers
    % Plot figure for each whisker
    figure(100+wid); clf; hold on
    set(gcf, 'Color','w', 'Units','centimeters', 'Position',[2 2 8 4])
    error_shade(y_awake(:,wid), err_awake(:,wid), x, 'Color','g');
    error_shade(y_anesth(:,wid), err_anesth(:,wid), x, 'Color','b');
    h(1) = plot(x, y_awake(:,wid), 'g'); % awake
    h(2) = plot(x, y_anesth(:,wid), 'b'); % anesth
    legend(h, {'Awake','Anesthetized'}, 'Location','northeastoutside', 'Box','off')
    
    % Format axes
    xlim([settings.dist_centers(1) settings.dist_centers(end)])
    ylim([-25 5])
    axis square
    title([whisker_names{wid} ' mean midpoint change'], 'FontWeight','normal', 'FontSize',6)
    xlabel('Wall Distance (mm)')
    ylabel(['\Delta' '\theta_{mid} (deg)'])
    set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)
    set(gca, 'XTick',settings.dist_centers, 'XTickLabelRotation',0)

    % Save cross-subject comparison figure
    pause(0.001)
    filename = ['FigS1_PanelG_grand_mean_midpoint_change_by_distance-' whisker_names{wid}];
    saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf');
    saveas(gcf, fullfile(settings.fig_dir, filename), 'png');
    pause(0.001)
end

%% Panel G: Awake vs anesthetized average PSTH for each distance (delta-kappa)

% Use same example mouse as Panel D
example_dir = 'E:\PrV_Wall_Recordings\20231213_sp\session_1';
input_file = 'whiskerTracking_session.mat';
example_mouse = 'beh102';
whisker_id = 2; % C1

baseline_period = settings.baseline_period;
wall_period = settings.wall_period;
psth_period = settings.psth_period;

colors = flipud(turbo(numel(settings.dist_centers)));
distance_names = arrayfun(@(x) [num2str(x) ' mm'], settings.dist_centers, 'Uni',0);

% Load session data
load(fullfile(example_dir, input_file));
session.dataTable.trialConditions{201} = 'anesth';

% Get curvature trace
settings.varname = 'kappa';
var = get_whisker_data(session, settings.varname);

% Sync frame times
frame_times = session.frame_times;
median_frame_duration = mean(diff(frame_times));
frame_times = frame_times(1) : median_frame_duration : (frame_times(1) + (size(var,1)-1)*median_frame_duration);

% Lengths
baseline_length = diff(baseline_period) * settings.frame_rate;
wall_length = diff(wall_period) * settings.frame_rate;
psth_length = diff(psth_period) * settings.frame_rate;

% Save psths as matrix (time, distance, condition)
[psth_mean, psth_sem] = deal(NaN(psth_length, numel(settings.dist_centers), 2));

for c = 1:2
    for d = 1:numel(settings.dist_centers)
        wd = unique(session.dataTable.wallDistanceMM);
        [~,wdex] = min(abs(wd-settings.dist_centers(d)));
        wdex = wd(wdex);

        if c==1
            passIdx = session.dataTable.wallDistanceMM==wdex & cellfun(@(x) isempty(x), session.dataTable.trialConditions, 'Uni',1);
        else
            passIdx = session.dataTable.wallDistanceMM==wdex & cellfun(@(x) strcmp(x,'anesth'), session.dataTable.trialConditions, 'Uni',1);
        end
        pass_start_times = session.pass_start_times(passIdx);
        
        [k0] = deal(NaN(numel(pass_start_times),1));
        for i = 1:numel(pass_start_times)
            % baseline
            idx = frame_times >= pass_start_times(i)+baseline_period(1) & frame_times < pass_start_times(i)+baseline_period(2);
            k0(i) = mean(var(idx, whisker_id), 'omitnan');
        end

        [baseline, psth] = deal(cell(numel(pass_start_times),1));
        for i = 1:numel(pass_start_times)
            % baseline
            idx = frame_times >= pass_start_times(i)+baseline_period(1) & frame_times < pass_start_times(i)+baseline_period(2);
            baseline{i} = var(idx, whisker_id);

            % psth
            idx = frame_times >= pass_start_times(i)+psth_period(1) & frame_times < pass_start_times(i)+psth_period(2);
            psth{i} = var(idx, whisker_id);

            baseline{i} = baseline{i} - mean(k0);
            psth{i} = psth{i} - mean(k0);
        end

        % Interpolate to fixed length
        for i = 1:numel(pass_start_times)
            baseline{i} = interp1(1:numel(baseline{i}), baseline{i}, 1:baseline_length, 'linear', NaN);
            psth{i}     = interp1(1:numel(psth{i}),     psth{i},     1:psth_length,     'linear', NaN);
        end

        tmp_psth = cell2mat(psth); % [nPass x psth_length]
        psth_mean(:,d,c) = mean(tmp_psth,1,'omitnan')';
        psth_sem(:,d,c)  = std(tmp_psth,0,1,'omitnan')' / sqrt(size(tmp_psth,1));
    end
end

% Time vector
x = shiftdim(psth_period(1) : (1/settings.frame_rate) : (psth_period(2) - 1/settings.frame_rate));

yl = [-1.5e-3 2e-4];

% Plot awake
figure(501); clf; hold on
set(gcf,'Color','w','Units','centimeters','Renderer','painters','Position',[2 2 9 6])
clear h
for d = 1:numel(settings.dist_centers)
    y   = squeeze(psth_mean(:,d,1));
    err = squeeze(psth_sem(:,d,1));
    error_shade(y, err, x, 'Color',colors(d,:));
    h(d) = plot(x, y, 'LineWidth',1, 'Color',colors(d,:));
end
legend(h, distance_names, 'Location','northeastoutside', 'Box','off')
axis square
xlim(psth_period)
ylim(yl)
xlabel('Time (s)')
ylabel('\Delta\kappa (mm^{-1}')
title([example_mouse ' awake, C1 whisker, \Delta\kappa PSTH'], 'FontWeight','normal', 'FontSize',6)
set(gca,'TickDir','out','Box','off','FontSize',6)

filename = ['FigS1_PanelG_awake_wallpass_deltaKappa_psth_by_distance_' example_mouse];
saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf')
saveas(gcf, fullfile(settings.fig_dir, filename), 'png')

% Plot anesth
figure(502); clf; hold on
set(gcf,'Color','w','Units','centimeters','Renderer','painters','Position',[2 2 9 6])
clear h
for d = 1:numel(settings.dist_centers)
    y   = squeeze(psth_mean(:,d,2));
    err = squeeze(psth_sem(:,d,2));
    error_shade(y, err, x, 'Color',colors(d,:));
    h(d) = plot(x, y, 'LineWidth',1, 'Color',colors(d,:));
end
legend(h, distance_names, 'Location','northeastoutside', 'Box','off')
axis square
xlim(psth_period)
ylim(yl)
xlabel('Time (s)')
ylabel('\Delta\kappa (mm^{-1})')
title([example_mouse ' anesth, C1 whisker, \Delta\kappa PSTH'], 'FontWeight','normal', 'FontSize',6)
set(gca,'TickDir','out','Box','off','FontSize',6)

filename = ['FigS1_PanelG_anesth_wallpass_deltaKappa_psth_by_distance_' example_mouse];
saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf')
saveas(gcf, fullfile(settings.fig_dir, filename), 'png')

%% Panel H: Mean delta-kappa by distance for one example mouse (like Panel E midpoint)

warning('off')

settings.varname = 'kappa';

base_dir = 'E:\PrV_Wall_Recordings';
example_dirs = {'20231213_sp\session_1';};
subjects = {'beh102'};
input_file = 'whiskerTracking_session.mat';

[mean_var_awake, mean_var_anesth] = deal(cell(numel(subjects),1));
for i = 1:numel(subjects)
    example_dir = fullfile(base_dir, example_dirs{i});
    disp(['(' num2str(i) '/' num2str(numel(subjects)) '): Loading whisker data from ' subjects{i}])
    load(fullfile(example_dir, input_file)); % 'session'
    [mean_var_awake{i}, mean_var_anesth{i}] = get_subject_data(session, settings);
end

% Plot for C1 only (wid=2)
x = shiftdim(settings.dist_centers);
y_awake  = shiftdim(mean(mean_var_awake{1},1));
err_awake = shiftdim(std(mean_var_awake{1},[],1) / sqrt(size(mean_var_awake{1},1)));
y_anesth  = shiftdim(mean(mean_var_anesth{1},1));
err_anesth = shiftdim(std(mean_var_anesth{1},[],1) / sqrt(size(mean_var_anesth{1},1)));

whisker_ids = {'Greek','C1','C2'};
wid = 2;

figure(601); clf; hold on
set(gcf, 'Color','w', 'Units','centimeters', 'Position',[2 2 8 4])
error_shade(y_awake(:,wid),  err_awake(:,wid),  x, 'Color','g');
error_shade(y_anesth(:,wid), err_anesth(:,wid), x, 'Color','b');
h(1) = plot(x, y_awake(:,wid),  'g');
h(2) = plot(x, y_anesth(:,wid), 'b');
legend(h, {'Awake','Anesthetized'}, 'Location','northeastoutside', 'Box','off')

axis square
xlim([settings.dist_centers(1) settings.dist_centers(end)])
xlabel('Wall Distance (mm)')
ylabel('\Delta\kappa (a.u.)')
title([subjects{1} ' ' whisker_ids{wid} ' mean \Delta\kappa (wall window)'], 'FontWeight','normal', 'FontSize',6)
set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)

filename = ['FigS1_PanelH_mean_deltaKappa_by_distance_' whisker_ids{wid} '_' subjects{1}];
saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf');
saveas(gcf, fullfile(settings.fig_dir, filename), 'png');

%% Panel Fk: Across mice mean delta-kappa by distance (like Panel F midpoint)

clearvars -except settings
warning('off')

settings.varname = 'kappa';

base_dir = 'E:\PrV_Wall_Recordings';
example_dirs = {
    '20231204_sp\session_1';
    '20231212_sp\session_1';
    '20231212_sp\session_3';
    '20231213_sp\session_1';
    };
subjects = {'beh104', 'beh103', 'beh101', 'beh102'};
input_file = 'whiskerTracking_session.mat';

% Collect per-subject data
[mean_var_awake, mean_var_anesth] = deal(cell(numel(subjects),1));
for i = 1:numel(subjects)
    example_dir = fullfile(base_dir, example_dirs{i});
    disp(['(' num2str(i) '/' num2str(numel(subjects)) '): Loading whisker data from ' subjects{i}])
    load(fullfile(example_dir, input_file)); % loads 'session'
    [mean_var_awake{i}, mean_var_anesth{i}] = get_subject_data(session, settings);
end

% Distance labels
distance_names = arrayfun(@(x) [num2str(x) ' mm'], settings.dist_centers, 'Uni',0);

% Build tidy table + save .mat in settings.data_dir
T_kappa = build_whisker_source_table(subjects, mean_var_awake, mean_var_anesth, ...
    distance_names, settings.data_dir, settings.varname);

% Save per-subject CSVs (same pattern as midpoint)
unique_subjects = categories(T_kappa.subject);
for i = 1:numel(unique_subjects)
    subj = char(unique_subjects(i));
    Tsubj = T_kappa(T_kappa.subject == subj, :);
    out_csv = fullfile(settings.data_dir, ...
        sprintf('whisker_retraction_source_%s_%s.csv', settings.varname, subj));
    writetable(Tsubj, out_csv);
end

% Save combined CSV too (optional but convenient)
out_csv_all = fullfile(settings.data_dir, sprintf('whisker_retraction_source_%s_ALL.csv', settings.varname));
writetable(T_kappa, out_csv_all);

disp('Saved delta-kappa source tables:')
disp(out_csv_all)

% Concatenate across mice
num_whiskers = 3;
[y_awake_all, y_anesth_all] = deal(zeros(numel(subjects), numel(settings.dist_centers), num_whiskers));
for i = 1:numel(subjects)
    y_awake_all(i,:,:)  = squeeze(mean(mean_var_awake{i},1));
    y_anesth_all(i,:,:) = squeeze(mean(mean_var_anesth{i},1));
end

% Mean +/- SEM across subjects
y_awake  = squeeze(mean(y_awake_all,1));
err_awake = squeeze(std(y_awake_all,[],1) / sqrt(size(y_awake_all,1)));
y_anesth  = squeeze(mean(y_anesth_all,1));
err_anesth = squeeze(std(y_anesth_all,[],1) / sqrt(size(y_anesth_all,1)));
x = shiftdim(settings.dist_centers);


%%
% Plot
yl = [-1.5e-3 2e-4];

whisker_names = {'Greek Arc','1st Arc','2nd Arc'};
for wid = 1:num_whiskers
    figure(700+wid); clf; hold on
    set(gcf, 'Color','w', 'Units','centimeters', 'Position',[2 2 8 4])
    error_shade(y_awake(:,wid),  err_awake(:,wid),  x, 'Color','g');
    error_shade(y_anesth(:,wid), err_anesth(:,wid), x, 'Color','b');
    h(1) = plot(x, y_awake(:,wid),  'g');
    h(2) = plot(x, y_anesth(:,wid), 'b');
    legend(h, {'Awake','Anesthetized'}, 'Location','northeastoutside', 'Box','off')

    axis square
    xlim([settings.dist_centers(1) settings.dist_centers(end)])
    ylim(yl)
    xlabel('Wall Distance (mm)')
    ylabel('\Delta\kappa (mm^{-1})')
    title([whisker_names{wid} ' mean \Delta\kappa (across mice)'], 'FontWeight','normal', 'FontSize',6)
    set(gca, 'TickDir','out', 'Box','off', 'FontSize',6)
    set(gca, 'XTick',settings.dist_centers, 'XTickLabelRotation',0)

    filename = ['FigS1_PanelI_grand_mean_deltaKappa_by_distance-' whisker_names{wid}];
    saveas(gcf, fullfile(settings.fig_dir, filename), 'pdf');
    saveas(gcf, fullfile(settings.fig_dir, filename), 'png');
end


%% Functions

function error_shade(y, err, varargin)
    % Handle user inputs
    ip = inputParser();
    ip.addOptional('x', []);
    ip.addParameter('color', 'k');
    ip.addParameter('alpha', 0.3, @isnumeric);
    ip.parse(varargin{:});
    x = ip.Results.x;
    
    if isempty(x)
        x = repmat((1 : size(y,1))', 1, size(y,2));
    end
    
    for ii = 1 : size(x,2)
        patch([ x(:,ii); flip(x(:,ii)) ], [ y(:,ii)+err(:,ii); flip(y(:,ii)-err(:,ii)) ], ip.Results.color, ...
            'FaceAlpha',ip.Results.alpha, ...
            'LineStyle','none');
    end
end

function [var, vars] = get_whisker_data(session, varname)
    vars = struct;

    % Get traces from whisker data for these frames
    vars.theta_raw = session.whiskerData.theta(:,1:3);
    for w = size(vars.theta_raw, 2):-1:1
        % Interpolate NaN in theta
        nans = find(isnan(vars.theta_raw(:,w)));
        vars.theta_raw(nans,w) = interp1(1:numel(vars.theta_raw(:,w)), vars.theta_raw(:,w), nans, 'linear', NaN);
    end

    % Smooth theta a couple times to smooth over tracking errors, noise
    vars.theta_med = smoothdata(vars.theta_raw, 1, 'movmedian',5);
    vars.theta_smooth = smoothdata(vars.theta_med, 1, 'sgolay',21, 'degree',3);
    vars.theta_smooth = smoothdata(vars.theta_smooth, 1, 'sgolay',21, 'degree',3);

    vars.midpoint = smoothdata(vars.theta_smooth, 1, 'sgolay',101, 'degree',3);
    
%     % Re-calculate midpoint, amplitude, setpoint
%     for w = size(vars.theta_raw, 2):-1:1
%         % Get lower envelope
%         [pks, locs] = findpeaks(vars.theta_smooth(:,w));
%         vars.upper_envelope(:,w) = interp1(locs, pks, 1:numel(vars.theta_smooth(:,w)), 'linear', NaN);
%         [vals, locs] = findpeaks(-vars.theta_smooth(:,w));
%         vars.setpoint(:,w) = interp1(locs, -vals, 1:numel(vars.theta_smooth(:,w)), 'linear', NaN);
%     
%         % Get Amplitude from hilbert transform
%         vars.amplitude(:,w) = (vars.upper_envelope(:,w) - vars.setpoint(:,w)) / 2;
%     end
%     vars.midpoint = mean(cat(3, vars.setpoint, vars.upper_envelope), 3, 'omitnan');

%     fs = 500; % fps
%     [b,a] = butter(3, [2.5 25]/(fs/2), 'stop'); % bandstop filter
%     vars.midpoint = filtfilt(b, a, vars.midpoint); 

    % Get delta-kappa (curvature)
    vars.kappa = session.whiskerData.kappa;

    % Interpolate NaNs in kappa (tracking dropouts)
    for w = size(vars.kappa,2):-1:1
        nans = isnan(vars.kappa(:,w));
        if any(nans)
            vars.kappa(nans,w) = interp1(find(~nans), vars.kappa(~nans,w), find(nans), 'pchip', 'extrap');
        end
    end
    vars.kappa = smoothdata(vars.kappa, 1, 'movmedian',5);
    vars.kappa = smoothdata(vars.kappa, 1, 'sgolay',51, 'degree',3);

    % Flip sign (negative = bending toward posterior)
    vars.kappa = -vars.kappa;

    % Output variable by name requested
    var = vars.(varname);
end

function [mean_var_awake, mean_var_anesth, mean_baseline_awake, mean_baseline_anesth] = get_subject_data(session, settings)
    % Hard-coded settings
    num_passes = 20;
    num_whiskers = 3;
    num_conditions = 2;

    % Fix wrong condition label in session data
    session.dataTable.trialConditions{201} = 'anesth';
    
    % Process whisker data
    var = get_whisker_data(session, settings.varname);

    % Synchronize ephys with video frame times
    frame_times = session.frame_times;
    median_frame_duration = mean(diff(frame_times));
    frame_times = frame_times(1) : median_frame_duration : (frame_times(1) + (size(var,1)-1)*median_frame_duration);
    
    % Get pass indices for each distance, awake (1) vs. anesth (2), Cut to even size
    baseline_length = diff(settings.baseline_period) * settings.frame_rate;
    wall_length = diff(settings.wall_period) * settings.frame_rate;
    psth_length = diff(settings.psth_period) * settings.frame_rate;
    
    % Save psths as matrix (time, distance, condition)
    [wall_mean, baseline_mean] = deal(NaN(num_passes, numel(settings.dist_centers), num_whiskers, num_conditions));
    [psth_mean] = deal(NaN(psth_length, numel(settings.dist_centers), num_whiskers, num_conditions));
    for w = 1:num_whiskers
        whisker_id = w;
        for c = 1:num_conditions
            for d = 1:numel(settings.dist_centers)
                wd = unique(session.dataTable.wallDistanceMM);
                [~,wdex] = min(abs(wd-settings.dist_centers(d)));
                wdex = wd(wdex);
        
                % Get passes from awake
                if c==1 % awake
                    passIdx = session.dataTable.wallDistanceMM==wdex & cellfun(@(x) isempty(x), session.dataTable.trialConditions, 'Uni',1);
                else % anesthetized
                    passIdx = session.dataTable.wallDistanceMM==wdex & cellfun(@(x) strcmp(x,'anesth'), session.dataTable.trialConditions, 'Uni',1);
                end
                pass_start_times = session.pass_start_times(passIdx);
        
                [baseline, wall, psth] = deal(cell(numel(pass_start_times),1));
                for i = 1:numel(pass_start_times)
                    % Get baseline period
                    idx = frame_times >= pass_start_times(i)+settings.baseline_period(1) & ...
                        frame_times < pass_start_times(i)+settings.baseline_period(2);
                    baseline{i} = var(idx, whisker_id);
                
                    % Get wall period data
                    idx = frame_times >= pass_start_times(i)+settings.wall_period(1) & ...
                        frame_times < pass_start_times(i)+settings.wall_period(2);
                    wall{i} = var(idx, whisker_id);
                
                    % Get full psth period
                    idx = frame_times >= (pass_start_times(i)+settings.psth_period(1)) & ...
                        frame_times < (pass_start_times(i)+settings.psth_period(2));
                    psth{i} = var(idx, whisker_id);
                end
        
                % Interpolate if necessary
                for i = 1:numel(pass_start_times)
                    baseline{i} = interp1(1:numel(baseline{i}), baseline{i}, 1:baseline_length, 'linear', NaN);
                    wall{i} = interp1(1:numel(wall{i}), wall{i}, 1:wall_length, 'linear', NaN);
                    psth{i} = interp1(1:numel(psth{i}), psth{i}, 1:psth_length, 'linear', NaN);
                end
    
                % Get raw psth for each whisker
                tmp_psth  = cell2mat(psth)';
                psth_mean(:,d,w,c) = mean(tmp_psth,2,'omitnan');
        
                % Insert mean wall fr (across time) into matrix
                tmp_wall = cell2mat(wall)';
                wall_mean(:,d,w,c) = mean(tmp_wall,1,'omitnan');
        
                % Insert mean baseline fr into matrix
                tmp_baseline = cell2mat(baseline)';
                baseline_mean(:,d,w,c) = mean(tmp_baseline,1,'omitnan');
            end
        end
    end
    
    % Get mean baseline period (average across time, repetitions, distances)
    mean_baseline = squeeze(mean(mean(mean(baseline_mean,1,'omitnan'),1,'omitnan')));
    
    % Get mean (baseline-subtracted) delta-theta for each distance and each whisker
    clear mean_fr_baseline_awake mean_fr_baseline_anesth mean_var_awake mean_var_anesth
    clear mean_baseline_awake mean_baseline_anesth
    for w = num_whiskers:-1:1
        mean_var_awake(:,:,w) = squeeze(wall_mean(:,:,w,1)) - mean_baseline(w,1);
        mean_var_anesth(:,:,w) = squeeze(wall_mean(:,:,w,2)) - mean_baseline(w,2);
        mean_baseline_awake(w) = mean_baseline(w,1);
        mean_baseline_anesth(w) = mean_baseline(w,2);
    end
end

function T = build_whisker_source_table(subjects, mean_var_awake, mean_var_anesth, distance_names, save_dir, varname)
    % BUILD_WHISKER_SOURCE_TABLE
    % Convert whisker retraction data into a tall / tidy table for stats.
    %
    % Inputs:
    %   subjects          : cellstr, {nSubjects x 1}
    %   mean_var_awake    : cell {nSubjects x 1}, each [reps x dist x whisk]
    %   mean_var_anesth   : cell {nSubjects x 1}, same size/shape as awake
    %   distance_names    : cellstr {1 x nDist} or {nDist x 1}, e.g. {'5 mm','7 mm',...}
    %   save_dir          : char, directory where to save .mat
    %   varname           : char, e.g. 'midpoint' (used in filename)
    %
    % Output:
    %   T : table with columns:
    %       subject (categorical)
    %       rep (double)
    %       distance (categorical, natural sorted)
    %       whiskerID (categorical)
    %       condition (categorical)
    %       [value] (double)

    if nargin < 6
        error('Need subjects, mean_var_awake, mean_var_anesth, distance_names, save_dir, varname');
    end

    % Collectors
    all_subject   = {};
    all_rep       = [];
    all_distance  = {};
    all_whisker   = [];
    all_condition = {};
    all_value     = [];

    % Ensure row vector distance labels
    dist_labels = distance_names(:)';

    for s = 1:numel(subjects)
        awake_dat  = mean_var_awake{s};
        anesth_dat = mean_var_anesth{s};

        if isempty(awake_dat) || isempty(anesth_dat)
            warning('Subject %s has empty data. Skipping.', subjects{s});
            continue
        end

        [nRep, nDist, nWhisk] = size(awake_dat);

        if numel(dist_labels) ~= nDist
            error('distance_names length (%d) does not match nDist (%d)', numel(dist_labels), nDist);
        end

        cond_names = {'awake','anesthetized'};
        cond_data  = {awake_dat, anesth_dat};

        for c = 1:numel(cond_names)
            this_dat = cond_data{c}; % [rep x dist x whisk]
            [rep_idx, dist_idx, whisk_idx] = ndgrid(1:nRep, 1:nDist, 1:nWhisk);

            all_subject   = [all_subject;   repmat(subjects(s),    numel(this_dat), 1)];
            all_rep       = [all_rep;       rep_idx(:)];
            all_distance  = [all_distance;  dist_labels(dist_idx(:))'];
            all_whisker   = [all_whisker;   whisk_idx(:)];
            all_condition = [all_condition; repmat(cond_names(c), numel(this_dat), 1)];
            all_value     = [all_value;     this_dat(:)];
        end
    end

    % Natural sort of distance labels
    try
        % Use natsort from MATLAB File Exchange if available
        [~, sort_idx] = natsort(distance_names);
        sorted_distances = distance_names(sort_idx);
    catch
        % Fallback: extract numeric portion and sort numerically
        nums = regexp(distance_names, '([\d\.]+)', 'tokens', 'once');
        nums = cellfun(@(x) str2double(x{1}), nums);
        [~, sort_idx] = sort(nums);
        sorted_distances = distance_names(sort_idx);
    end

    % Build table
    T = table;
    T.subject    = categorical(string(all_subject));
    T.rep        = all_rep;
    T.distance   = categorical(string(all_distance), sorted_distances);
    T.whiskerID  = categorical(all_whisker);
    T.condition  = categorical(string(all_condition));

    % Name the value column based on varname
    switch lower(varname)
        case {'midpoint','theta','theta_mid','theta_midpoint'}
            value_col = 'deltaMidpoint';
        case {'kappa','curvature','delta_kappa','deltakappa'}
            value_col = 'deltaKappa';
        otherwise
            value_col = ['delta_' matlab.lang.makeValidName(varname)];
    end
    
    T.(value_col) = all_value;


    % Save
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    out_file = fullfile(save_dir, ['whisker_retraction_source_' varname '.mat']);
    save(out_file, 'T');

    % Print brief summary
    disp('Summary (mean deltaAngle by subject/condition/distance):');
    summary_tbl = groupsummary(T, {'subject','condition','distance'}, "mean", value_col);

    disp(summary_tbl);
end







