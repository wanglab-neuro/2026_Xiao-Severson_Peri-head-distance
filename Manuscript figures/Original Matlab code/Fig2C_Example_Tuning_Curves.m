%% NAIVE PrV silicon probe units
clear; clc
settings.region = 'PrV';
settings.dataset = 'naive';
settings.output_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData';
settings.fig_dir = 'E:\PrV_Wall_Recordings\Figures';
settings.data_dir = 'F:\PrV_Wall_Recordings';
settings.rec_dirs = {
    '20221212_sp\session_2'
    '20221213_sp\session_2'
    '20230131_sp\session_2'
    '20230202_sp\session_1'
    '20230202_sp\session_3'
    '20230203_sp\session_2'
    '20230208_sp\session_1'
    '20230209_sp\session_2'
    '20230210_sp\session_2'
    '20230215_sp\session_1'
    '20230215_sp\session_2'
    '20230216_sp\session_2'
    '20230217_sp\session_1'
    '20230217_sp\session_2'
    '20230222_sp\session_2'
    '20230223_sp\session_2'
    '20230225_sp\session_2'
    '20230226_sp\session_2'
    '20230226_sp\session_4'
    '20230227_sp\session_1'
    '20230227_sp\session_2'
    '20230228_sp\session_1'
    '20230228_sp\session_2'
    '20230321_sp\session_1'
    '20230321_sp\session_2'
    '20230322_sp\session_1'
    '20230322_sp\session_2'
    '20230323_sp\session_1'
    '20230323_sp\session_2'
    '20230326_sp\session_2'
    '20230327_sp\session_2'
    '20230328_sp\session_1'
    '20230328_sp\session_3'
    '20230329_sp\session_1'
    '20230329_sp\session_2'
    '20230329_sp\session_4'
    '20230330_sp\session_1'
    '20250320_sp\session_1'
    '20250320_sp\session_2'
    '20250324_sp\session_1'
    '20250326_sp\session_1'
    '20250327_sp\session_1'
    };

%% NAIVE PrV Whisker Cutting silicon probe Units:
 clear; clc
 settings.region = 'PrV';
 settings.dataset = 'naive';
 settings.output_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\Whisker_Cutting_naive';
 settings.fig_dir = 'E:\PrV_Wall_Recordings\Figures\Whisker_Cutting_naive';
 settings.data_dir = 'F:\PrV_Wall_Recordings';
 settings.rec_dirs = {
    '20230525_sp\session_1'
    '20230602_sp\session_1'
    '20230608_sp\session_1'
    '20230615_sp\session_1'
    '20230626_sp\session_1'
    '20230627_sp\session_1'
    '20230711_sp\session_1'
    '20230725_sp\session_1'
 };

%% NAIVE PrV One-Whisker Cutting silicon probe Units:
clear; clc
settings.region = 'PrV';
settings.dataset = 'onewhisker_naive';
settings.output_path = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\One_Whisker_naive';
settings.fig_dir = 'E:\PrV_Wall_Recordings\Figures\One_Whisker_naive';
settings.data_dir = 'F:\PrV_Wall_Recordings';
settings.rec_dirs = {
     '20240522_sp\session_1'
     '20240624_sp\session_1'
     '20240626_sp\session_1'
     '20240702_sp\session_1'
     '20240703_sp\session_1'
     '20240724_sp\session_1'
     '20240726_sp\session_1'
     '20240814_sp\session_1'
     '20240815_sp\session_1'
     '20240905_sp\session_1'
     '20240909_sp\session_1'
      '20240923_sp\session_1' 
     '20240924_sp\session_1' 
     '20240925_sp\session_1' 
     '20240926_sp\session_1'  
     '20250325_sp\session_2'
     '20250327_sp\session_2'
};

%% Get Baseline & Wall FR, Plot PSTH (prerequisite for tuning)

% Load uniform distance set (7:2:23 mm)
dist_centers = matfile('E:\Code\Analysis\NeuralDecoder\dist_centers.mat').dist_centers; pause(0.001);
dist_edges = (dist_centers(1) - 1) : 2 : (dist_centers(end) + 1);

if strcmpi(settings.dataset, 'tungsten')
    % Bin distances empirically into 1-mm wide bins
    use_distance_centers = false; % 
    dist_centers = 1:23;
    dist_edges = 0.5:23.5;
else
    % Bin distances into 2-mm wide bins
    use_distance_centers = true; % false; % 
end

% Pre-allocate output structure
psth = cell(numel(settings.rec_dirs),1);
for f = 1:numel(settings.rec_dirs)
    close all
    clearvars -except settings psth f use_distance_centers dist_centers dist_edges

    disp(['Rec: ' num2str(f)])
    session = matfile(fullfile(settings.data_dir,settings.rec_dirs{f},'session.mat')).session;
    recID = strrep(session.recID,'_',' ');
    mouseID = strrep(session.mouseID,'_',' ');
    disp(session.recID)

    % Fix random seed
    rng('Default')
    rng(42)

    % Settings
    show_plots = 1; % faster if you don't plot
    save_plots = 1;
    save_table = 1;
    bin_size = 0.2; % seconds 
    baseline_period = [-4.5 0]; % changed from [-4.5 1.5] 20250728 to exclude Zaber artifacts
    wall_period = [3 14];
    boot_num = 1000;
    ci = [2.5 97.5]; % [1 99];

    data = struct;
    if contains(settings.rec_dirs{f},'_t')
        for p = 1:size(session.dataTable.spikeTimes,1)
            session.dataTable.spikeTimes{p} = {session.dataTable.spikeTimes{p}};
        end
    end
    if ~iscell(session.sortedSpikeTimes)
        session.sortedSpikeTimes = {session.sortedSpikeTimes};
    end
    for u = 1:numel(session.dataTable.spikeTimes{1})
        pause(0.5)
        tic
        disp(u)
        st = session.sortedSpikeTimes;
        time_edges = arrayfun(@(x) (x + baseline_period(1)) : bin_size : (x + wall_period(2)), session.pass_start_times, 'Uni',0);
        fr = cellfun(@(e) histcounts(st{u},e) / bin_size, time_edges, 'Uni',0);
        [min_edges,min_ind] = min(cell2mat(cellfun(@(x) numel(x)-1, time_edges, 'Uni',0)),[],'omitnan');

        % time = time_edges{min_ind}(1:min_edges+1);
        time = baseline_period(1) : bin_size : wall_period(2);
        base_time = time >= baseline_period(1) & time <= baseline_period(2);
        wall_time = time >= wall_period(1) & time <= wall_period(2);

        time = mean([time(1:end-1); time(2:end)]);
        base_time = max([base_time(1:end-1); base_time(2:end)]);
        wall_time = max([wall_time(1:end-1); wall_time(2:end)]);

        base_time = base_time(1:min_edges);
        wall_time = wall_time(1:min_edges);

        % If a perturbation experiment was performed (e.g., whisker cutting 
        % or laser stimulation), parse trial conditions
        if isfield(session,'trial_conditions')
            null = cell2mat(cellfun(@(x) isempty(x), session.dataTable.trialConditions, 'Uni',0));
            fr_condition = fr(~null);
            fr_shuf_condition = datasample(fr_condition, numel(fr_condition), 'Replace',false);
            for i = numel(fr_condition):-1:1
                fr_cut_condition(i,:) = fr_condition{i}(1:min_edges);
                fr_cut_shuf_condition(i,:) = fr_shuf_condition{i}(1:min_edges);
            end

            % Compute z-scored firing rate
            fr_base = fr_cut_condition(:,base_time);
            zfr_cut_condition = (fr_cut_condition - mean(fr_base(:),'omitnan')) / std(fr_base(:),'omitnan');
            
            % Mask trials for null period
            fr = fr(null);
        else
            null = true(size(fr));
        end

        fr_shuf = datasample(fr, numel(fr), 'Replace',false);
        for i = numel(fr):-1:1
            fr_cut(i,:) = fr{i}(1:min_edges);
            fr_cut_shuf(i,:) = fr_shuf{i}(1:min_edges);
        end

        % Compute z-scored firing rate
        fr_base = fr_cut(:,base_time);
        fr_all = fr_cut(:,base_time | wall_time);
        zfr_cut = (fr_cut - mean(fr_base(:),'omitnan')) / std(fr_base(:),'omitnan');

        % For each distance, compute mean firing rate (tuning) over each period
        wallDist = session.dataTable.wallDistanceMM;
        % wallDist = session.dataTable.wallVelocityMMperSec; 

        if strcmpi(settings.dataset, 'laser')
            % Use empirical presented distances directly
            wd = sort(unique(wallDist(~isnan(wallDist))), 'ascend');
            dist_centers_used = wd;
        else
            wd = NaN(1,numel(dist_centers));
            distances = unique(wallDist);
            for d = numel(dist_centers):-1:1
                if any(distances > dist_edges(d)) && any(distances < dist_edges(d+1))
                    [~, nearest_neighbor] = min(abs(distances - dist_centers(d)));
                    wd(d) = distances(nearest_neighbor);
                else
                    wd(d) = NaN;
                end
            end
            nan_dist = isnan(wd);
            wd = wd(~nan_dist);
            dist_centers_used = dist_centers(~nan_dist);
        end

        clear mean_fr n_fr std_fr sem_fr
        mean_fr_boot = NaN(numel(wd),min_edges,boot_num);
        for w = numel(wd):-1:1
            passInd = wallDist(null)==wd(w);
            % Pool PSTHs into distance bins
            % passInd = (wallDist(null) > (dist_centers_used(w) - 1)) & ...
            %     (wallDist(null) < (dist_centers_used(w) + 1));
            passIdx = find(passInd);
            n_fr(w,1) = sum(passInd);
            nx = sqrt(n_fr(w));

            % Compute PSTH from true responses
            mean_fr(w,:) = mean(fr_cut(passInd,:),'omitnan');
            sem_fr(w,:) = std(fr_cut(passInd,:),[],'omitnan') / nx;

            % Compute PSTH from shuffled responses
            mean_fr_shuf(w,:) = mean(fr_cut_shuf(passInd,:),'omitnan');
            sem_fr_shuf(w,:) = std(fr_cut_shuf(passInd,:),[],'omitnan') / nx;

            % Compute bootstrapped mean to compute confidence intervals
            resamp = randi(numel(passIdx), [numel(passIdx) boot_num]);
            for b = 1:boot_num
                mean_fr_boot(w,:,b) = mean(fr_cut(passIdx(resamp(:,b)),:),'omitnan');
                % mean_fr_shuf_boot(w,:,b) = mean(fr_cut_shuf(pIdx(resamp(:,b)),:),'omitnan');
            end
            ci_fr(w,:,:) = prctile(mean_fr_boot(w,:,:),ci,3);

            if isfield(session,'trial_conditions')
                passInd = wallDist(~null)==wd(w);
                passIdx = find(passInd);
                n_fr(w,1) = sum(passInd);
                nx = sqrt(n_fr(w));

                % Compute PSTH from true responses
                mean_fr_condition(w,:) = mean(fr_cut_condition(passInd,:),'omitnan');
                sem_fr_condition(w,:) = std(fr_cut_condition(passInd,:),[],'omitnan') / nx;
    
                % Compute PSTH from shuffled responses
                mean_fr_shuf_condition(w,:) = mean(fr_cut_shuf_condition(passInd,:),'omitnan');
                sem_fr_shuf_condition(w,:) = std(fr_cut_shuf_condition(passInd,:),[],'omitnan') / nx;
            end
        end

        % Compute overall trial average and bootstrapped samples
        mean_fr_all(:,1) = mean(fr_cut,'omitnan');
        mean_fr_boot_all = squeeze(mean(mean_fr_boot,'omitnan'));
        % mean_fr_shuf_boot_all = squeeze(mean(mean_fr_shuf_boot,'omitnan'));

        % Gather data for baseline and wall periods
        mean_fr_baseline = mean_fr(:, base_time);
        sem_fr_baseline = sem_fr(:, base_time);
        mean_fr_wall = mean_fr(:, wall_time);
        sem_fr_wall = sem_fr(:, wall_time);

        if isfield(session, 'trial_conditions')
            % Gather data for baseline and wall periods
            mean_fr_baseline_condition = mean_fr_condition(:, base_time);
            sem_fr_baseline_condition = sem_fr_condition(:, base_time);
            mean_fr_wall_condition = mean_fr_condition(:, wall_time);
            sem_fr_wall_condition = sem_fr_condition(:, wall_time);
        end

        % Compute confidence intervals
        mean_fr_all_baseline = mean(mean_fr_all(base_time),1,'omitnan');
        mean_fr_all_wall = mean(mean_fr_all(wall_time),1,'omitnan');
        ci_fr_all_baseline = squeeze(prctile(mean(mean_fr_boot_all(base_time,:),1,'omitnan'),ci,2));
        ci_fr_all_wall = squeeze(prctile(mean(mean_fr_boot_all(wall_time,:),1,'omitnan'),ci,2));
        ci_fr_baseline = squeeze(prctile(mean(mean_fr_boot(:,base_time,:),2,'omitnan'),ci,3));
        ci_fr_wall = squeeze(prctile(mean(mean_fr_boot(:,wall_time,:),2,'omitnan'),ci,3));

        if show_plots || save_plots
            colors = flipud(turbo(size(mean_fr,1)));
            figure(u); clf;
            set(gcf,'Color','w', 'Renderer','painters', 'Position',[200 200 800 400])
            condition_names = {'none'};
            pause(0.001)
            if isfield(session, 'trial_conditions')
                condition_names = unique(session.trial_conditions(~null));
                tl = tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact');
                if strcmp(condition_names{1},'cut')
                    layout_title = 'intact';
                else 
                    layout_title = ['no ' condition_names{1}];
                end
                pause(0.1)
            else
                layout_title = '';
                tl = tiledlayout(1,1, 'TileSpacing','compact', 'Padding','compact');
            end
            pause(0.001)
            b = nexttile;
            title(b, layout_title);
            hold on
            for w = size(mean_fr,1):-1:1
                MPlot.ErrorShade(mean_fr(w,:), sem_fr(w,:), time, 'Color',colors(w,:)); % shade SEM
%                 MPlot.ErrorShade(mean_fr(w,:), mean_fr(w,:) - squeeze(ci_fr(w,:,1)), time, 'Color',colors(w,:)); % shade CI
                h2(w) = plot(time, mean_fr(w,:), 'Color',colors(w,:), 'LineWidth',1);
            end
            set(gca, 'TickDir','out', 'box','off', 'Color','none')
            if ~isfield(session,'channelID')
                session.channelID = 1;
            end
            if u <= numel(session.channelID)
                main_title = strjoin({mouseID, recID, 'Unit', num2str(u), 'Ch', num2str(session.channelID(u))}, ' ');
            else
                disp('Invalid index for session.channelID');
                main_title = strjoin({mouseID, recID, 'Unit', num2str(u)}, ' ');
            end
            ylabel('Firing rate (sp/s)')
            xlabel('Time (s)')
            xlim([0, max(time)])
            max_sem = max(max(sem_fr,[],'omitnan'));

            ymin1 = min(min(mean_fr,[],'omitnan')) - max_sem;
            ymax1 = max(max(mean_fr,[],'omitnan')) + max_sem;
            ymin1 = max([0 ymin1]);
            ylim([ymin1 ymax1])

            % Plot PSTHs from perturbation trials in side-by-side subplot
            if isfield(session, 'trial_conditions')
                % Set same y-limits for each subplot
                ymin1 = max([0 ymin1]);
                
                ymin2 = max([0 min(min(mean_fr_condition,[],2,'omitnan')) - max(max(sem_fr_condition,[],'omitnan'),[],'omitnan')]);
                ymax2 = max(max(mean_fr_condition,[],'omitnan'),[],'omitnan') + max(max(sem_fr_condition,[],'omitnan'),[],'omitnan');
                
                ymin = min([ymin1 ymin2],[],'omitnan');
                ymax = max([ymax1 ymax2],[],'omitnan');

                ylim([ymin max([ymin+0.001,ymax],[],'omitnan')])

                b = nexttile;
                hold on
                for w = size(mean_fr_condition,1):-1:1
                    MPlot.ErrorShade(mean_fr_condition(w,:), sem_fr_condition(w,:), time, 'Color',colors(w,:)); % shade SEM
                    h2(w) = plot(time, mean_fr_condition(w,:), 'Color',colors(w,:), 'LineWidth',1);
                end
                set(gca, 'TickDir','out', 'box','off', 'YTickLabels',[], 'Color','none')
                title(b, condition_names{1})
                try
                    if u <= numel(session.channelID)
                        title(b, condition_names{1})
                    else
                        disp('Invalid index for session.channelID');
                    end
                end
                % title(strjoin({session.mouseID,recID,'Unit',num2str(u),'Ch',num2str(session.channelID(u))},' '))
                xlabel('Time (s)')
                xlim([0, max(time)])
                ylim([ymin max([ymin+0.001,ymax],[],'omitnan')])
                
                pause(0.1)
            end

            title(tl, main_title)

            h = copyobj(h2, gca);
            set(h, 'XData',NaN, 'YData',NaN, 'LineWidth',2)
            if strcmpi(settings.dataset, 'laser')
                distance_labels = wd;
            elseif use_distance_centers
                distances_used = find(histcounts(wd, dist_edges)>0);
                distance_labels = dist_centers(distances_used);
            else
                distance_labels = floor(wd);
            end
            legend(h, num2str(shiftdim(distance_labels)), ...
                'Location','northeastoutside', ...
                'Position',[0.9 0.3 0.1 0.7], ...
                'Color','None', ...
                'EdgeColor','None') % 
            pause(0.1)

            % Save plots
            if save_plots
                filename = [session.recID '_unit' num2str(u) '_PSTH'];
                if strcmp(condition_names{1}, 'laser')
                    filefolder = fullfile(settings.fig_dir,'Laser','PSTH');
                elseif strcmp(condition_names{1}, 'cut')
                    filefolder = fullfile(settings.fig_dir,'Whisker_Cutting','PSTH');
                else
                    filefolder = fullfile(settings.fig_dir,'PSTH');
                end
                filepath = fullfile(filefolder, filename);
                if ~isfolder(filefolder); mkdir(filefolder); end
                savefig(gcf, filepath)
                saveas(gcf, filepath, 'png')
                saveas(gcf, filepath, 'pdf')
                pause(0.1)
            end
        end

        % Compile unit psth data
        data.time(u,:) = time;
        data.baseline_period(u,:) = baseline_period;
        data.wall_period(u,:) = wall_period;
        data.wall_dist(u,:) = wd;
        data.ci(u,:) = ci;
        data.boot_num(u,:) = boot_num;

        data.fr_cut(u,:,:) = fr_cut;
        data.zfr_cut(u,:,:) = zfr_cut;

        data.mean_fr_all_baseline(u,:) = mean_fr_all_baseline;
        data.ci_fr_all_baseline(u,:) = ci_fr_all_baseline;
        data.mean_fr_all_wall(u,:) = mean_fr_all_wall;
        data.ci_fr_all_wall(u,:) = ci_fr_all_wall;

        data.mean_fr(u,:,:) = mean_fr;
        data.sem_fr(u,:,:) = sem_fr;
        data.mean_fr_shuf(u,:,:) = mean_fr_shuf;
        data.sem_fr_shuf(u,:,:) = sem_fr_shuf;

        data.fr_baseline(u,:,:) = mean_fr_baseline;
        data.sem_fr_baseline(u,:,:) = sem_fr_baseline;

        data.mean_fr_baseline(u,:) = mean(mean_fr_baseline,2,'omitnan');
        data.ci_fr_baseline(u,:,:) = ci_fr_baseline;

        data.fr_wall(u,:,:) = mean_fr_wall;
        data.sem_fr_wall(u,:,:) = sem_fr_wall;

        data.mean_fr_wall(u,:) = mean(mean_fr_wall,2,'omitnan');
        data.ci_fr_wall(u,:,:) = ci_fr_wall;

        if isfield(session, 'trial_conditions')
            % Compile data if perturbation was performed
            data.fr_cut_condition(u,:,:) = fr_cut_condition;
            data.zfr_cut_condition(u,:,:) = zfr_cut_condition;
    
            data.mean_fr_condition(u,:,:) = mean_fr_condition;
            data.sem_fr_condition(u,:,:) = sem_fr_condition;
            data.mean_fr_shuf_condition(u,:,:) = mean_fr_shuf_condition;
            data.sem_fr_shuf_condition(u,:,:) = sem_fr_shuf_condition;
        end

        toc
    end

    % Compile recording psth data
    psth{f}.recID = session.recID;
    psth{f}.data = data;
    pause(0.001)
end

disp('Saving PSTH, Z-Score, and Bootstrap Data...')

% Concatenate PSTH to same dimensions as tuning
if exist('psth','var')==1
    psth2 = psth;
    psth_cell = psth2;
    psth = struct;
    recID = [];
    fields = fieldnames(psth_cell{1}.data);
    cell_list = {'time','wall_dist','fr_cut','zfr_cut', ...
        'mean_fr','sem_fr','mean_fr_shuf','sem_fr_shuf','fr_baseline', ...
        'sem_fr_baseline','mean_fr_baseline','ci_fr_baseline', ...
        'fr_wall','sem_fr_wall','mean_fr_wall','ci_fr_wall'};
    if isfield(session, 'trial_conditions')
        cell_list = cat(2, ...
            cell_list,...
            { ...
            'fr_cut_condition','zfr_cut_condition','mean_fr_condition', ...
            'sem_fr_condition','mean_fr_shuf_condition','sem_fr_shuf_condition', ...
            });
    end
    for f = 1:numel(fields)
        psth.(fields{f}) = [];
    end
    for r = 1:numel(psth_cell)
        num_units = size(psth_cell{r}.data.(fields{1}),1);
        recID = cat(1,recID,repmat({psth_cell{r}.recID},num_units,1));
        for f = 1:numel(fields)
            tmp = psth_cell{r}.data.(fields{f});
            if any(strcmp(cell_list,fields{f}))
                % psth.(fields{f}) = cat(1, psth.(fields{f}), num2cell(tmp,[2:ndims(tmp)]));
                clear tmp2
                for u = size(tmp,1):-1:1
                    tmp2{u,1} = squeeze(tmp(u,:,:));
                end
                psth.(fields{f}) = cat(1, psth.(fields{f}), tmp2);
            else
                psth.(fields{f}) = cat(1, psth.(fields{f}), tmp);
            end
        end
    end
    psth.recID = recID;
    if ~isfolder(settings.output_path)
        mkdir(settings.output_path)
    end

    if save_table
        save(fullfile(settings.output_path,'psth_cat_PrV'), 'psth', '-v7.3')
        save(fullfile(settings.output_path,'psth_cell_PrV'), 'psth_cell', '-v7.3')
    end
end

disp('Done!')

%% Generate Tuning Curve Data and Plots

clearvars -except settings

% Load uniform distance set (7:2:23 mm)
dist_centers = matfile('E:\Code\Analysis\NeuralDecoder\dist_centers.mat').dist_centers;
dist_edges = (dist_centers(1) - 1) : 2 : (dist_centers(end) + 1);

if strcmpi(settings.dataset, 'laser')
    % For laser experiments, use the empirical distances directly
    use_distance_centers = false;
elseif strcmpi(settings.dataset, 'tungsten')
    % Bin distances empirically into 1-mm wide bins
    use_distance_centers = false; % 
    dist_centers = 1:23;
    dist_edges = 0.5:23.5;
else
    % Bin distances into 2-mm wide bins
    use_distance_centers = true; % false; % 
end

% Fix random seed
rng('Default')
rng(42)

% Settings 
settings.base_time = [-4.5 0];
settings.time_win = [3 11];
settings.show_plots = 0;
settings.save_plots = 1;
settings.save_table = 1;
settings.compute_ci = 0; % 1; %
settings.plot_ci = 0; % 1; % else plot SEM as error bars
settings.nboot = 1e4;
settings.ci = [2.5 97.5]; % [1 99]; % 
cnt = 1;

psth = matfile(fullfile(settings.output_path,'psth_cat_PrV')).psth;
tuning = struct;

for f = 1:numel(settings.rec_dirs)
    tic
    close all
    clearvars -except settings psth f cnt tuning dist_centers dist_edges use_distance_centers

    disp(f)
    session = matfile(fullfile(settings.data_dir,settings.rec_dirs{f},'session.mat')).session;
    recID = strrep(session.recID,'_',' ');

    rec_name = strrep(strrep(settings.rec_dirs{f}, filesep, '_'), '/', '_');
    rec_rows = find(strcmp(psth.recID, rec_name));

    num_conditions = 1;
    if isfield(session, 'trial_conditions')
        num_conditions = 2;
        null = cell2mat(cellfun(@(x) isempty(x), session.dataTable.trialConditions, 'Uni',0));
        condition_name = unique(session.trial_conditions(~null));
        if strcmp(condition_name,'laser')
            condition_names = cat(1, {'no laser'}, condition_name);
        elseif strcmp(condition_name,'cut')
            condition_names = cat(1, {'intact'}, condition_name);
        elseif strcmp(condition_name,'anesth')
            condition_names = cat(1, {'awake'}, condition_name);
        end
    else
        null = true(size(session.dataTable.trialNum));
        condition_names = cell(1);
    end
    
    % Handle cell for single electrode (Tungsten) recordings
    if ~iscell(session.dataTable.spikeTimes{1})
        session.dataTable.spikeTimes = cellfun(@(x) {x}, session.dataTable.spikeTimes, 'Uni',0);
    end

    for u = 1:numel(session.dataTable.spikeTimes{1})
        for c = 1:num_conditions
            if c==1
                condition = null;
            elseif c==2 
                condition = ~null;
            end

            % Plot forward pass wall distance tuning (average spike rate during pass period)
            st = session.dataTable.spikeTimes;
            edges = settings.time_win;
            fr = cell2mat(cellfun(@(x) histcounts(x{u},edges) / diff(edges), st, 'Uni',0));
            fr_shuf = datasample(fr, numel(fr), 'Replace',false);
            fr = fr(condition);
            
            wallDist = session.dataTable.wallDistanceMM(condition);
            
            % Get time bins for baseline
            edges = settings.base_time;
            fr_baseline = cell2mat(cellfun(@(st) histcounts(st{u},edges) / diff(edges), st, 'Uni',0));
            time = psth.time{rec_rows(u)};
    
            % Get z-scored PSTH
            if c==1
                zfr = squeeze(psth.zfr_cut{rec_rows(u)});
            else
                zfr = squeeze(psth.zfr_cut_condition{rec_rows(u)});
            end
            zfr_shuf = datasample(zfr(:), numel(zfr), 'Replace',false);
            zfr_shuf = reshape(zfr_shuf, size(zfr,1), size(zfr,2));

            % Compute area-under-curve of z-scored PSTH
            % preserves time-dependent modulation of firing rates
            % auc = mean(abs(zfr(:,time >= settings.time_win(1) & time <= settings.time_win(2))),2,'omitnan');
            % auc_shuf = mean(abs(zfr_shuf(:,time >= settings.time_win(1) & time <= settings.time_win(2))),2,'omitnan');

            % Compute z-scored tuning curves from z-scored PSTH
            zfr = mean(zfr(:,time >= settings.time_win(1) & time <= settings.time_win(2)),2,'omitnan');
            zfr_shuf = mean(zfr_shuf(:,time >= settings.time_win(1) & time <= settings.time_win(2)),2,'omitnan');

            if strcmpi(settings.dataset, 'laser')
                % Use empirical presented distances directly
                wd = sort(unique(wallDist(~isnan(wallDist))), 'ascend');
                dist_centers_used = wd;
            else
                wd = NaN(1,numel(dist_centers));
                distances = unique(wallDist);
                for d = numel(dist_centers):-1:1
                    if any(distances > dist_edges(d)) && any(distances < dist_edges(d+1))
                        [~, nearest_neighbor] = min(abs(distances - dist_centers(d)));
                        wd(d) = distances(nearest_neighbor);
                    else
                        wd(d) = NaN;
                    end
                end
                nan_dist = isnan(wd);
                wd = wd(~nan_dist);
                dist_centers_used = dist_centers(~nan_dist);
            end

            % Fix random seed
            rng('Default')
            rng(42)
            
            if settings.compute_ci
                [fr_base_boot, fr_boot, zfr_boot] = deal(zeros(numel(wd),settings.nboot));
            end

            % max_pass = histcounts(wallDist);
            for w = numel(wd):-1:1
                passInd = wallDist==wd(w);
                % passInd = (wallDist > (dist_centers_used(w) - 1)) & ...
                %     (wallDist < (dist_centers_used(w) + 1));
                n_fr(w) = sum(passInd);
                nx = sqrt(n_fr(w));

                % Get mean firing rate for each distance during baseline
                fr_base{w} = fr_baseline(passInd);

                % Compute mean baseline firing
                mean_fr_base(w) = mean(fr_baseline(passInd),'omitnan');
                sem_fr_base(w) = std(fr_baseline(passInd),[],'omitnan') / nx;
    
                % Get mean firing rate for each distance during wall pass
                pass_ind{w} = find(passInd);
                fr_rep{w} = fr(passInd);
                
                % Compute tuning curve from true responses
                mean_fr(w) = mean(fr(passInd),'omitnan');
                sem_fr(w) = std(fr(passInd),[],'omitnan') / nx;
    
                % Compute tuning curve from shuffled responses
                mean_fr_shuf(w) = mean(fr_shuf(passInd),'omitnan');
                sem_fr_shuf(w) = std(fr_shuf(passInd),[],'omitnan') / nx;

                % Get z-scored sample of each distance
                zfr_rep{w} = zfr(passInd);
    
                % Compute z-scored tuning curve from true responses
                mean_zfr(w) = mean(zfr(passInd),'omitnan');
                sem_zfr(w) = std(zfr(passInd),[],'omitnan') / nx;
    
                % Compute z-scored tuning curve from shuffled responses
                mean_zfr_shuf(w) = mean(zfr_shuf(passInd),'omitnan');
                sem_zfr_shuf(w) = std(zfr_shuf(passInd),[],'omitnan') / nx;
                
                % Get auc sample of each distance
                % auc_rep{w} = auc(passInd);

                % Compute z-scored tuning curve from true responses
                % mean_auc(w) = mean(auc(passInd),'omitnan');
                % sem_auc(w) = std(auc(passInd),[],'omitnan') / nx;
    
                % Compute z-scored tuning curve from shuffled responses
                % mean_auc_shuf(w) = mean(auc_shuf(passInd),'omitnan');
                % sem_auc_shuf(w) = std(auc_shuf(passInd),[],'omitnan') / nx;

                % Fix random seed
                rng('Default')
                rng(42)

                % Compute confidence intervals for baseline and pass
                if settings.compute_ci
                    pind = find(passInd);
                    nfr = n_fr(w);
                    [fr_base_boots, fr_boots, zfr_boots, auc_boots] = ...
                        deal(zeros(nfr,settings.nboot));
                    for n = 1:settings.nboot
                        frn = nfr;
                        frb = fr_baseline;
                        frw = fr;
                        frz = zfr;
                        % fra = auc;
                        idx = datasample(pind, frn, 'Replace',true);
                        fr_boots(:,n) = frw(idx);
                        zfr_boots(:,n) = frz(idx);
                        % auc_boots(:,n) = fra(idx);
                        fr_base_boots(:,n) = frb(idx);
                        fr_base_boot(w,:) = mean(fr_base_boots,'omitnan');
                    end
                    
                    fr_boot(w,:) = mean(fr_boots,'omitnan');
                    zfr_boot(w,:) = mean(zfr_boots,'omitnan');
                    % auc_boot(w,:) = mean(auc_boots,'omitnan');
                end
            end

            if settings.compute_ci
                % Fix random seed
                rng('Default')
                rng(42)

                % Average across all wall distances for baseline
                nfr = numel(wallDist);
                fr_base_boots = zeros(nfr,settings.nboot);
                for n = 1:settings.nboot
                    frn = nfr;
                    frb = fr_baseline;
                    idx = datasample(1:nfr, frn, 'Replace',true);
                    fr_base_boots(:,n) = frb(idx);
                end
                fr_base_boot = mean(fr_base_boots,'omitnan');
                ci_base = prctile(fr_base_boot, settings.ci, 2);

                % Average each wall distance for wall period
                ci_fr = prctile(fr_boot, settings.ci, 2);
                ci_zfr = prctile(zfr_boot, settings.ci, 2);
                % ci_auc = prctile(auc_boot, settings.ci, 2);
            end
    
            if settings.show_plots || settings.save_plots
                if settings.plot_ci
                    err = ci_fr - shiftdim(mean_fr);
                    err_zfr = ci_zfr - shiftdim(mean_zfr);
                else
                    err = [-sem_fr; sem_fr]';
                    err_zfr = [-sem_zfr; sem_zfr]';
                end

                % Plot mean firing rate for each wall distance
                if c==1
                    figure(u); clf
                    figpos = [10 10 10 10];
                    set(gcf, 'Color','w', 'Renderer','painters', 'Units','Centimeters', 'Position',figpos)
                    if ~settings.show_plots
                        set(gcf, 'visible','off')
                    end
                    pos = get(gcf,'Position');
                    tl = tiledlayout(1,num_conditions, 'TileSpacing','compact', 'Padding','compact');
                    yl(1) = 0;
                    yl(2) = 0.001 + max([0 max(mean_fr + err(:,2))],[],'omitnan');
                    xlim([min(dist_centers_used)-1 max(dist_centers_used)+1])
                else
                    figure(u);
                    yl = ylim;
                    yl(1) = max([0 yl(1)],[],'omitnan');
                    yl(2) = max([yl(2) 0.001 + max(mean_fr + err(:,2))],[],'omitnan');
                    ylim(yl)
                    xlim([min(dist_centers_used)-1 max(dist_centers_used)+1])
                end
                nexttile
                hold on

                % Plot (smoothed) line for tuning curve
                plot(wd, mean_fr, 'k', 'LineWidth',1)

                % Plot dot +/- SEM for each point in "turbo" color scheme
                colors = flipud(turbo(numel(wd)));
                for d = 1:numel(wd)
                    % MPlot.ErrorShade(mean_fr, sem_fr, wd); hold on
                    plot(wd(d), mean_fr(d), '.', ...
                        'Color',colors(d,:), 'MarkerSize',15)
                    plot([wd(d) wd(d)], [mean_fr(d)+err(d,1) mean_fr(d)+err(d,2)], ...
                        'Color',colors(d,:), 'LineWidth',1)
                end

                if strcmpi(settings.dataset, 'laser')
                    xt = wd;
                    xl = [min(wd)-1 max(wd)+1];
                else
                    xt = dist_centers;
                    xl = [min(dist_centers)-1 max(dist_centers)+1];
                end

                set(gca, 'TickDir','out', 'box','off', 'XTick',xt, 'XTickLabelRotation',0)
                axis tight square
                title(gca, string(strjoin({session.mouseID,recID,'unit',num2str(u),'wall distance tuning'},' ')), 'FontSize',8)
                title(condition_names{c}, 'Units','normalized', 'Position',[0.5, 0.95, 0], 'FontWeight','normal')
                if c==1
                    ylabel('Firing rate (sp/s)')
                else
                    set(gca, 'YTickLabel',[])
                end
                xlabel('Wall distance (mm)')
                if isnan(yl(1))
                    yl(1) = 0;
                end
                if isnan(yl(2))
                    yl(2) = yl(1) + 0.001;
                end
                if yl(1) <= yl(2)
                    yl(2) = yl(2) + 0.001;
                end
                ylim(yl)
                xlim(xl)
                pause(0.001)
                
                % Plot z-scored firing rate for each wall distance
                if c==1
                    figure(100+u); clf
                    pos(1) = pos(1) + pos(3);
                    set(gcf, 'Color','w', 'Renderer','painters', 'Units','Centimeters', 'Position',pos)
                    tl = tiledlayout(1,num_conditions, 'TileSpacing','compact', 'Padding','compact');
                    yls = [min(min(mean_zfr' + err_zfr)) max(max(mean_zfr' + err_zfr))];
                else
                    figure(100+u);
                    yls = ylim;
                    yls(1) = min([yls(1) min(min(mean_zfr' + err_zfr))],[],'omitnan');
                    yls(2) = max([yls(2) max(max(mean_zfr' + err_zfr))],[],'omitnan');
                    ylim(yls)
                end
                b = nexttile;
                hold on
                % Plot (smoothed) line for tuning curve
                plot(wd, mean_zfr, 'k', 'LineWidth',1)

                % Plot dot +/- SEM for each point in "turbo" color scheme
                colors = flipud(turbo(numel(wd)));
                for d = 1:numel(wd)
                    % MPlot.ErrorShade(mean_fr, sem_fr, wd); hold on
                    plot(wd(d), mean_zfr(d), '.', ...
                        'Color',colors(d,:), 'MarkerSize',15)
                    plot([wd(d) wd(d)], [mean_zfr(d)+err_zfr(d,1) mean_zfr(d)+err_zfr(d,2)], ...
                        'Color',colors(d,:), 'LineWidth',1)
                end

                % MPlot.ErrorShade(mean_zfr, sem_zfr, wd); hold on
                % plot(wd, mean_zfr,'k')

                set(gca, 'TickDir','out', 'box','off', 'XTick',dist_centers, 'XTickLabelRotation',0)
                t = title(tl,strjoin({session.mouseID,recID,'unit',num2str(u),'z-scored wall distance tuning'},' '), ...
                    'FontSize',8);
                if c==1
                    ylabel('Z-Scored Firing rate')
                else
                    set(gca, 'YTickLabel',[])
                end
                xlabel('Wall distance (mm)')
                axis tight square
                title(b, condition_names{c}, 'Units','normalized', 'Position',[0.5, 0.95, 0])

                if isnan(yls(1))
                    yls(1) = 0;
                end
                if isnan(yls(2))
                    yls(2) = yls(1) + 0.001;
                end
                if yls(1) <= yls(2)
                    yls(2) = yls(2) + 0.001;
                end
                ylim(yls)
                xlim([min(dist_centers)-1 max(dist_centers)+1])
                pause(0.001)
            end
            
            if c==1
                % Save tuning data
                tuning.recID{cnt,1} = session.recID;
                tuning.unit(cnt,1) = u;
                tuning.distMM{cnt,1} = shiftdim(wd);
                tuning.fr{cnt,1} = shiftdim(fr_rep); % {distances}(n replicates)
                tuning.fr_base{cnt,1} = shiftdim(fr_base); 
                tuning.fr_base_mean{cnt,1} = shiftdim(mean_fr_base); 
                tuning.fr_base_sem{cnt,1} = shiftdim(sem_fr_base); 
                tuning.mean{cnt,1} = shiftdim(mean_fr);
                tuning.sem{cnt,1} = shiftdim(sem_fr);
                tuning.mean_shuf{cnt,1} = shiftdim(mean_fr_shuf);
                tuning.sem_shuf{cnt,1} = shiftdim(sem_fr_shuf);
                tuning.zfr{cnt,1} = shiftdim(zfr_rep); % {distances}(n replicates)
                tuning.mean_z{cnt,1} = shiftdim(mean_zfr);
                tuning.sem_z{cnt,1} = shiftdim(sem_zfr);
                tuning.mean_shuf_z{cnt,1} = shiftdim(mean_zfr_shuf);
                tuning.sem_shuf_z{cnt,1} = shiftdim(sem_zfr_shuf);
                % tuning.auc{cnt,1} = shiftdim(auc_rep); % {distances}(n replicates)
                % tuning.mean_auc{cnt,1} = shiftdim(mean_auc);
                % tuning.sem_auc{cnt,1} = shiftdim(sem_auc);
                % tuning.mean_shuf_auc{cnt,1} = shiftdim(mean_auc_shuf);
                % tuning.sem_shuf_auc{cnt,1} = shiftdim(sem_auc_shuf);
                
                if settings.compute_ci
                    tuning.ci_base{cnt,1} = shiftdim(ci_base);
                    tuning.ci{cnt,1} = shiftdim(ci_fr);
                    tuning.ci_z{cnt,1} = shiftdim(ci_zfr);
                    % tuning.ci_auc{cnt,1} = shiftdim(ci_auc);
                end
            else
                % Save tuning data for perturbation trials
                tuning.fr_condition{cnt,1} = shiftdim(fr_rep); % {distances}(n replicates)
                tuning.fr_base_condition{cnt,1} = shiftdim(fr_base); 
                tuning.mean_condition{cnt,1} = shiftdim(mean_fr);
                tuning.sem_condition{cnt,1} = shiftdim(sem_fr);
                tuning.mean_shuf_condition{cnt,1} = shiftdim(mean_fr_shuf);
                tuning.sem_shuf_condition{cnt,1} = shiftdim(sem_fr_shuf);
                tuning.mean_z_condition{cnt,1} = shiftdim(mean_zfr);
                tuning.sem_z_condition{cnt,1} = shiftdim(sem_zfr);
                tuning.mean_shuf_z_condition{cnt,1} = shiftdim(mean_zfr_shuf);
                tuning.sem_shuf_z_condition{cnt,1} = shiftdim(sem_zfr_shuf);
                % tuning.mean_auc_condition{cnt,1} = shiftdim(mean_auc);
                % tuning.sem_auc_condition{cnt,1} = shiftdim(sem_auc);
                % tuning.mean_shuf_auc_condition{cnt,1} = shiftdim(mean_auc_shuf);
                % tuning.sem_shuf_auc_condition{cnt,1} = shiftdim(sem_auc_shuf);

                if settings.compute_ci
                    tuning.ci_base_condition{cnt,1} = shiftdim(ci_base);
                    tuning.ci_condition{cnt,1} = shiftdim(ci_fr);
                    tuning.ci_z_condition{cnt,1} = shiftdim(ci_zfr);
                    % tuning.ci_auc_condition{cnt,1} = shiftdim(ci_auc);
                end
            end
            if c == num_conditions
                cnt = cnt + 1;
            end
        end

        % Save plots
        if settings.save_plots
            figure(u)
            filename = [rec_name '_unit' num2str(u) '_tuning'];
            filepath = fullfile(settings.fig_dir, 'Tuning Curves', filename);
            if ~isfolder(fullfile(settings.fig_dir, 'Tuning Curves'))
                mkdir(fullfile(settings.fig_dir, 'Tuning Curves'))
                pause(0.001)
            end
            % savefig(gcf, filepath)
            saveas(gcf, filepath, 'png')
            saveas(gcf, filepath, 'pdf')
            pause(0.001)

            figure(100+u)
            filename = [rec_name '_unit' num2str(u) '_zscore_tuning'];
            filepath = fullfile(settings.fig_dir, 'Tuning Curves', filename);
            % savefig(gcf, filepath)
            saveas(gcf, filepath, 'png')
            saveas(gcf, filepath, 'pdf')
            pause(0.001)
        end
    end
    toc
end

% Save tuning table and PSTH data
disp('Saving tuning curve datatable ...')
tt = struct2table(tuning);
dataTable = table(tt.recID, tt.unit, tt.distMM, ...
    tt.fr, tt.mean, tt.sem, ...
    tt.mean_shuf, tt.sem_shuf, ...
    tt.fr_base, tt.fr_base_mean, tt.fr_base_sem, ...
    tt.zfr, tt.mean_z, tt.sem_z, tt.mean_shuf_z, tt.sem_shuf_z, ...
    'VariableNames',{...
    'recID','unit','distanceMM', ...
    'fr_rep','fr_mean','fr_sem', ...
    'fr_mean_shuf','fr_sem_shuf', ...
    'fr_base','fr_base_mean','fr_base_sem' ...
    'zfr_rep','zfr_mean','zfr_sem','zfr_mean_shuf','zfr_sem_shuf'});
dataTable.Properties.VariableUnits = {...
    '','','mm', ...
    'sp/s','sp/s','sp/s', ...
    'sp/s','sp/s', ...
    'sp/s','sp/s','sp/s', ...
    '','','','','', ...
    };

% Append bootstrap data (if computed)
if settings.compute_ci
    dataTable.fr_ci = tt.ci; 
    dataTable.Properties.VariableUnits{end} = 'sp/s';
    dataTable.fr_base_ci = tt.ci_base; 
    dataTable.Properties.VariableUnits{end} = 'sp/s';
    dataTable.zfr_ci = tt.ci_z; 
    dataTable.Properties.VariableUnits{end} = '';
end

% Append condition data (if perturbation performed)
if isfield(session, 'trial_conditions')
    dataTable = addvars(dataTable, ...
        tt.fr_base_condition, tt.fr_condition, ...
        tt.mean_condition, tt.sem_condition, ...
        tt.mean_shuf_condition, tt.sem_shuf_condition, ...
        tt.mean_z_condition, tt.sem_z_condition, ...
        tt.mean_shuf_z_condition, tt.sem_shuf_z_condition, ...
        'NewVariableNames', ...
        {'fr_base_cond', 'fr_rep_cond', ...
        'fr_mean_cond', 'fr_sem_cond', ...
        'fr_mean_shuf_cond', 'fr_sem_shuf_cond', ...
        'zfr_mean_cond', 'zfr_sem_cond', ...
        'zfr_mean_shuf_cond', 'zfr_sem_shuf_cond', ...
        } ...
    );

    dataTable.Properties.VariableUnits(end-9:end) = ...
        {'sp/s','sp/s','sp/s','sp/s','sp/s','sp/s','','','',''}; 

    % Append bootstrap data (if computed)
    if settings.compute_ci
        dataTable.fr_ci_condition = tt.ci_condition; 
        dataTable.Properties.VariableUnits{end} = 'sp/s';
        dataTable.fr_base_ci_condition = tt.ci_base_condition; 
        dataTable.Properties.VariableUnits{end} = 'sp/s';
        dataTable.zfr_ci_condition = tt.ci_z_condition; 
        dataTable.Properties.VariableUnits{end} = '';
    end
end

clear tt

if settings.save_table
    filename = fullfile(settings.output_path,'tuningTable_PrV');
    disp(['Saving tuning data to ' filename ' ... '])
    save(filename,'dataTable')
end

disp('Done!')