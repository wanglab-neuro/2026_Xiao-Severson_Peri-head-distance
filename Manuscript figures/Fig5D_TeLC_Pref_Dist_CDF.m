%% Plot CDFs to compare average firing rates for TeLC vs. Naive
% Compare baseline/wall firing between naive and telc populations

% We use CDFs and a two-sample statistical test (K-S test) here 
% to compare distributions because data from Control (wild-type/GFP) and TeLC 
% must be collected in separate animals due to chronic, irreversible manipulation 
% (tetanus toxin chemogenetic silencing)

% tl;dr 
% in TeLC condition, baseline and wall pass firing rates were
% signifcantly elevated for untuned cells but not for distance-tuned cells

clear; clc; close all
settings.region = 'PrV';
settings.fig_dir = 'E:\PrV_Wall_Recordings\Figures\TeLC\CDFs';
settings.save_plots = true; % false; % 

% Load naive psth data
settings.wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_naive';
sst_naive = matfile(fullfile(settings.wv_stats, 'PrV_spike_stats_table_naive.mat')).sst;

data_folder = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData';
dataTable_naive = matfile(fullfile(data_folder,'tuningTable_PrV_stats.mat')).dataTable;

% Load telc psth data
settings.dataset = 'telc';
settings.wv_stats = 'E:\PrV_Wall_Recordings\WaveformStats_TeLC';
sst_telc = matfile(fullfile(settings.wv_stats, 'PrV_spike_stats_table_telc.mat')).sst;

% Load tuning curve data combined from naive and whisker cutting experiments
data_folder = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\TeLC';
dataTable_telc = matfile(fullfile(data_folder,'tuningTable_PrV_stats.mat')).dataTable;

% Initialize stats table to save for easy reference in text
variableNames = {'G1', 'G2', 'Data G1 (Sp/s)', 'Mean G1 (Sp/s)', 'CI G1 (Sp/s)', 'Data G2 (Sp/s)', 'Mean G2 (Sp/s)', 'CI G2 (Sp/s)', 'p-value'};
rowNames = {
    'Baseline_Tuned_Naive_vs_TeLC', ... % one-tailed
    'Wall_Tuned_Naive_vs_TeLC', ... % one-tailed
    'Baseline_Untuned_Naive_vs_TeLC', ... % one-tailed
    'Wall_Untuned_Naive_vs_TeLC', ... % one-tailed
    'TeLC_Untuned_Baseline_vs_Wall', ... % one-tailed
    'TeLC_Tuned_Baseline_vs_Wall', ... % one-tailed
    'TeLC_Baseline_Untuned_vs_Tuned', ... % **two-tailed
    'TeLC_Wall_Untuned_vs_Tuned' ... % one-tailed
    };
tails = { ...
    'larger', ... % expect TeLC > Naive (disinhibited)
    'larger', ... % expect TeLC > Naive (disinhibited)
    'larger', ... % expect TeLC > Naive (disinhibited)
    'larger', ... % expect TeLC > Naive (disinhibited)
    'larger', ... % Wall > baseline
    'larger', ... % Wall > baseline
    'unequal', ... % No expectation
    'larger' ... % Wall > baseline
    };

% Initialize source data table
T = table('RowNames',rowNames); 
[T.G1, T.G2, T.data1] = deal(cell(numel(rowNames),1));
[T.mean1, T.ci1] = deal(NaN(numel(rowNames),1));
T.data2 = cell(numel(rowNames),1);
[T.mean2, T.ci2, T.p] = deal(NaN(numel(rowNames),1));

%% Panel D: Plot CDF comparing distribution of preferred wall distances Control vs. TeLC

dist_edges = 6:2:24;
dist_centers = 7:2:23;
[pref_distance, pref_distance_binned] = deal(cell(2,1));

% Get preferred wall distances for Control wall tuned and activated
sst1 = sst_naive;
dt1 = dataTable_naive;
mask1 = sst1.acceptable & ~dt1.low_firing & ~dt1.stats_err & dt1.activated; % sst1.mean_fr>1
[pref_distance{1}, pref_distance_binned{1}] = get_preferred_distance(mask1, dt1, dist_edges, dist_centers);

% Get preferred wall distances for TeLC wall tuned and activated
sst2 = sst_telc;
dt2 = dataTable_telc;
mask2 = sst2.acceptable & ~dt2.low_firing & ~dt2.stats_err & dt2.activated;
[pref_distance{2}, pref_distance_binned{2}] = get_preferred_distance(mask2, dt2, dist_edges, dist_centers);

% Plot CDFs and run stats for Preferred Distance 
figname = 'Pref_Dist_Naive_vs_TeLC';
name1 = 'Control';
name2 = 'TeLC';
data1 = pref_distance_binned{1}; % pref_distance{1}; % 
data2 = pref_distance_binned{2}; % pref_distance{2}; % 
color1 = [0 0 0]; % Black Naive
color2 = [0.5 0.85 0.5]; % Green TeLC

% Make figure for cdf comparing two datasets
figure(randi(1000)); clf;
set(gcf, 'Color','w', 'Renderer','painters', 'Units','centimeters', 'Position',[2 2 5 5])
hold on
h1 = cdfplot(data1);
h2 = cdfplot(data2);
title('')
set(h1, 'Color',color1, 'LineWidth',0.5);
set(h2, 'Color',color2, 'LineWidth',0.5);
MPlot.FormatAxes(gca)
% set(gca, 'XScale','log')
axis square
grid off
xlim([dist_edges(1) dist_edges(end)])
set(gca, 'XTick',dist_centers, 'XTickLabelRotation',0)
ylabel('Frac. Units')
xlabel('Pref. Distance (mm)')
% title(strrep(['CDF ' figname],'_',' '))
legend({name1,name2}, 'Location','southeast', 'Box','off')

if settings.save_plots
    filename = fullfile(settings.fig_dir, ['cdf_' figname]);
    saveas(gcf,filename,'pdf');
    exportgraphics(gcf,[filename '.png'], 'Resolution',600);
    pause(0.01)
end

% Run two-sample Kolmogorov-Smirnov (KS) Test for significance
[~,p,~] = kstest2(data1, data2, 'Tail','smaller');
disp(['P-value from KS-test: ' num2str(p)])
disp(newline)

%% Save source/results datatable

% Format table
T.Comparison = T.Properties.RowNames;
T = movevars(T, {'Comparison'}, 'Before','G1'); % add rownames as a variable
T = removevars(T, {'data1', 'data2'});

% Save table
filename = 'cdf_telc_source_data_results_table';
writetable(T, fullfile(settings.fig_dir, [filename '.csv']))

%% Functions

function [pref_distance, pref_distance_binned] = get_preferred_distance(mask, dataTable, dist_edges, dist_centers)
    % Initialize
    num_units = sum(mask);
    unit2use = find(mask);
    [pref_distance, pref_distance_binned] = deal(NaN(num_units,1));

    for u = 1:num_units
        uu = unit2use(u);
        wd = dataTable.distanceMM{uu};
        fr = dataTable.fr_mean{uu};
    
        % Bin wall distances
        [~,~,wd_bins] = histcounts(wd, dist_edges);
        actual_dist_centers = wd(wd_bins);
        fr2 = NaN(size(dist_centers));
        wd_unique = unique(wd_bins);
        for w = 1:numel(wd_unique)
            fr2(wd_bins==wd_unique(w)) = mean(fr(wd_bins==wd_unique(w)));
        end
        if strcmp(dataTable.map{uu},'suppressed')
            [~,pd] = min(fr2,[],2,'omitnan');
        else
            [~,pd] = max(fr2,[],2,'omitnan');
        end
        pref_distance(u) = actual_dist_centers(pd);
        pref_distance_binned(u) = dist_centers(pd);
    end
end

