%% ----------------------------
%  Load all CSVs and stack
%  ----------------------------
% EITHER: point to a folder with the 4 CSVs:
clear; clc
dataDir = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\Whisker-Wall Tracking';   % <-- change me
files = dir(fullfile(dataDir, 'whisker_retraction*.csv'));

T = table();
for k = 1:numel(files)
    Tk = readtable(fullfile(files(k).folder, files(k).name), 'TextType','string');
    T = [T; Tk]; %#ok<AGROW>
end

%% Make types explicit
T.subject   = categorical(T.subject);
T.condition = categorical(T.condition, {'anesthetized','awake'});  % baseline first
% whiskerID is 1/2/3; keep numeric. If it is text, do: T.whiskerID = double(categorical(T.whiskerID));

%% (Optional) give human-friendly whisker labels in output/plots
% Define your mapping here if you want pretty labels in results/plots:
whiskerLabels = containers.Map([1 2 3], {'Greek','C1','C2'});  % <-- edit if needed

%% ----------------------------
%  Collapse to mouse×whisker×condition
%  (averaging over reps & distances)
%  ----------------------------
G = groupsummary(T, {'subject','whiskerID','condition'}, 'mean', 'deltaAngle');
W = unstack(G, 'mean_deltaAngle', 'condition');  % columns: anesthetized, awake

% Ensure both condition columns exist
if ~ismember('awake', W.Properties.VariableNames),        W.awake = NaN(height(W),1); end
if ~ismember('anesthetized', W.Properties.VariableNames), W.anesthetized = NaN(height(W),1); end

% Save the per-mouse means (useful “source data”)
writetable(W, fullfile(dataDir, 'mouse_whisker_condition_means.csv'));

%% ----------------------------
%  Paired TWO-SIDED tests per whisker
%  H1: awake ≠ anesthetized
%  ----------------------------
uw = unique(W.whiskerID);
results = table('Size',[0 8], ...
    'VariableTypes', {'double','string','double','double','double','double','double','double'}, ...
    'VariableNames', {'whiskerID','whiskerLabel','n_mice','meanDiff_awake_minus_anesth','t','df','p_two_sided','sig_Bonf'});

p_all = nan(numel(uw),1);

for i = 1:numel(uw)
    w = uw(i);
    Wi = W(W.whiskerID == w, :);
    Wi = rmmissing(Wi, 'DataVariables', {'awake','anesthetized'});  % require both conds per mouse

    if height(Wi) < 2
        mDiff = NaN; tstat = NaN; df = NaN; p2 = NaN;
        % ci = [NaN NaN]; dz = NaN;
    else
        diffs = Wi.awake - Wi.anesthetized; % negative => awake more retracted
        [~, p2, ci, stats] = ttest(Wi.awake, Wi.anesthetized);  % TWO-SIDED
        mDiff = mean(diffs);
        tstat = stats.tstat;
        df = stats.df;
        % dz = mDiff / std(diffs);  % Cohen's dz (paired)
    end

    wl = string(w);
    if isKey(whiskerLabels, double(w)), wl = whiskerLabels(double(w)); end

    results = [results; {double(w), wl, height(Wi), mDiff, tstat, df, p2, NaN}]; %#ok<AGROW>
    p_all(i) = p2;
end

% Holm–Bonferroni across the 3 whiskers (two-sided p-values)
% results.p_corrected = holm_bonferroni(p_all);
% results.p_corrected = mafdr(p_all, 'BHFDR',true); % one p-value 

alpha = 0.05/3
results.sig_Bonf = results.p_two_sided < alpha;
disp(results(:, {'whiskerLabel','p_two_sided','sig_Bonf'}))

% (Optional) if you captured ci and dz above, you can add columns like:
% results.CI_lo = CI_lo_vec; results.CI_hi = CI_hi_vec; results.cohens_dz = dz_vec;

% Save stats table
writetable(results, fullfile(dataDir, 'optionA_twoSided_awake_vs_anesth_stats_by_whisker.csv'));

disp(results)


%% ----------------------------
%  JUST TO OUTPUT THE AWAKE AND ANESTH ANGLES FOR TABLE S6
% Paired TWO-SIDED tests per whisker
%  H1: awake ≠ anesthetized
%  ----------------------------
uw = unique(W.whiskerID);

results = table('Size',[0 12], ...
    'VariableTypes', {'double','string','double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'whiskerID','whiskerLabel','n_mice', ...
                      'mean_awake','sem_awake','mean_anesth','sem_anesth', ...
                      'meanDiff_awake_minus_anesth','t','df','p_two_sided','sig_Bonf'});

p_all = nan(numel(uw),1);

for i = 1:numel(uw)
    w = uw(i);
    Wi = W(W.whiskerID == w, :);
    Wi = rmmissing(Wi, 'DataVariables', {'awake','anesthetized'});  % require both conds per mouse
    n  = height(Wi);

    if n < 2
        mAw = NaN; seAw = NaN; mAn = NaN; seAn = NaN;
        mDiff = NaN; tstat = NaN; df = NaN; p2 = NaN;
    else
        % condition summaries
        mAw = mean(Wi.awake, 'omitnan');
        seAw = std(Wi.awake, 'omitnan')/sqrt(n);
        mAn = mean(Wi.anesthetized, 'omitnan');
        seAn = std(Wi.anesthetized, 'omitnan')/sqrt(n);

        % paired test and paired-diff summary
        diffs = Wi.awake - Wi.anesthetized;
        [~, p2, ~, stats] = ttest(Wi.awake, Wi.anesthetized);  % TWO-SIDED
        mDiff = mean(diffs, 'omitnan');
        tstat = stats.tstat;
        df = stats.df;
    end

    wl = string(w);
    if isKey(whiskerLabels, double(w)), wl = whiskerLabels(double(w)); end

    results = [results; {double(w), wl, n, ...
                         mAw, seAw, mAn, seAn, ...
                         mDiff, tstat, df, p2, NaN}]; %#ok<AGROW>
    p_all(i) = p2;
end

% Bonferroni across whiskers
alpha = 0.05/3;
results.sig_Bonf = results.p_two_sided < alpha;

% Display concise summary
disp(results(:, {'whiskerLabel','mean_awake','sem_awake','mean_anesth','sem_anesth','meanDiff_awake_minus_anesth','p_two_sided','sig_Bonf'}))

% Save stats table
writetable(results, fullfile(dataDir, 'optionA_twoSided_awake_vs_anesth_stats_by_whisker.csv'));

% (Optional) also save per-mouse values + deltas for each whisker
perMouse = table();
for i = 1:numel(uw)
    w = uw(i);
    Wi = W(W.whiskerID == w, :);
    Wi = rmmissing(Wi, 'DataVariables', {'awake','anesthetized'});
    if isempty(Wi), continue; end
    wl = string(w);
    if isKey(whiskerLabels, double(w)), wl = whiskerLabels(double(w)); end
    Tpm = Wi(:, {'subject','whiskerID','awake','anesthetized'});
    Tpm.whiskerLabel = repmat(wl, height(Tpm), 1);
    Tpm.delta_awake_minus_anesth = Tpm.awake - Tpm.anesthetized;
    perMouse = [perMouse; Tpm]; %#ok<AGROW>
end
writetable(perMouse, fullfile(dataDir, 'perMouse_awake_anesth_and_delta_by_whisker.csv'));





%% ----------------------------
%  Quick diagnostic plots
%  ----------------------------
figure('Color','w'); tl = tiledlayout(1, numel(uw), 'Padding','compact','TileSpacing','compact');
for i = 1:numel(uw)
    nexttile(tl,i);
    w = uw(i);
    Wi = W(W.whiskerID == w, :);
    Wi = rmmissing(Wi, 'DataVariables', {'awake','anesthetized'});
    diffs = Wi.awake - Wi.anesthetized;

    % per-mouse paired dots
    plot([0 1], [Wi.anesthetized, Wi.awake]', '-o'); hold on
    yline(0,'--');
    titleStr = sprintf('Whisker %s', results.whiskerLabel(i));
    title(titleStr);
    xlim([-0.25 1.25]); set(gca,'XTick',[0 1],'XTickLabel',{'Anesth','Awake'});
    ylabel('\Delta\theta_{mid} (deg)');
end

%% ----------------------------
%  %% After you build W (mouse×whisker×condition means) and before/after the t-tests

% --- Mouse-level condition summaries (Mean ± SEM across mice) ---
% G already has per-mouse means; reuse it to summarize across mice
S = groupsummary(G, {'whiskerID','condition'}, {'mean','std','numel'}, 'mean_deltaAngle');
% Compute SEM across mice
S.SEM = S.std_mean_deltaAngle ./ sqrt(S.numel_mean_deltaAngle);

% Make wide tables for means and SEMs
Means = unstack(S(:, {'whiskerID','condition','mean_mean_deltaAngle'}), 'mean_mean_deltaAngle', 'condition'); % cols: anesthetized, awake
SEMs  = unstack(S(:, {'whiskerID','condition','SEM'}),                 'SEM',                 'condition');

% (Optional) add labels
Means.whiskerLabel = strings(height(Means),1);
for i = 1:height(Means)
    w = Means.whiskerID(i);
    if exist('whiskerLabels','var') && isKey(whiskerLabels, double(w))
        Means.whiskerLabel(i) = whiskerLabels(double(w));
    else
        Means.whiskerLabel(i) = string(w);
    end
end

% --- Attach to your 'results' table in matching order ---
results = sortrows(results, 'whiskerID');  % ensure order by ID
Means   = sortrows(Means,   'whiskerID');
SEMs    = sortrows(SEMs,    'whiskerID');

results.awake_mean  = Means.awake;
results.awake_sem   = SEMs.awake;
results.anesth_mean = Means.anesthetized;
results.anesth_sem  = SEMs.anesthetized;

% --- Pretty print (Mean ± SEM) per whisker ---
fprintf('\nMouse-level means (mean ± SEM, deg):\n');
for i = 1:height(Means)
    lbl = char(Means.whiskerLabel(i));
    am = Means.awake(i);  as = SEMs.awake(i);
    nm = Means.anesthetized(i); ns = SEMs.anesthetized(i);
    fprintf('  %s — Awake: %.2f ± %.2f; Anesth: %.2f ± %.2f\n', lbl, am, as, nm, ns);
end

% (Your Bonferroni flag remains valid for two-sided tests)
alpha = 0.05/3;
results.sig_Bonf = results.p_two_sided < alpha;

% Save & display
writetable(results, fullfile(dataDir, 'optionA_twoSided_awake_vs_anesth_stats_by_whisker.csv'));
disp(results)

%% Run stats for delta Kappa (curvature)

% ----------------------------
%  Load ONLY delta-kappa CSVs and stack
%  ----------------------------
clear; clc
dataDir = 'E:\PrV_Wall_Recordings\Analysis\IntermediateData\Whisker-Wall Tracking';

% If you saved per-subject CSVs like:
%   whisker_retraction_source_kappa_beh104.csv, etc.
files = dir(fullfile(dataDir, 'whisker_retraction_source_kappa_*.csv'));

T = table();
for k = 1:numel(files)
    Tk = readtable(fullfile(files(k).folder, files(k).name), 'TextType','string');
    T  = [T; Tk]; 
end

% Make types explicit
T.subject   = categorical(T.subject);
T.condition = categorical(T.condition, {'anesthetized','awake'});  % baseline first

% Optional pretty whisker labels
whiskerLabels = containers.Map([1 2 3], {'Greek','C1','C2'});

% ----------------------------
%  Collapse to mouse×whisker×condition
%  (averaging over reps & distances)
%  ----------------------------
valueVar = 'deltaKappa';   % <-- KEY CHANGE (was deltaAngle)

G = groupsummary(T, {'subject','whiskerID','condition'}, 'mean', valueVar);
W = unstack(G, ['mean_' valueVar], 'condition');  % columns: anesthetized, awake

% Ensure both condition columns exist
if ~ismember('awake', W.Properties.VariableNames),        W.awake = NaN(height(W),1); end
if ~ismember('anesthetized', W.Properties.VariableNames), W.anesthetized = NaN(height(W),1); end

% Save per-mouse means (source data)
writetable(W, fullfile(dataDir, 'mouse_whisker_condition_means_deltaKappa.csv'));

% ----------------------------
%  Paired TWO-SIDED tests per whisker
%  H1: awake ≠ anesthetized
%  ----------------------------
uw = unique(W.whiskerID);

results = table('Size',[0 12], ...
    'VariableTypes', {'double','string','double', ...
                      'double','double','double','double', ...
                      'double','double','double','double','double'}, ...
    'VariableNames', {'whiskerID','whiskerLabel','n_mice', ...
                      'mean_awake','sem_awake','mean_anesth','sem_anesth', ...
                      'meanDiff_awake_minus_anesth','t','df','p_two_sided','sig_Bonf'});

p_all = nan(numel(uw),1);

for i = 1:numel(uw)
    w  = uw(i);
    Wi = W(W.whiskerID == w, :);
    Wi = rmmissing(Wi, 'DataVariables', {'awake','anesthetized'});  % require both
    n  = height(Wi);

    if n < 2
        mAw = NaN; seAw = NaN; mAn = NaN; seAn = NaN;
        mDiff = NaN; tstat = NaN; df = NaN; p2 = NaN;
    else
        % condition summaries across mice
        mAw = mean(Wi.awake, 'omitnan');
        seAw = std(Wi.awake, 'omitnan')/sqrt(n);
        mAn = mean(Wi.anesthetized, 'omitnan');
        seAn = std(Wi.anesthetized, 'omitnan')/sqrt(n);

        % paired test
        diffs = Wi.awake - Wi.anesthetized;
        [~, p2, ~, stats] = ttest(Wi.awake, Wi.anesthetized);  % TWO-SIDED
        mDiff = mean(diffs, 'omitnan');
        tstat = stats.tstat;
        df    = stats.df;
    end

    wl = string(w);
    if exist('whiskerLabels','var') && isKey(whiskerLabels, double(w))
        wl = whiskerLabels(double(w));
    end

    results = [results; {double(w), wl, n, ...
                         mAw, seAw, mAn, seAn, ...
                         mDiff, tstat, df, p2, NaN}];
    p_all(i) = p2;
end

% Bonferroni across the 3 whiskers
alpha = 0.05/3;
results.sig_Bonf = results.p_two_sided < alpha;

disp(results(:, {'whiskerLabel','mean_awake','sem_awake','mean_anesth','sem_anesth', ...
                 'meanDiff_awake_minus_anesth','p_two_sided','sig_Bonf'}))

% Save
outStats = fullfile(dataDir, 'twoSided_awake_vs_anesth_stats_by_whisker_deltaKappa.csv');
writetable(results, outStats);
disp(['Saved: ' outStats])

% ----------------------------
%  Optional: per-mouse values + deltas for each whisker (Δκ)
%  ----------------------------

perMouse = table('Size',[0 6], ...
    'VariableTypes', {'categorical','double','double','double','string','double'}, ...
    'VariableNames', {'subject','whiskerID','awake','anesthetized','whiskerLabel','delta_awake_minus_anesth'});

for i = 1:numel(uw)
    w  = uw(i);
    Wi = W(W.whiskerID == w, :);
    Wi = rmmissing(Wi, 'DataVariables', {'awake','anesthetized'});
    if isempty(Wi), continue; end

    if isKey(whiskerLabels, double(w))
        wl = string(whiskerLabels(double(w)));
    else
        wl = string(w);
    end

    Tpm = Wi(:, {'subject','whiskerID','awake','anesthetized'});
    Tpm.subject = categorical(Tpm.subject);
    Tpm.whiskerLabel = repmat(wl, height(Tpm), 1);
    Tpm.delta_awake_minus_anesth = Tpm.awake - Tpm.anesthetized;

    perMouse = [perMouse; Tpm(:, perMouse.Properties.VariableNames)];
end


outPM = fullfile(dataDir, 'perMouse_awake_anesth_and_delta_by_whisker_deltaKappa.csv');
writetable(perMouse, outPM);
disp(['Saved: ' outPM])

% ----------------------------
%  Quick paired-dot diagnostic plots (Δκ)
%  ----------------------------
figure('Color','w'); tl = tiledlayout(1, numel(uw), 'Padding','compact','TileSpacing','compact');
for i = 1:numel(uw)
    nexttile(tl,i);
    w  = uw(i);
    Wi = W(W.whiskerID == w, :);
    Wi = rmmissing(Wi, 'DataVariables', {'awake','anesthetized'});

    plot([0 1], [Wi.anesthetized, Wi.awake]', '-o'); hold on
    yline(0,'--');
    title(sprintf('Whisker %s', results.whiskerLabel(i)));
    xlim([-0.25 1.25]); set(gca,'XTick',[0 1],'XTickLabel',{'Anesth','Awake'});
    ylabel('\Delta\kappa (mm^{-1})');
end


