%% Tuning curve shape analysis script

clear; clc

% Settings
dl = 'Y:';
base_dir = fullfile(dl, 'Kyle\Data\WallPassing\PrV_Wall_Recordings\Analysis\IntermediateData');
fig_dir = fullfile(dl, 'Kyle\Data\WallPassing\Figures\Tuning Curve Shape');

% Load tuning curves from naive and telc datasets
dataTable = matfile(fullfile(base_dir, 'Sim', 'tuningTable_PrV_stats_sim_wildtype.mat')).dataTable;
dataTable_telc = matfile(fullfile(base_dir, 'Sim', 'tuningTable_PrV_stats_sim_telc.mat')).dataTable;

% Format into tables
dataTable = struct2table(dataTable);
dataTable_telc = struct2table(dataTable_telc);

% Format rownames
clear rn
for u = size(dataTable,1):-1:1
    rn{u,1} = [dataTable.recID{u} '_' num2str(dataTable.unit(u))]; % dataTable.properties.RowNames
end
dataTable.Properties.RowNames = rn;

clear rn
for u = size(dataTable,1):-1:1
    rn{u,1} = [dataTable_telc.recID{u} '_' num2str(dataTable_telc.unit(u))]; % dataTable.properties.RowNames
end
dataTable_telc.Properties.RowNames = rn;

binCenters = 7:2:23;

% Get average tuning curves
condExp = dataTable.fr_mean;
condExp_telc = dataTable_telc.fr_mean;

% Compute tuning-curve shape metrics for naive and telc populations
dataTable = add_distance_tuning_metrics(dataTable, binCenters);
dataTable_telc = add_distance_tuning_metrics(dataTable_telc, binCenters);

%% Panel G (right): Plot FWHM CDF for Naive vs. TeLC

x_wt = dataTable.FWHM;
x_telc = dataTable_telc.FWHM;

mask_wt = get_wall_tuned_mask(dataTable, x_wt);
mask_telc = get_wall_tuned_mask(dataTable_telc, x_telc);

plot_group_cdf(x_wt(mask_wt), x_telc(mask_telc), ...
    'Tuning curve width (mm)', ...
    90211, ...
    'SIM_CDF_full_width_half_max_WT_vs_TeLC_sim', ...
    fig_dir, ...
    'XLim',[0 18]);

%% Panel H (right): Plot modulation depth CDF for Naive vs. TeLC

x_wt = dataTable.modDepth;
x_telc = dataTable_telc.modDepth;

mask_wt = get_wall_tuned_mask(dataTable, x_wt);
mask_telc = get_wall_tuned_mask(dataTable_telc, x_telc);

fig = plot_group_cdf( ...
    x_wt(mask_wt), ...
    x_telc(mask_telc), ...
    'Modulation depth (A.U.)', ...
    90212, ...
    'SIM_CDF_modulation_depth_WT_vs_TeLC', ...
    fig_dir, ...
    'XLim',[0 max([x_wt(mask_wt); x_telc(mask_telc)])]);

%% Panel I (right): Plot signed monotonicity index histogram for Naive vs. TeLC side-by-side

x_wt = dataTable.DI;
x_telc = dataTable_telc.DI;

mask_wt = get_wall_tuned_mask(dataTable, x_wt);
mask_telc = get_wall_tuned_mask(dataTable_telc, x_telc);

plot_group_cdf(x_wt(mask_wt), x_telc(mask_telc), ...
    'Signed monotonicity index', ...
    90210, ...
    'SIM_CDF_signed_monotonicity_index_WT_vs_TeLC_sim', ...
    fig_dir, ...
    'XLim',[-1 1]);

%% Functions

% Wrapper to compute all tuning metrics from tuning curve table
function dataTable = add_distance_tuning_metrics(dataTable, binCenters, varargin)
    % Add distance-tuning shape metrics to a table.

    p = inputParser;
    p.addParameter('TuningVar', 'fr_mean', @(x)ischar(x) || isstring(x));
    p.addParameter('BaselineVar', 'fr_base_mean', @(x)ischar(x) || isstring(x));
    p.addParameter('BaselineDefault', NaN, @(x)isnumeric(x) && isscalar(x));
    p.parse(varargin{:});

    tuningVar = char(p.Results.TuningVar);
    baselineVar = char(p.Results.BaselineVar);
    baselineDefault = p.Results.BaselineDefault;

    nUnits = height(dataTable);

    MI = nan(nUnits,1);
    DI = nan(nUnits,1);
    bestDist = nan(nUnits,1);
    peakRate = nan(nUnits,1);
    baselineRate = nan(nUnits,1);

    modDepth = nan(nUnits,1);      % NEW: gain metric
    respMin = nan(nUnits,1);       % optional helper
    respRange = nan(nUnits,1);     % optional helper

    FWHM = nan(nUnits,1);
    left50 = nan(nUnits,1);
    right50 = nan(nUnits,1);
    halfMaxRate = nan(nUnits,1);
    leftHW = nan(nUnits,1);
    rightHW = nan(nUnits,1);
    AI = nan(nUnits,1);
    nValidBins = nan(nUnits,1);
    edgeSide = strings(nUnits,1);

    hasBaselineVar = ismember(baselineVar, dataTable.Properties.VariableNames);

    for u = 1:nUnits
        try
            tc = get_table_curve(dataTable.(tuningVar), u);
            nValidBins(u) = sum(isfinite(tc));

            if hasBaselineVar
                b = get_table_scalar(dataTable.(baselineVar), u, baselineDefault);
            else
                b = baselineDefault;
            end
            baselineRate(u) = b;

            [MI(u), DI(u), bestDist(u), peakRate(u), edgeSide(u)] = ...
                compute_distance_monotonicity_index(tc, binCenters);

            [modDepth(u), respMin(u), respRange(u)] = ...
                compute_distance_modulation_depth(tc, b);

            [FWHM(u), left50(u), right50(u), halfMaxRate(u)] = ...
                compute_distance_fwhm(tc, binCenters, b);

            [AI(u), leftHW(u), rightHW(u)] = ...
                compute_distance_asymmetry_index(bestDist(u), left50(u), right50(u));

        catch ME
            fprintf('Unit %d failed in add_distance_tuning_metrics: %s\n', u, ME.message);
        end
    end

    % Append metrics to table
    dataTable.MI = MI;
    dataTable.DI = DI;
    dataTable.bestDist = bestDist;
    dataTable.peakRate = peakRate;
    dataTable.edgeSide = edgeSide;
    dataTable.baselineRate = baselineRate;

    dataTable.modDepth = modDepth;     % NEW
    dataTable.respMin = respMin;       % optional
    dataTable.respRange = respRange;   % optional

    dataTable.FWHM = FWHM;
    dataTable.left50 = left50;
    dataTable.right50 = right50;
    dataTable.halfMaxRate = halfMaxRate;
    dataTable.leftHW = leftHW;
    dataTable.rightHW = rightHW;
    dataTable.AI = AI;
    dataTable.nValidBins = nValidBins;
end


% Tuning metrics
function [MI, DI, bestDist, peakRate, edgeSide] = compute_distance_monotonicity_index(tc, binCenters)
    % Compute simple distance-tuning monotonicity metrics.
%
% MI = unsigned edge-vs-peak metric in [0,1]
% DI = signed edge difference metric in [-1,1]
%
% INPUTS
%   tc         : vector of mean firing rates across ordered distance bins
%   binCenters : vector of distance bin centers, same length as tc
%
% OUTPUTS
%   MI         : monotonicity index in [0,1]
%                1 = edge-peaked / monotonic-like
%                0 = interior-peaked / map-like
%   DI         : signed monotonicity index in [-1,1]
%                <0 = near-preferring
%                >0 = far-preferring
%                 0 = symmetric edge responses or centered peak
%   bestDist   : distance at peak firing rate
%   peakRate   : peak firing rate
%   edgeSide   : which edge has larger response ("near" or "far")

    tc = tc(:)';
    binCenters = binCenters(:)';

    valid = isfinite(tc) & isfinite(binCenters);
    tc = tc(valid);
    binCenters = binCenters(valid);

    MI = NaN;
    DI = NaN;
    bestDist = NaN;
    peakRate = NaN;
    edgeSide = "";

    if isempty(tc)
        return
    end

    [peakRate, idxPeak] = max(tc);
    bestDist = binCenters(idxPeak);

    if tc(1) >= tc(end)
        edgeSide = "near";
    else
        edgeSide = "far";
    end

    % Min-subtracted curve so the metric reflects shape rather than offset
    tc0 = tc - min(tc);
    denom = max(tc0);

    if denom <= 0
        return
    end

    edgeResp = max(tc0(1), tc0(end));
    MI = edgeResp / denom;
    DI = (tc(end) - tc(1)) / denom;

    MI = max(0, min(1, MI));
    DI = max(-1, min(1, DI));
end

function [FWHM, left50, right50, halfMaxRate] = compute_distance_fwhm(tc, binCenters, baselineRate)
% Compute baseline-referenced full width at half maximum (FWHM) for a
% distance tuning curve using linear interpolation.
%
% Half-max is defined as:
%   baselineRate + 0.5 * (peakRate - baselineRate)
%
% FWHM is the distance span over which tc >= halfMaxRate.
%
% If the curve does not drop below half-max before one sampled edge,
% the width is truncated to that sampled edge.

    tc = tc(:)';
    binCenters = binCenters(:)';

    valid = isfinite(tc) & isfinite(binCenters);
    tc = tc(valid);
    binCenters = binCenters(valid);

    FWHM = NaN;
    left50 = NaN;
    right50 = NaN;
    halfMaxRate = NaN;

    if isempty(tc) || ~isfinite(baselineRate)
        return
    end

    [peakRate, idxPeak] = max(tc);
    peakMod = peakRate - baselineRate;

    if ~isfinite(peakMod) || peakMod <= 0
        return
    end

    halfMaxRate = baselineRate + 0.5 * peakMod;

    % Left crossing
    leftVals = tc(1:idxPeak);
    leftBins = binCenters(1:idxPeak);

    idxBelowLeft = find(leftVals < halfMaxRate, 1, 'last');

    if isempty(idxBelowLeft)
        left50 = leftBins(1); % truncated at left sampled edge
    elseif idxBelowLeft == numel(leftVals)
        left50 = NaN;
    else
        x1 = leftBins(idxBelowLeft);
        x2 = leftBins(idxBelowLeft+1);
        y1 = leftVals(idxBelowLeft);
        y2 = leftVals(idxBelowLeft+1);
        left50 = interp_crossing(x1, y1, x2, y2, halfMaxRate);
    end

    % Right crossing
    rightVals = tc(idxPeak:end);
    rightBins = binCenters(idxPeak:end);

    idxBelowRight = find(rightVals < halfMaxRate, 1, 'first');

    if isempty(idxBelowRight)
        right50 = rightBins(end); % truncated at right sampled edge
    elseif idxBelowRight == 1
        right50 = NaN;
    else
        x1 = rightBins(idxBelowRight-1);
        x2 = rightBins(idxBelowRight);
        y1 = rightVals(idxBelowRight-1);
        y2 = rightVals(idxBelowRight);
        right50 = interp_crossing(x1, y1, x2, y2, halfMaxRate);
    end

    if isfinite(left50) && isfinite(right50)
        FWHM = right50 - left50;
    end
end

function [AI, leftHW, rightHW] = compute_distance_asymmetry_index(bestDist, left50, right50)
    % Compute asymmetry index from half-max widths.
    %
    % AI = (leftHW - rightHW) / (leftHW + rightHW)
    %
    % OUTPUTS
    %   AI      : asymmetry index in [-1,1]
    %             <0 broader on near/left side
    %             >0 broader on far/right side
    %   leftHW  : half-width on near/left side
    %   rightHW : half-width on far/right side

    AI = NaN;
    leftHW = NaN;
    rightHW = NaN;

    if ~isfinite(bestDist) || ~isfinite(left50) || ~isfinite(right50)
        return
    end

    leftHW = bestDist - left50;
    rightHW = right50 - bestDist;

    if ~isfinite(leftHW) || ~isfinite(rightHW) || leftHW < 0 || rightHW < 0
        return
    end

    denom = leftHW + rightHW;
    if denom <= 0
        return
    end

    AI = (leftHW - rightHW) / denom;
    AI = max(-1, min(1, AI));
end

function [modDepth, respMin, respRange] = compute_distance_modulation_depth(tc, baselineRate)
    % Compute tuning-curve modulation depth relative to baseline.
    %
    % ModDepth = peak response - baseline
    %
    % OUTPUTS
    %   modDepth : gain above baseline
    %   respMin  : minimum firing rate in the tuning curve
    %   respRange: peak-to-trough range across the tuning curve

    tc = tc(:)';
    tc = tc(isfinite(tc));

    modDepth = NaN;
    respMin = NaN;
    respRange = NaN;

    if isempty(tc) || ~isfinite(baselineRate)
        return
    end

    peakRate = max(tc);
    respMin = min(tc);
    respRange = peakRate - respMin;
    modDepth = peakRate - baselineRate;
end


% Helper functions
function xCross = interp_crossing(x1, y1, x2, y2, yCross)
% Linear interpolation for x at which y reaches yCross between two points.

    if ~isfinite(x1) || ~isfinite(x2) || ~isfinite(y1) || ~isfinite(y2) || ~isfinite(yCross)
        xCross = NaN;
        return
    end

    if y2 == y1
        xCross = mean([x1 x2]);
        return
    end

    xCross = x1 + (yCross - y1) * (x2 - x1) / (y2 - y1);
end

function tc = get_table_curve(col, u)
    % Robustly extract one tuning curve from a table variable.
    % Supports cell arrays, numeric matrices, and row vectors.

    if iscell(col)
        tc = col{u};
    elseif isnumeric(col) || islogical(col)
        if isvector(col)
            tc = col(u);
        else
            tc = col(u,:);
        end
    else
        error('Unsupported tuning-curve variable type: %s', class(col));
    end

    tc = double(tc(:)');
end

function x = get_table_scalar(col, u, defaultValue)
    % Robustly extract one scalar baseline from a table variable.
    % If the entry is a vector, returns its mean omitting NaNs.

    if nargin < 3
        defaultValue = NaN;
    end

    if iscell(col)
        val = col{u};
    elseif isnumeric(col) || islogical(col)
        val = col(u,:);
    else
        x = defaultValue;
        return
    end

    if isempty(val)
        x = defaultValue;
        return
    end

    val = double(val(:));
    val = val(isfinite(val));

    if isempty(val)
        x = defaultValue;
    else
        x = mean(val, 'omitnan');
    end
end

function mask = get_wall_tuned_mask(tbl, x)
    % Return analysis mask for wall-tuned, responsive, non-suppressed,
    % acceptable, non-low-firing units with finite metric x.

    mask = ismember(tbl.map,{'proximity','map'}) & isfinite(x);
end


% Plotting Helpers
function Mnorm = normalize_rows_minmax(M)
    Mnorm = nan(size(M));
    for r = 1:size(M,1)
        x = M(r,:);
        if all(~isfinite(x))
            continue
        end
        xmin = min(x, [], 'omitnan');
        xmax = max(x, [], 'omitnan');
        if isfinite(xmin) && isfinite(xmax) && xmax > xmin
            Mnorm(r,:) = (x - xmin) ./ (xmax - xmin);
        else
            Mnorm(r,:) = zeros(size(x));
        end
    end
end

% Helper for consistent two-group CDF panels
function fig = plot_group_cdf(x_wt, x_telc, xlab, figNum, filename, fig_dir, varargin)

    p = inputParser;
    p.addParameter('XLim',[]);
    p.addParameter('YLim',[0 1]);
    p.addParameter('XTick',[]);
    p.addParameter('YTick',0:0.2:1);
    p.addParameter('WTColor','k');
    p.addParameter('TeLCColor',[0 0.6 0]);
    p.parse(varargin{:});

    x_wt = x_wt(isfinite(x_wt));
    x_telc = x_telc(isfinite(x_telc));
    
    % Format figure
    fig = figure(figNum); clf
    Paperize(gcf, 'Position',[2 2 5 5])
    hold on; clear h

    % WT
    if ~isempty(x_wt)
        [f, x] = ecdf(x_wt);
        h(1) = plot(x, f, '-', 'Color',p.Results.WTColor, 'LineWidth',1);
    end

    % TeLC
    if ~isempty(x_telc)
        [f, x] = ecdf(x_telc);
        h(2) = plot(x, f, '-', 'Color',p.Results.TeLCColor, 'LineWidth',1);
    end

    % Insert legend
    legend(h, {['WT (n=' num2str(numel(x_wt)) ')'],['TeLC (n=' num2str(numel(x_telc)) ')']}, 'Location','southeast', 'Box','off')

    % Format axes
    xlim(p.Results.XLim); ylim(p.Results.YLim);
    set(gca, 'YTick',p.Results.YTick)
    if ~isempty(p.Results.XTick), set(gca, 'XTick',p.Results.XTick); end
    axis square
    xlabel(xlab)
    ylabel('Cum. frac')
    FormatAxes(gca);
    
    % Save figure
    saveas(gcf, fullfile(fig_dir, filename), 'pdf')
    saveas(gcf, fullfile(fig_dir, filename), 'png')
end

function Paperize(fig, varargin)
    % "Paperize" figure formatting to standard settings

    % Handle varargin (user inputs)
    p = inputParser();
    p.addParameter("Position",[5 5 10 10]);
    p.addParameter("Color","w")
    p.addParameter("Renderer","painters")
    p.parse(varargin{:});
    set(fig, ...
        "Color",p.Results.Color, ...
        "Units","centimeters", ...
        "Renderer",p.Results.Renderer, ...
        "Position",p.Results.Position...
        )
end

function FormatAxes(ax, varargin)
    % Formats axes to standard (or user-defined) settings

    % Handle varargin (user inputs)
    p = inputParser();
    p.addParameter("Color","k")
    p.addParameter("LineWidth",0.5)
    p.addParameter("FontName","arial")
    p.addParameter("FontSize",8)
    p.parse(varargin{:});

    set(ax, ...
        "YColor",p.Results.Color, ...
        "XColor",p.Results.Color, ...
        "LineWidth",p.Results.LineWidth, ...
        "FontName",p.Results.FontName, ...
        "FontSize",p.Results.FontSize, ...
        "TickDir","out", ...
        "Box","off")

    ax.XAxis.TickLabelColor = p.Results.Color;
    ax.YAxis.TickLabelColor = p.Results.Color;
    ax.XAxis.Label.Color = p.Results.Color;
    ax.YAxis.Label.Color = p.Results.Color;
end

