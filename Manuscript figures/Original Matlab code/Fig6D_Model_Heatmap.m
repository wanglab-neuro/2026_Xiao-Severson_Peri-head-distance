%% Panel D: Quick tuning heatmap from modeled neurons CSV (rows = neurons, cols = activity bins)
cd('C:\Code\PrV_simulation\results_20251205_171134')
csv_files = {'heatmap.csv', 'heatmap_telc.csv'};
for f = 1:numel(csv_files)
    csv_file = csv_files{f};
    xvals      = [0:98];             % leave empty to use 1:nBins; or set to e.g., [6 8 10 ...]
    upsample   = 10;              % set >1 to interpolate columns (e.g., 20)
    save_plots = 1;
    save_base  = csv_file(1:end-4);
    
    % Load
    M = readmatrix(csv_file);           % rows = neurons, cols = bins
    valid_rows = any(~isnan(M),2);      % drop all-NaN rows defensively
    M = M(valid_rows,:);
    [nUnits, nBins] = size(M);
    
    % x-axis values
    if isempty(xvals), xvals = 1:nBins; end
    assert(numel(xvals)==nBins, 'xvals length must match number of columns in CSV');
    
    % Per-neuron normalization to [0,1]
    minv   = min(M,[],2,'omitnan');
    maxv   = max(M,[],2,'omitnan');
    rangev = maxv - minv;
    rangev(rangev==0) = 1; % avoid div-by-zero for flat rows
    Mnorm  = (M - minv) ./ rangev;
    
    % Optional column upsampling (purely cosmetic)
    if upsample > 1
        xi = linspace(1, nBins, nBins*upsample);
        Mup = nan(size(Mnorm,1), numel(xi));
        for i = 1:size(Mnorm,1)
            Mup(i,:) = interp1(1:nBins, Mnorm(i,:), xi, 'linear', 'extrap');
        end
        Mnorm = Mup;
        xvals = linspace(xvals(1), xvals(end), size(Mnorm,2));
    end
    
    % Sort neurons by preferred bin (argmax)
    [~, pref_idx] = max(Mnorm,[],2,'omitnan');
    [~, ord]      = sort(pref_idx,'ascend');
    H             = Mnorm(ord,:);
    bin_edges = linspace(0,max(xvals)+1,10);
    bin_centers = mean([bin_edges(1:end-1); bin_edges(2:end)]);
    
    % Plot
    figure(1); clf
    set(gcf, 'Color','w', 'Units','centimeters', 'Position',[2 2 6 16])
    imAlpha = ~isnan(H);
    imagesc(xvals, 1:size(H,1), H, 'AlphaData', imAlpha);
    colormap('parula'); axis tight
    set(gca, 'YDir','reverse', 'TickDir','out', 'Box','off', 'Color',[1 1 1], ...
        'XTick',bin_centers, 'XTickLabels',{'0','1','2','3','4','5','6','7','8'})
    cb = colorbar; set(cb,'TickDir','out'); ylabel(cb, 'Normalized activity', 'Rotation', 270)
    xlabel('Distance bin')  % change label if xvals are physical units
    ylabel('Neurons (far \rightarrow close)')
    title('Population tuning heat map')

    pause(0.001)
    
    if save_plots
        savefig(gcf, [save_base '.fig']);
        saveas(gcf, [save_base '.png']);
        saveas(gcf, [save_base '.pdf']);
    end

    pause(0.001)
end