close all
clear all
clc
%need to cd to this folder: C:\Users\kssev\Desktop\Manuscript\Figure S5-single neuron decoder
cd('C:\Users\kssev\Desktop\Manuscript\Figure S5-single neuron decoder')
load('single_neuron_decoder_results_table.mat');

%% Plot example unit's confusion matrix w/ acc values in diag 
X = 7:2:23;
Y = X;
 cm = T.Decoder(1276,1).all.data.confusion;   % fig S5 A unit number 142; fig S5 B 1223
%cm = results(1276).all.data.confusion;

figure(222); clf
imagesc('XData',X, 'YData',Y, 'CData',cm);
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

% Add text labels only along the diagonal with a larger font size
for i = 1:length(X)
    val = cm(i,i);  % diagonal element: row=i, col=i
    text(X(i), Y(i), sprintf('%.2f', val), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color','k', 'FontSize',12);  % Increase font size here
end

%% Plot example unit accuracy
x = X; 
y = diag(T.Decoder(1223,1).all.data.confusion); 

shuffled_y = diag(T.Decoder(1223,1).all.shuffled.confusion);
shuffled_ci = T.Decoder(1223,1).all.shuffled.accuracy_ci_perstimulus;  % 9x2: [lower, upper]

figure;
scatter(x, y, 50, 'b', 'filled');  % scatter actual
hold on;

% plot shuffled points
plot(x, shuffled_y, 'ko', 'LineWidth',1.5);

% calculate asymmetrical error lengths
err_lower = shuffled_y - shuffled_ci(:,1);  
err_upper = shuffled_ci(:,2) - shuffled_y;  

% plot error bars manually
for i = 1:length(x)
    plot([x(i) x(i)], [shuffled_y(i)-err_lower(i), shuffled_y(i)+err_upper(i)], 'k', 'LineWidth',1.2);
    cap_width = 0.3;  
    plot([x(i)-cap_width x(i)+cap_width], [shuffled_y(i)-err_lower(i), shuffled_y(i)-err_lower(i)], 'k', 'LineWidth',1.2); % bottom cap
    plot([x(i)-cap_width x(i)+cap_width], [shuffled_y(i)+err_upper(i), shuffled_y(i)+err_upper(i)], 'k', 'LineWidth',1.2); % top cap
end

xlabel('Distance');
ylabel('Accuracy');

% Set ticks explicitly
set(gca, 'XTick', X);  % ensures ticks are exactly at 7,9,11,…,23
xlim([min(X)-1 max(X)+1])
ylim([0 1]);  % or adjust based on your data
