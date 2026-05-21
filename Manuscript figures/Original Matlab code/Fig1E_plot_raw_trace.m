%% Load ephys data
% 
clear; clc
[recInfo,recordings,spikes,...
        videoTTL,trialTTL] = LoadData('20230208_sp_ephys2_session1.ns6','F:\PrV_Wall_Recordings\20230208_sp\session_1\ephys',[]);

% Channel num
channel_num = 1;

% 

% wall_onset = wall_onset(wall_onset>startPoint &  wall_onset<endPoint);
% wall_onset = wall_onset - startPoint;

% wall_offset = TTLtable.Var2;
% wall_offset = wall_offset(wall_offset>startPoint &  wall_offset<endPoint);
% wall_offset = wall_offset - startPoint;

%%
figure(1); clf
set(gcf, 'Color','w')
plot(x, filtered_trace);
title('High-Pass Filtered Ephys Trace');
xlabel('Time (s)');
ylabel('Amplitude');
xtick = gca().XTick;
set(gca, 'TickDir','out', 'box','off', 'XTickLabels',xtick+startPoint)
hold on
