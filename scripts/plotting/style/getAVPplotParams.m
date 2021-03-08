% Some parameter values which are common across plots, for ease of making
% consistent changes 

%% Colors %%

lightGreyCol = [0.9 0.9 0.9];
darkGreyCol = [0.5 0.5 0.5];


% Colors (trying to fit with Omoto 2017 diagrams)
%
% 
% TuBuCols(1,:) = [109 148 207]./255;   % inferior bulb 
% TuBuCols(2,:) = [57 80 163]./255;     % superior bulb
% TuBuCols(3,:) = [119 104 169]./255;   % anterior bulb
% TuBuCols(4,:) = [0 158 173]./255;     % TuTu
% 
% MeTuCols(1,:) = [84 207 96]./255;     % lighter green 
% MeTuCols(2,:) = [51 153 52]./255;     % darker green 
% 

objCols.R49E09_AOTU = [109 148 207]./255;   % inferior bulb 
objCols.R88A06_AOTU = [57 80 163]./255;     % superior bulb
objCols.R34H10_AOTU = [119 104 169]./255;   % anterior bulb
objCols.R49E09_Bu = [109 148 207]./255;     % inferior bulb 
objCols.R88A06_Bu = [57 80 163]./255;       % superior bulb
objCols.R88A06_Bu_ant = [87 74 130]./255;   % anterior bulb (in 88A06)
objCols.R34H10_Bu = [119 104 169]./255;     % anterior bulb
objCols.R73C04_AOTU = [84 207 96]./255;     % lighter green 
objCols.R56F07_AOTU = [51 153 52]./255;     % darker green 
objCols.R17F12_AOTU = [0 158 173]./255;     % turquoise 
objCols.R7R8 = [190 50 160]./255;           % purple/pink
objCols.DmDRA = [165 60 190]./255;          % purple/lilac
objCols.R34D03_Bu = [217 82 25]./255;       % red/orange
objCols.R34D03_EB = [200 62 55]./255;       % red/orange
objCols.R19C08_Bu = [220 25 34]./255;       % red
objCols.R78B06_Bu = [192 27 72]./255;       % red/purple
objCols.SS00096_PB = [203 158 56]./255;     % yellow 

fig3s1ROIcols(1,:) = [32 96 33]./255;       % dark green
fig3s1ROIcols(2,:) = [51 153 52]./255;      % mid green
fig3s1ROIcols(3,:) = [159 223 160]./255;    % light green
fig3s1ROIcols(4,:) = [33 32 96]./255;       % dark blue
fig3s1ROIcols(5,:) = [65 64 191]./255;      % mid blue
fig3s1ROIcols(6,:) = [160 159 223]./255;    % light blue
%% Plot aesthetics ** 
defaultAxisHeight_cm = 2;
defaultImageHeight_cm = 2.5;
defaultAxisWidth_cm = 1.5;

axisTickLength = 0.05; % set ticklength  = axisTickLength/ axisDim(max)  in centimeters
axisLabelFontSize = 5;
defaultMarkerSize = 3;
defaultLineWidth = 0.5;
defaultAxesWidth = 0.3;