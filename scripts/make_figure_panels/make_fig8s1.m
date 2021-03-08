% Generate plots for fig8s1 (TuBu/R4m polarotopy scatter plots, anterior bulb)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig8s1path),'dir')
    try rmdir(fullfile(fig8s1path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TuBu_a, anterior bulb

pathsAVP

lineStr = 'R34H10_Bu';
printpath = fig8s1path;

savePlotStr = {'horiz_all';'vert_all';'circ_all'};
plot_TuBu_a_R4m_EPG_ROIstruct_polarotopy 
% script requires lineStr (for ROIstruct), savePlotStr (selection to save)
% and printpath (for printAVP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% R4m, anterior bulb

pathsAVP

lineStr = 'R34D03_Bu';
printpath = fig8s1path;

savePlotStr = {'horiz_all';'vert_all';'circ_all'};
plot_TuBu_a_R4m_EPG_ROIstruct_polarotopy 
% script requires lineStr (for ROIstruct), savePlotStr (selection to save)
% and printpath (for printAVP)
