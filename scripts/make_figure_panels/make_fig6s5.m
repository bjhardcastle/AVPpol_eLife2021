% Generate plots for Fig6s5 (TuBu/R4m polarotopy scatter plots, anterior bulb)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig6s5path),'dir')
    try rmdir(fullfile(fig6s5path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TuBu_a, anterior bulb

pathsAVP

lineStr = 'R34H10_Bu';
printpath = fig6s5path;

savePlotStr = {'horiz_all';'vert_all';'circ_all'};
plot_TuBu_a_R4m_EPG_ROIstruct_polarotopy 
% script requires lineStr (for ROIstruct), savePlotStr (selection to save)
% and printpath (for printAVP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% R4m, anterior bulb

pathsAVP

lineStr = 'R34D03_Bu';
printpath = fig6s5path;

savePlotStr = {'horiz_all';'vert_all';'circ_all'};
plot_TuBu_a_R4m_EPG_ROIstruct_polarotopy 
% script requires lineStr (for ROIstruct), savePlotStr (selection to save)
% and printpath (for printAVP)
