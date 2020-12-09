% Generate plots for Fig4s4 (blue flash responses in TuBu drivers in AOTU & BU)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig4s1path),'dir')
    try rmdir(fullfile(fig4s1path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
make_panel_TuBu_flash_response