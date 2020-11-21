%% Make all figures
% All plots with numerical data will be generated as pdf files, close to
% their appearance in the multi-panel figures. Object arrays for each Gal4
% will be loaded when needed - some are >1GB and can take a few seconds to
% open.

% some stats are printed in command window: keep a log file 
clear all 
pathsAVP
if ~exist(plotpathAVP,'file')
    mkdir(plotpathAVP)
end
delete(fig_command_window_output)
diary(fig_command_window_output)
diary on 

make_all_timer = tic;

make_fig1
make_fig1s1
make_fig1s2
make_fig2
make_fig2s3
make_fig3 
make_fig4
make_fig4s4
make_fig5
make_fig6
make_fig6s5
make_fig7
make_fig8
make_fig8s6

toc(make_all_timer)
diary off 

% when finished clear the last object array loaded 
clear x 