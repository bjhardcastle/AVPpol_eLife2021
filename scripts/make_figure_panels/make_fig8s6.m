% Generate plots for Fig8s6 (E-PG responses in protocerebral bridge)
% Each cell can be run independently

pathsAVP
thisfigpath = fig8s6path;
if exist(thisfigpath,'dir')
    try rmdir(thisfigpath,'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Glomerulus pairing scheme
loadSS00096_PB
printpath = fig8s6path;

plot_PB_glom_pairing_crosscorr
prefix = 'glom_pairing_crosscorr';
savename = '';
suffix = '';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Autocorrelation peak shift from stimulus period
pathsAVP
printpath = fig8s6path;
plot_R4m_EPG_autocorr_peak_timeshift

prefix = 'autocorr_peak_timeshift';
savename = 'R4m_EPG';
suffix = '';
printAVP

