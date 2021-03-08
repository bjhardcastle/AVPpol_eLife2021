% Generate plots for Fig2 (two MeTu drivers in AOTU)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig3path),'dir')
    try rmdir(fullfile(fig3path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Selectivity histograms and boxplots for pol vs no pol median values
pathsAVP
printpath = fig2path;
savename = 'MeTu_AOTU';
prefix = 'psi_';

% abs values
b = boxplotPSI({'R56F07_AOTU';'R73C04_AOTU'});
suffix = 'box';
printAVP

% mean-of-controls-subtracted values
b = boxplotPSI({'R56F07_AOTU';'R73C04_AOTU'},'pol','compare');
suffix = 'diff';
printAVP

pdfplotPSI({'R56F07_AOTU';'R73C04_AOTU'},'add','none');
suffix = 'pdf';
printAVP