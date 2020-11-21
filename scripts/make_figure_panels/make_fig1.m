% Generate plots for Fig1 (DmDRA in MEDRA)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig1path),'dir')
    try rmdir(fullfile(fig1path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Timeseries from 3 ROIs in similar positions in 3 flies (exp2 is with polarizer removed)
loadDmDRA
printpath = fig1path;

prefix = 'polMapping_exp4_';
plotPolMapTrialManual( x([selectObj(1),selectObj(4),selectObj(6)]) , [1:3] , 0,[0 3])
suffix = '_manualROIs_R';
printAVP


prefix = 'polMapping_exp2_';
plotPolMapTrialManual( x([selectObj(3),selectObj(5),selectObj(7)]), [2:4] , 0,[0 3])
suffix = '_manualROIs_R';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example image showing corresponding ROIs for timeseries above
loadDmDRA
printpath = fig1path;

prefix = 'avgImg_exp4_';
plotCombPolImgManual( x(selectObj(6)).MIP ,0,[1:3])  % no cols
suffix = '_manualROIs_R';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSI maps for example recording above, including layer mask
loadDmDRA
printpath = fig1path;

prefix = 'polSelectivity_exp4_';
plotPolSelImg( x(selectObj(6)).MIP ,1,-1) 
suffix = 'R_MIP_mask';
printAVP


prefix = 'polSelectivity_exp2_';
plotPolSelImg( x(selectObj(7)).MIP ,1,-1) 
suffix = 'R_MIP_mask';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Selectivity histograms and boxplots for pol vs no pol median values
pathsAVP
printpath = fig1path;
savename = 'DmDRA';
prefix = 'psi_';

% abs values
b = boxplotPSI('DmDRA'); 
suffix = 'box';
printAVP

% abs values (Ctrl)
b = boxplotPSI('DmDRA','pol','ctrl'); 
suffix = 'box_ct';
printAVP

% mean-of-controls-subtracted values
b = boxplotPSI('DmDRA','pol','compare');
suffix = 'diff';
printAVP

% probability density
pdfplotPSI('DmDRA','add','ctrlCell');
% legend off
suffix = 'pdf';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RHS DmDRA polarization tuning maps 
loadDmDRA
printpath = fig1path;

prefix = 'polTuning_';

superUseMSP(x,1)
loadLayerMasks([x.MIP])
superPolThreshold(x,-1); 

plotCombPolImg( x(selectObj(1)).MIP ,1)  % with mask applied
suffix = 'R_MIP_mask';
printAVP


%{
plotCombPolImgManual( x(selectObj(1)).MIP ,1,[],1,1)  % no mask
suffix = 'R_MIP_nomask';
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ant-post position in medulla vs pol tuning scatter plot. All recordings. Both colormaps
loadDmDRA
printpath = fig1path;

prefix = 'scatter_horiz_';

superUseMSP(x,1)
superPolThreshold(x,-1);

superXY(x,1,[],[],[],0,[],[],2)
suffix = '_sel';
printAVP


superXY(x,1,[],[],[],1,[],[],2)
suffix = '_tuning';
printAVP




