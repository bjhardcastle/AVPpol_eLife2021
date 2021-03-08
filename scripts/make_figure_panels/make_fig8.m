% Generate plots for fig8 (TuBu/R4m comparison in anterior bulb)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig8path),'dir')
    try rmdir(fullfile(fig8path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSI values in ring neuron drivers tested
objStr = {'R19C08_Bu';'R34D03_Bu'};

pathsAVP
printpath = fig8path;
savename = 'Ring';
prefix = 'psi_';


% abs values 
b = boxplotPSI(objStr,'dist','mask');
suffix = 'box';
printAVP

% probability density
pdfplotPSI(objStr,'add','none');
% legend off
suffix = 'pdf';
printAVP


%% PSI values in R4m alone, compared to controls
pathsAVP
printpath = fig8path;
savename = 'R4m';
prefix = 'psi_';

% mean-of-controls-subtracted values
b = boxplotPSI('R34D03_Bu','pol','compare');
suffix = 'diff';
printAVP

% probability density
pdfplotPSI('R34D03_Bu','add','ctrlCell');
% legend off
suffix = 'pdf';
printAVP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tuning map in TuBu_a in anterior bulb 
loadR34H10_Bu
printpath = fig8path;
superUseMSP(x,1)
superPolThreshold(x,-1); 

useTuningCols = 1;
useLayerMask = 1;

% PSI map
prefix = 'polSel';
plotPolSelImg( x(selectObj(1)).MIP ,1,-1) 
printAVP


% tuning map with mask applied
prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(1)).MIP ,useTuningCols,-1,~useLayerMask,1) 
printAVP

%{
% tuning map without mask applied
prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(1)).MIP ,useTuningCols,-1,0,1) 
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tuning map in R4m in anterior bulb 
loadR34D03_Bu
printpath = fig8path;
superUseMSP(x,1)
superPolThreshold(x,-1); 

useTuningCols = 1;
useLayerMask = 1;


prefix = 'polSel';
plotPolSelImg( x(selectObj(1)).MIP ,1,-1) 
printAVP


prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(1)).MIP ,useTuningCols,[],~useLayerMask,1) 
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalized response curves
% make_panel_TuBu_a_R4m_ROI_tuning_curve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tuning map in R2 in superior bulb (Gal4 #1 - individual layer)
loadR19C08_Bu
printpath = fig8path;
superUseMSP(x,1)
superPolThreshold(x,-1); 

useTuningCols = 1;
useLayerMask = 1;

%{
% PSI map - individual layer
prefix = 'polSel';
plotPolSelImg( x(selectObj(1)).Layers(selectLayer(1)) ,1,-1) 
suffix = '_R_Layer_mask';
printAVP
%}

% PSI map - MIP
prefix = 'polSel';
plotPolSelImg( x(selectObj(1)).MIP ,1,-1)
suffix = '_R_MIP_mask';
printAVP
   
%{
% tuning map with mask applied - individual layer
prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(1)).Layers(selectLayer(1)) ,useTuningCols,-1,~useLayerMask,1) 
suffix = '_R_Layer_mask';
printAVP
%}

% tuning map with mask applied -MIP
prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(1)).MIP ,useTuningCols,-1,~useLayerMask,1) 
suffix = '_R_MIP_mask';
printAVP

%{
% tuning map without mask applied
prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(2)).MIP ,useTuningCols,-1,0,1) 
printAVP
%}

