% Generate plots for Fig3 (pol+nonpol responses single MeTu driver in AOTU, TuTu characterization)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig3path),'dir')
    try rmdir(fullfile(fig3path),'s'),end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Three example left side recordings, showing retinotopy in non-pol responses
% ROI(1:3) lateral pol ROIs, ventral to dorsal
% ROI(4:6) mid bar ROIs, ventral to dorsal
% ROI(7) dorsal ROI inhibited by everything
% ROI(8) dorsal ROI, thrust resp

loadR73C04_AOTU
 
printpath = fig3path;


superUseMSP(x); 

% pol exp responses
prefix = 'polMapping_exp4_';
plotPolMapTrialManual(x(selectPBTobjects_L),[1:3],0)
suffix = '_manualROIs_L';
printAVP

prefix = 'polPolar_';
plotSinglePolarResp(x(selectPBTobjects_L),[1:3],0,0)
suffix = '_manualROIs_L';
printAVP



% non-pol ROI responses to non-pol stimuli:
prefix = 'barMapping_exp5_';
plotBarMapTrial(x(selectPBTobjects_L),[4:6],0)
suffix = '_manualROIs_L';
printAVP

%{
prefix = 'thrustMapping_exp6_';
plotThrustTrial(x(selectPBTobjects_L),[4:6],0)
suffix = '_manualROIs_L';
printAVP
%}



% pol ROI responses to non-pol stimuli:

prefix = 'thrustMapping_exp6_';
plotThrustTrial(x(selectPBTobjects_L),[1:3],0)
suffix = '_manualROIs_L_pol';
printAVP

%{
prefix = 'barMapping_exp5_';
plotBarMapTrial(x(selectPBTobjects_L),[1:3],0)
suffix = '_manualROIs_L_pol';
printAVP
%}

%{
% special case with inhibition/antagonism
prefix = 'thrustMapping_exp6_';
plotThrustTrial(x(selectPBTobjects_L),[7:8],0)
suffix = '_manualROIs_L_dorsal_inhib';
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example image showing corresponding ROIs for timeseries above
loadR73C04_AOTU
 
printpath = fig3path;

superUseMSP(x); 

prefix = 'avgImg_';

plotCombPolImgManual( x(selectPBTobjects_L(1)).MIP ,0,[1:6])  % no cols
suffix = '_manualROIs_L';
printAVP

%{

plotCombPolImgManual( x(selectPBTobjects_L(1)).MIP ,0,[1:8])  % no cols
suffix = '_manualROIs_L1';
printAVP

plotCombPolImgManual( x(selectPBTobjects_L(2)).MIP ,0,[1:8])  % no cols
suffix = '_manualROIs_L2';
printAVP

plotCombPolImgManual( x(selectPBTobjects_L(3)).MIP ,0,[1:8])  % no cols
suffix = '_manualROIs_L3';
printAVP
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example bar/pol polar map (left hand AOTU)

loadR73C04_AOTU
 
printpath = fig3path;

prefix = 'barpolImg_';

superUseMSP(x); 

plotBarPolImg( x(selectPBTobjects_L(1)).MIP ,1,1)  % no layermask, noise filtered
suffix = '_L';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% right side recordings from the same flies as above, showing pol responses
% ROI(1:3) lateral pol ROIs, ventral to dorsal
% ROI(4:6) mid bar ROIs, ventral to dorsal
% ROI(7) medial dorsal strong pol ROI
loadR73C04_AOTU
 
printpath = fig3path;

superUseMSP(x); 

% pol experiment responses
prefix = 'polMapping_exp4_';
plotPolMapTrialManual(x(selectPBTobjects_R),[1:3],0)
suffix = '_manualROIs_R';
printAVP

prefix = 'polPolar_';
plotSinglePolarResp(x(selectPBTobjects_R),[1:3],0,0)
suffix = '_manualROIs_R';
printAVP


% non-pol ROI responses to non-pol stimuli
prefix = 'barMapping_exp5_';
plotBarMapTrial(x(selectPBTobjects_R),[4:6],0)
suffix = '_manualROIs_R';
printAVP

%{
prefix = 'thrustMapping_exp6_';
plotThrustTrial(x(selectPBTobjects_R),[4:6],0)
suffix = '_manualROIs_R';
printAVP
%}


% pol ROI responses to non-pol stimuli:

prefix = 'thrustMapping_exp6_';
plotThrustTrial(x(selectPBTobjects_R),[1:3],0)
suffix = '_manualROIs_R_pol';
printAVP

%{
prefix = 'barMapping_exp5_';
plotBarMapTrial(x(selectPBTobjects_R),[1:3],0)
suffix = '_manualROIs_R_pol';
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example image showing corresponding ROIs for timeseries above
loadR73C04_AOTU
 
printpath = fig3path;

prefix = 'avgImg_';

superUseMSP(x); 

plotCombPolImgManual( x(selectPBTobjects_R(1)).MIP ,0,[1:6])  % no cols
suffix = '_manualROIs_R';
printAVP

%{
% Same images showing ROIs for all recordings

plotCombPolImgManual( x(selectPBTobjects_R(1)).MIP ,0,[1:7])  % no cols
suffix = '_manualROIs_R1';
printAVP

plotCombPolImgManual( x(selectPBTobjects_R(2)).MIP ,0,[1:7])  % no cols
suffix = '_manualROIs_R2';
printAVP

plotCombPolImgManual( x(selectPBTobjects_R(3)).MIP ,0,[1:7])  % no cols
suffix = '_manualROIs_R3';
printAVP
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example bar/pol polar map (right hand AOTU)

loadR73C04_AOTU
 
printpath = fig3path;

prefix = 'barpolImg_';

superUseMSP(x); 

plotBarPolImg( x(selectPBTobjects_R(1)).MIP ,1,1)  % no layermask, noise filtered
suffix = '_R';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Various plots for TuTu driver in AOTU
clear objStr
% objStr{1} = 'R49E09_AOTU'; % TuBu inf driver
% objStr{2} = 'R88A06_AOTU'; % TuBu sup + ant driver
% objStr{3} = 'R34H10_AOTU'; % TuBu ant driver
objStr{4} = 'R17F12_AOTU'; % TuTu driver
% objStr{5} = 'R48B06_AOTU'; % pan-TuBu driver

selectIdx = [ 1, 4, 1, 2, 1];

for tIdx = 4
    
    % load object array
    eval(['load' objStr{tIdx}]) 
    
    pathsAVP
    printpath = fig3path;

    superUseMSP(x,1)
    superPolThreshold(x,-1)
    
    
    % PSI map
    %{
    prefix = 'polSelectivity_';
       
    plotPolSelImg( x(selectObj(selectIdx(tIdx))).MIP ,1,-1)
    suffix = '_R_MIP';
    printAVP
    %}
    
    
    % Tuning map
    prefix = 'polTuning_';
    
    plotCombPolImgManual( x(selectObj(selectIdx(tIdx))).MIP ,1,[1:3],0) 
    suffix = '_R_MIP_mask';
    printAVP
    
     plotCombPolImgManual( x(selectObj(selectIdx(tIdx))).Layers(selectLayer(selectIdx(tIdx))) ,1,[1:3],0) 
    suffix = '_R_Layer_mask';
    printAVP
 
    % Scatter plot
    prefix = 'scatter_vert_';

    superXY(x,0,1,[],[],1,1,0,2)
    suffix = '_tuning';
    printAVP
    
%     superXY(x,0,1,[],[],0,1,0,2)
%     suffix = '_sel';
%     printAVP
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TuTu recording, bar map trial 
loadR17F12_AOTU

printpath = fig3path;

superUseMSP(x); 

prefix = 'barMapping_exp5_';

plotBarMapTrial(x,[1:3],0) % one recording had massive inhibition after pol mapping, so timeseries in bar mapping showed large increase during 1st trial set, so we just plot 2nd set alone
suffix = '_manualROIs_R';

printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Selectivity histograms and boxplots for pol vs no pol median values
clear objStr
objStr{1} = 'R17F12_AOTU';
% objStr{2} = 'R49E09_AOTU';
% objStr{3} = 'R88A06_AOTU';
% objStr{4} = 'R34H10_AOTU';

pathsAVP
printpath = fig3path;
savename = 'TuTu_AOTU';

prefix = 'psi_';

% abs values
b = boxplotPSI(objStr);
suffix = 'box';
printAVP

% probability density
pdfplotPSI(objStr,'add','none');
% legend off
suffix = 'pdf';
printAVP

% mean-of-controls-subtracted values
b = boxplotPSI(objStr,'pol','compare');
suffix = 'diff';
printAVP

