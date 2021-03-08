% Generate plots for fig4 (pol+nonpol responses single MeTu driver in AOTU)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig4path),'dir')
    try rmdir(fullfile(fig4path),'s'),end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Three example left side recordings, showing retinotopy in non-pol responses
% ROI(1:3) lateral pol ROIs, ventral to dorsal
% ROI(4:6) mid bar ROIs, ventral to dorsal
% ROI(7) dorsal ROI inhibited by everything
% ROI(8) dorsal ROI, thrust resp

loadR73C04_AOTU
 
printpath = fig4path;


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
 
printpath = fig4path;

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
 
printpath = fig4path;

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
 
printpath = fig4path;

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
 
printpath = fig4path;

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
 
printpath = fig4path;

prefix = 'barpolImg_';

superUseMSP(x); 

plotBarPolImg( x(selectPBTobjects_R(1)).MIP ,1,1)  % no layermask, noise filtered
suffix = '_R';
printAVP

