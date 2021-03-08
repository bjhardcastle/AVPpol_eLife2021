% Generate plots for Fig5 (TuTu characterization)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig5path),'dir')
    try rmdir(fullfile(fig5path),'s'),end
end

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
    printpath = fig5path;

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

printpath = fig5path;

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
printpath = fig5path;
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

