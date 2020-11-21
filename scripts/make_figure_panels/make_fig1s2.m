% Generate plots for Fig1s2 (DmDRA opponency)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig1s2path),'dir')
    try rmdir(fullfile(fig1s2path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example image showing an auto ROI with 90deg pref tuning
loadDmDRA
printpath = fig1s2path;

prefix = 'avgImg_exp4_';
plotCombPolImgManual( x(selectObj(6)).MIP ,0,90)  % no cols
suffix = '_ROI90_R';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example image and corresponding layermask used for flash response as control
loadDmDRA
printpath = fig1s2path;

% average activity image
prefix = 'avgImg_exp2_';
plotCombPolImgManual( x(selectObj(7)).MIP ,0,0)  % no cols
suffix = '';
printAVP

% layermask
mask = x(selectObj(7)).MIP.ROI(1).mask;
figure('color','w')
imagesc(imrotate(mask,90))
axis 'image'
axis 'off'
colormap(gca,[1 1 1;1 0 0]) % white background will be transparent in pdf
getAVPplotParams
setAVPaxes(gca,defaultImageHeight_cm)
tightfig(gcf)
suffix = 'cellmask_R';
printAVP

%% polar tuning map (exp4) for auto-ROIs with 90deg tuning
loadDmDRA
printpath = fig1s2path;

plotSinglePolarResp(x,[90],[],0)
prefix = 'polPolar_exp4';
suffix = '_ROI90_R';
printAVP

%% polar tuning map (exp2) for entire cellmask as controls
loadDmDRA
printpath = fig1s2path;

plotSinglePolarResp(x, 1,[],0,1)
prefix = 'polPolar_exp2';
suffix = '_cellmask_R';
printAVP

%% flash responses (exp3) for auto-ROIs with 90 deg tuning
loadDmDRA
printpath = fig1s2path;

outputData = plotFlashResp( x , [ 90 ]);
prefix = 'polFlash_exp3_';
pol_ylim = ylim; % apply ylim to control plot below
patch_ylim = outputData(1,1).patch(1).Vertices([1,2],2);
% two plots produced
suffix = '_ROI90_pol90';
printAVP
suffix = '_ROI90_pol0';
printAVP

unpolData = plotFlashRespManual( x , 1, 0, 1); % controls,cell mask, unknown tuning
% apply lims from above, and adjust patch 
ylim(pol_ylim);
unpolData(1).patch.Vertices([1,4],2) = [1,1].*patch_ylim(1);
unpolData(1).patch.Vertices([2,3],2) = [1,1].*patch_ylim(2);
prefix = 'polFlash_exp1_';
suffix = '_cellmask';
printAVP


for ROIidx = 1:size(outputData,1)
    fprintf('ROI(%d): 0deg flash: max=%1.3f, min=%1.3f (of mean trial)\n',ROIidx,outputData(ROIidx,1).maxOfMean,outputData(ROIidx,1).minOfMean);
    fprintf('ROI(%d): 90deg flash: max=%1.3f, min=%1.3f (of mean trial)\n',ROIidx,outputData(ROIidx,2).maxOfMean,outputData(ROIidx,2).minOfMean);
    [p] = ranksum(outputData(ROIidx,1).maxOfFlies, -outputData(ROIidx,2).minOfFlies);
    fprintf('max(0) vs -min(90): p=%1.6f (ranksum), N=%d\n',p,outputData(1,1).populationSize)
end

fprintf('unpol flash [cellmask]: max=%1.3f, min=%1.3f (of mean trial, N=%d)\n',unpolData(ROIidx,1).maxOfMean,unpolData(ROIidx,1).minOfMean,unpolData(1,1).populationSize);



