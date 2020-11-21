% Generate plots for Fig1s1 (R7/R8 + DmDRA in MEDRA)
% Each cell can be run independently
pathsAVP
if exist(fullfile(fig1s1path),'dir')
    try rmdir(fullfile(fig1s1path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RHS pol average intensity mas no mask
loadR7R8
 
printpath = fig1s1path;

prefix = 'polTuning_wholeME_bkgthresh';

superUseMSP(x,1)
superPolThreshold(x,-1);

plotCombPolImg( x(selectObj(4)).MIP ,1,[],1)  % no ROIs
suffix = '_noMask_R';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RHS single example R7/R8 polarization tuning map
loadR7R8
printpath = fig1s1path;

prefix = 'polTuning_';

superUseMSP(x,1)
superPolThreshold(x,-1); 

plotCombPolImgManual( x(selectObj(4)).MIP ,1,[],0,1)  % with mask applied
suffix = '_noROIs_R';
printAVP


% plotCombPolImgManual( x(selectObj(4)).MIP ,1,[],1,1)  % no mask
% suffix = '_noROIs_noMask_R';
% printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ant-post position vs pol tuning scatter plot. example R7/R8 recording. MSP
loadR7R8
printpath = fig1s1path;

prefix = 'scatter_horiz_';

superUseMSP(x,1)
superPolThreshold(x,-1);

superXY(x(selectObj(4)),1,[],[],[],0,[],[],2)
suffix = '_sel';
printAVP


superXY(x(selectObj(4)),1,[],[],[],1,[],[],2)
suffix = '_tuning';
printAVP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RHS pol selectivity map
%{
loadR7R8
 
printpath = fig1s1path;

plotPolSelImg( x(selectObj(4)).MIP ,1,[])
prefix = 'polSel_wholeME';
suffix = '_noROIs_R';
printAVP

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RHS pol tuning map, mid-medra, no pol vs pol, no mask, bkgthresh
loadR7R8
 
printpath = fig1s1path;

prefix = 'polTuning_MidMEDRA_bkgthresh';

superUseMSP(x,1)
superPolThreshold(x,-1);

% activity image
plotCombPolImgManual( x(selectObj(2)).MIP ,1,[],1,1)  % no pol
suffix = '_no_pol_R';
printAVP

plotCombPolImgManual( x(selectObj(3)).MIP ,1,[],1,1)  % pol
suffix = '_pol_R';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RHS pol tuning map, mid-medra Fly2, no pol vs pol, no threshold
loadR7R8
 
printpath = fig1s1path;

prefix = 'polTuning_MidMEDRA1_nothresh_';

superUseMSP(x,1)
superPolThreshold(x,0);

% activity image
plotCombPolImgManual( x(selectObj(2)).MIP ,1,[],1,1)  % no pol
suffix = '_no_pol_R';
printAVP

plotCombPolImgManual( x(selectObj(3)).MIP ,1,[],1,1)  % pol
suffix = '_pol_R';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RHS pol selectivity map, mid-medra Fly1, no pol vs pol
loadR7R8
 
printpath = fig1s1path;

prefix = 'polSel_MidMEDRA1';

superUseMSP(x,1)
superPolThreshold(x,0);

% selectivity image
plotPolSelImg( x(selectObj(2)).MIP ,1,[])  % no pol
suffix = '_no_pol_R';
printAVP

plotPolSelImg( x(selectObj(3)).MIP ,1,[])  % pol
suffix = '_pol_R';
printAVP



%% RHS pol tuning map, mid-medra Fly2, no pol vs pol, bkg thresh
loadR7R8
 
printpath = fig1s1path;

prefix = 'polTuning_MidMEDRA2_bkgthresh';

superUseMSP(x,1)
superPolThreshold(x,-1);

% activity image
plotCombPolImgManual( x(selectObj(5)).MIP ,1,[],1,1)  % no pol
suffix = '_no_pol_R';
printAVP

plotCombPolImgManual( x(selectObj(6)).MIP ,1,[],1,1)  % pol
suffix = '_pol_R';
printAVP

%% RHS pol tuning map, mid-medra Fly2, no pol vs pol, no threshold
loadR7R8
 
printpath = fig1s1path;

prefix = 'polTuning_MidMEDRA2_nothresh';

superUseMSP(x,1)
superPolThreshold(x,0);

% activity image
plotCombPolImg( x(selectObj(5)).MIP ,1,[],1)  % no pol
suffix = '_no_pol_R';
printAVP

plotCombPolImg( x(selectObj(6)).MIP ,1,[],1)  % pol
suffix = '_pol_R';
printAVP

%% RHS time-series, anterior and posterior columns, R7/R8 pairs, Fly2
loadR7R8
 
printpath = fig1s1path;

superUseMSP(x,1)

suffix = '_manualROIs_R_post';
plotPolMapTrialManual( x(selectObj(6)) , [1,2] , 0,[0 3])
prefix = 'polMapping_exp4_MidMEDRA2';
printAVP

plotSinglePolarResp( x(selectObj(6)) , [1,2],[],0)
prefix = 'polPolar_exp4_MidMEDRA2';
printAVP

suffix = '_manualROIs_R_ant';
plotPolMapTrialManual( x(selectObj(6)) , [3,4] , 0,[0 3])
prefix = 'polMapping_exp4_MidMEDRA2';
printAVP

plotSinglePolarResp( x(selectObj(6)) , [3,4],[],0)
prefix = 'polPolar_exp4_MidMEDRA2';
printAVP

%% RHS time-series, anterior and posterior columns, R7/R8 pairs, whole medulla example
loadR7R8
 
printpath = fig1s1path;

superUseMSP(x,1)

prefix = 'polMapping_exp4_wholeME';

times = [20:50];

figure('color','w')
hold on

fps = (x(selectObj(4)).MIP.expStop(4)-x(selectObj(4)).MIP.expStart(4))/132;
F = x(selectObj(4)).MIP.ROI(1).response;
F0 =  sqrt(nanmean(F.^2));
T = linspace(0,length(F)/fps,length(F));
p = plot(T(times),F(times)/F0-1);

F = x(selectObj(4)).MIP.ROI(2).response;
F0 =  sqrt(nanmean(F.^2));
p(2) = plot(T(times),F(times)/F0-1);

daspect(gca,[30,1,1])
ylim([-0.8 0.5])
set(gca,'Clipping','off')
setAVPaxes(gca,3,3)
scalebarT(gca)
scalebarF(gca)
tightfig(gcf)

figure('color','w')
hold on

F = x(selectObj(4)).MIP.ROI(3).response;
F0 =  sqrt(nanmean(F.^2));
p(3) = plot(T(times),F(times)/F0-1);

F = x(selectObj(4)).MIP.ROI(4).response;
F0 =  sqrt(nanmean(F.^2));
p(4) = plot(T(times),F(times)/F0-1);

ylim([-0.8 0.5])
set(gca,'Clipping','off')
daspect(gca,[30,1,1])
setAVPaxes(gca,3,3)
scalebarT(gca)
scalebarF(gca)
tightfig(gcf)

for pidx = 1:length(p)
    p(pidx).LineWidth = 1;
    p(pidx).Color = x(selectObj(4)).MIP.ROI(pidx).color;
end
suffix = '_manualROIs_R_ant';
printAVP

suffix = '_manualROIs_R_post';
printAVP

prefix = 'polPolar_exp4_wholeME';
plotSinglePolarResp( x(selectObj(4)) , [3,4],[],0)
suffix = '_manualROIs_R_ant';
printAVP
plotSinglePolarResp( x(selectObj(4)) , [1,2],[],0)
suffix = '_manualROIs_R_post';
printAVP

prefix = 'polTuning_wholeME_';
plotCombPolImgManual(x(selectObj(4)).MIP,1,1:4,1,1)
suffix = 'manualROIs_R';
printAVP
prefix = 'avgActivity_wholeME_';
plotCombPolImg(x(selectObj(4)).MIP,0)
suffix = '_R';
printAVP
prefix = 'polTuning_wholeME_';
plotCombPolImg(x(selectObj(4)).MIP,1)
suffix = '_R';
printAVP

%% RHS pol selectivity map, mid-medra, no pol vs pol
loadR7R8
 
printpath = fig1s1path;

prefix = 'polSel_MidMEDRA2_';

% selectivity image
plotPolSelImg( x(selectObj(5)).MIP ,1,[])  % no pol
suffix = '_no_pol_R';
printAVP

plotPolSelImg( x(selectObj(6)).MIP ,1,[])  % pol
suffix = '_pol_R';
printAVP


%% RHS pol tuning map, ant,mid,post-medra Fly1 R
loadR7R8
 
printpath = fig1s1path;

superPolThreshold(x,0.3);

% activity image
prefix = 'zoomed_';

plotCombPolImgManual( x(selectObj(7)).MIP ,1,[],1,1)  % ant
suffix = '_1_R_ant_tuning';
printAVP


plotCombPolImgManual( x(selectObj(8)).MIP ,1,[],1,1)  % mid
suffix = '_1_R_mid_tuning';
printAVP


plotCombPolImgManual( x(selectObj(9)).MIP ,1,[],1,1)  % post
suffix = '_1_R_post_tuning';
printAVP


% selectivity image.
prefix = 'zoomed_';

plotPolSelImg( x(selectObj(7)).MIP ,1,[])  % pol
suffix = '_1_R_ant_sel';
printAVP


plotPolSelImg( x(selectObj(8)).MIP ,1,[])  % pol
suffix = '_1_R_mid_sel';
printAVP


plotPolSelImg( x(selectObj(9)).MIP ,1,[])  % pol
suffix = '_1_R_post_sel';
printAVP


%% LHS pol tuning map, ant,mid,post-medra, Fly1 L
loadR7R8
 
printpath = fig1s1path;

superPolThreshold(x,-1);

% activity image
prefix = 'zoomed_';

plotCombPolImgManual( x(selectObj(10)).MIP ,1,[],1,1)  % ant
suffix = '_1_L_ant_tuning';
printAVP


plotCombPolImgManual( x(selectObj(11)).MIP ,1,[],1,1)  % mid
suffix = '_1_L_mid_tuning';
printAVP


plotCombPolImgManual( x(selectObj(12)).MIP ,1,[],1,1)  % post
suffix = '_1_L_post_tuning';
printAVP


% selectivity image.
prefix = 'zoomed_';

plotPolSelImg( x(selectObj(10)).MIP ,1,[])  % pol
suffix = '_1_L_ant_sel';
printAVP


plotPolSelImg( x(selectObj(11)).MIP ,1,[])  % pol
suffix = '_1_L_mid_sel';
printAVP


plotPolSelImg( x(selectObj(12)).MIP ,1,[])  % pol
suffix = '_1_L_post_sel';
printAVP


%% RHS pol tuning map, ant,mid,post-medra Fly2 R
loadR7R8
 
printpath = fig1s1path;

superPolThreshold(x,-1);

% activity image
prefix = 'zoomed_';

plotCombPolImgManual( x(selectObj(13)).MIP ,1,[],1,1)  % ant
suffix = '_2_R_tuning_ant';
printAVP


plotCombPolImgManual( x(selectObj(14)).MIP ,1,[],1,1)  % mid
suffix = '_2_R_mid_tuning';
printAVP


plotCombPolImgManual( x(selectObj(15)).MIP ,1,[],1,1)  % post
suffix = '_2_R_post_tuning';
printAVP


% selectivity image.
prefix = 'zoomed_';

plotPolSelImg( x(selectObj(13)).MIP ,1,[])  % pol
suffix = '_2_R_ant_sel';
printAVP


plotPolSelImg( x(selectObj(14)).MIP ,1,[])  % pol
suffix = '_2_R_mid_sel';
printAVP


plotPolSelImg( x(selectObj(15)).MIP ,1,[])  % pol
suffix = '_2_R_post_sel';
printAVP

%% LHS DmDRA polarization tuning maps 
loadDmDRA
printpath = fig1s1path;

prefix = 'polTuning_';

superUseMSP(x,1)
loadLayerMasks([x.MIP])
superPolThreshold(x,-1); 

plotCombPolImg( x(selectObj(2)).MIP ,1)  % with mask applied
suffix = 'L_MIP_mask';
printAVP

%{
plotCombPolImgManual( x(selectObj(2)).MIP ,1,[],1,1)  % no mask
suffix = 'L_MIP_nomask';
printAVP
%}

