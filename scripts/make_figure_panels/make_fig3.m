% Generate plots for Fig3 (two MeTu drivers in AOTU)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig3path),'dir')
    try rmdir(fullfile(fig3path),'s'),end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 56F07 MeTu scatter plots
% recordings with predominant polarotopic organization
%{
loadR56F07_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

prefix = 'scatter_vert_';

superXY(x(),0,1,0,[],1,0,0,2)
suffix = '_all_tuning';
printAVP
%}
%% 56F07 MeTu scatter plots
% recordings with predominant polarotopic organization

loadR56F07_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

prefix = 'scatter_vert_';

pos_set = setdiff(1:length(x),[selectSetOutliers,selectSetTuTuCompare,selectSetNeg]);
superXY(x(pos_set),0,1,0,[],1,0,1,2)
suffix = '_pos_tuning';
printAVP
% 
% superXY(x(pos_set),0,1,0,[],0,1,0,2)
% suffix = '_pos_sel';
% printAVP

disp(['[R56F07] Positive polarotopy: ' num2str(length([x(pos_set).containsPolMapExp]==4)) '/' num2str(length([x().containsPolMapExp]==4)) ' recordings'])

neg_set = [selectSetTuTuCompare,selectSetNeg];
disp(['[R56F07] Negative polarotopy: ' num2str(length([x(neg_set).containsPolMapExp]==4)) '/' num2str(length([x().containsPolMapExp]==4)) ' recordings'])
%% 56F07 MeTu scatter plots
% recordings with less-frequent negative polarotopic organization
%{
loadR56F07_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

prefix = 'scatter_vert_';

neg_set = [selectSetTuTuCompare,selectSetNeg];
superXY(x(neg_set),0,1,0,[],1,0,0,2)
suffix = '_neg_tuning';
printAVP

superXY(x(neg_set),0,1,0,[],0,1,0,2)
suffix = '_neg_sel';
printAVP
disp(['[R56F07] Negative polarotopy: ' num2str(length([x(neg_set).containsPolMapExp]==4)) '/' num2str(length([x().containsPolMapExp]==4)) ' recordings'])
%}

%% 56F07 MeTu scatter plots
% all recordings, horizontal organization
%{
loadR56F07_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

prefix = 'scatter_horiz_';

superXY(x(),1,0,0,[],0,0,0,2)
suffix = '_all_sel';
printAVP
%}
%% 56F07 tuning map, example1: predominant positive polarotopy
loadR56F07_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

plotCombPolImgManual( x(selectObj(1)).MIP ,1,[],0)  % with  mask applied
prefix = 'polTuning_';
suffix = '_R_MIP_mask_pos';
printAVP

%{
plotCombPolImgManual( x(selectObj(1)).Layers(selectLayer(1)) ,1,[],0)  % with  mask applied
prefix = 'polTuning_';
suffix = '_R_Layer_mask_pos';
printAVP
%}
%% 56F07 tuning map, example2: less-frequent negative polarotopy
%{
loadR56F07_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

plotCombPolImgManual( x(selectObj(3)).MIP ,1,[],0)  % with  mask applied
prefix = 'polTuning_';
suffix = '_R_MIP_mask_neg';
printAVP
%}
%% 56F07 PSI map, for example1 above, including layer mask
loadR56F07_AOTU
 
printpath = fig3path;

plotPolSelImg( x(selectObj(1)).MIP ,1,-1)  % with  mask applied
prefix = 'polSelectivity_';
suffix = '_R_MIP_mask_pos';
printAVP

%% 56F07 PSI map, for example2 above, including layer mask
%{
loadR56F07_AOTU
 
printpath = fig3path;

plotPolSelImg( x(selectObj(3)).MIP ,1,-1)  % with  mask applied
prefix = 'polSelectivity_';
suffix = '_R_MIP_mask_neg';
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 73C04 MeTu scatter plots
% recordings with predominant polarotopic organization
%{
loadR73C04_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

prefix = 'scatter_vert_';

superXY(x(),0,1,0,[],1,0,0,2)
suffix = '_all_tuning';
printAVP
%}
%% 73C04 MeTu scatter plots
% recordings with predominant polarotopic organization

loadR73C04_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

prefix = 'scatter_vert_';

pos_set = setdiff(1:length(x),[selectSetOutliers,selectSetTuTuCompare]);
superXY(x(pos_set),0,1,0,[],1,0,1,2)
suffix = '_pos_tuning';
printAVP

% superXY(x(pos_set),0,1,0,[],0,1,0,2)
% suffix = '_pos_sel';
% printAVP
% 
disp(['[R73C04] Positive polarotopy: ' num2str(length([x(pos_set).containsPolMapExp]==4)) '/' num2str(length([x().containsPolMapExp]==4)) ' recordings'])

other_set = [selectSetOutliers,selectSetTuTuCompare];
disp(['[R73C04] Other polarotopy (outliers): ' num2str(length([x(other_set).containsPolMapExp]==4)) '/' num2str(length([x().containsPolMapExp]==4)) ' recordings'])

%% 73C04 MeTu scatter plots
% recordings with negative polarotopic organization, an organization of tunings that doesn't fit the
% common positive one
%{
loadR73C04_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

prefix = 'scatter_vert_';

other_set = [selectSetOutliers,selectSetTuTuCompare];
superXY(x(other_set),0,1,0,[],1,0,0,2)
suffix = '_other_tuning';
printAVP

superXY(x(other_set),0,1,0,[],0,1,0,2)
suffix = '_other_sel';
printAVP

disp(['[R73C04] Other polarotopy (outliers): ' num2str(length([x(other_set).containsPolMapExp]==4)) '/' num2str(length([x().containsPolMapExp]==4)) ' recordings'])
%}
%% 73C04 MeTu scatter plots
% all recordings, horizontal organization
%{
loadR73C04_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

prefix = 'scatter_horiz_';

superXY(x(),1,0,0,[],0,0,0,2)
suffix = '_all_sel';
printAVP
%}
%% 73C04 tuning map, example1: predominant positive polarotopy
loadR73C04_AOTU
 
printpath = fig3path;

superUseMSP(x,1)
superPolThreshold(x,-1);

plotCombPolImgManual( x(selectObj(1)).MIP ,1,[],0)  % with  mask applied
prefix = 'polTuning_';
suffix = '_R_MIP_mask_pos';
printAVP

%{
plotCombPolImgManual( x(selectObj(1)).Layers(selectLayer(1)) ,1,[],0)  % with  mask applied
prefix = 'polTuning_';
suffix = '_R_Layer_mask_pos';
printAVP
%}

%% 73C04 PSI map, for example1 above, including layer mask
loadR73C04_AOTU
 
printpath = fig3path;

plotPolSelImg( x(selectObj(1)).MIP ,1,-1)  % with  mask applied
prefix = 'polSelectivity_';
suffix = '_R_MIP_mask_pos';
printAVP
