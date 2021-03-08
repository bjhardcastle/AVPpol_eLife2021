% Generate plots for Fig5 (various PSI distributions in TuBu drivers in BU)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig5path),'dir')
    try rmdir(fullfile(fig5path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Selectivity histograms and boxplots for pol vs no pol median values
% First we looked at the average pol selectivity index in each compartment,
% exploiting the fact that 88A06 labels TuBu neurons projecting to both the 
% superior and anterior bulb

clear objStr
objStr{1} = 'R49E09_Bu';
objStr{2} = 'R88A06_Bu';
objStr{3} = 'R88A06_Bu_ant';
objStr{4} = 'R34H10_Bu';

pathsAVP
printpath = fig5path;
savename = 'TuBu_Bu';
prefix = 'psi_';

% abs values - whole layer mask
b = boxplotPSI(objStr,'dist','mask');
suffix = 'box_layer';
printAVP

% mean-of-controls-subtracted values  - whole layer mask
b = boxplotPSI(objStr(2:4),'pol','compare');
suffix = 'diff';
printAVP


% probability density - whole layer mask
pdfplotPSI(objStr,'add','none','dist','mask');
suffix = 'pdf_layer';
printAVP

% abs values 
b = boxplotPSI(objStr);
suffix = 'box';
printAVP

% mean-of-controls-subtracted values
b = boxplotPSI(objStr(2:4),'pol','compare');
suffix = 'diff';
printAVP

% probability density
pdfplotPSI(objStr,'add','none');
suffix = 'pdf';
printAVP

% probability density
pdfplotPSI(objStr,'add','none');
suffix = 'pdf';
printAVP

%
for tIdx = 2:length(objStr)

% mean-of-controls-subtracted values
b = boxplotPSI(objStr{tIdx},'pol','compare');
savename = objStr{tIdx};
suffix = 'diff';
printAVP

% probability density
pdfplotPSI(objStr{tIdx},'add','ctrlCell');
% legend off
savename = objStr{tIdx};
suffix = 'pdf';
printAVP

end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spatial maps of PSI and distribution of tunings - superior bulb

% We then looked at spatial distribution of pol responses in sup bulb: 
% Although mostly unresonsive, some glomeruli showed strong responses,
% shown by high PSI:

loadR88A06_Bu
printpath = fig5path;
superUseMSP(x,1)
superPolThreshold(x,-1); 

useTuningCols = 1;
useLayerMask = 1;

suffix = 'R';
prefix = 'polSel';
plotPolSelImg( x(selectObj(3)).MIP ,1,-1) 
printAVP

prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(3)).MIP ,useTuningCols,[],~useLayerMask,1) 
printAVP

suffix = 'L';
prefix = 'polSel';
plotPolSelImg( x(selectObj(4)).MIP ,1,-1) 
printAVP

prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(4)).MIP ,useTuningCols,[],~useLayerMask,1) 
printAVP

% When we look at the range of pol angles represented in the sup bulb, we
% find an over-representation of near-vertical angles, not uniform
% distribution. This distribution is mirrored in the left and right sides
% of the brain:

useTuningCols = 1;
usePSIthreshold = 1;

prefix = 'polTuning';

suffix = 'hist_R';
superPolAngHist(x(inclR),useTuningCols,usePSIthreshold)
printAVP

suffix = 'hist_L';
superPolAngHist(x(inclL),useTuningCols,usePSIthreshold)
printAVP

%% - and tuning maps for the anterior bulb, using the alternative layer
% masks. These are from the same recordings, by will have to be stitched
% together later. Saves making a whole new object with masks drawn around
% sup+ant bulb just for this figure.
loadR88A06_Bu_ant
printpath = fig5path;
superUseMSP(x,1)
superPolThreshold(x,-1);

useTuningCols = 1;
useLayerMask = 1;

suffix = 'R';
prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(3)).MIP ,useTuningCols,[],~useLayerMask,1) 
printAVP

suffix = 'L';
prefix = 'polTuning';
plotCombPolImgManual( x(selectObj(4)).MIP ,useTuningCols,[],~useLayerMask,1) 
printAVP

%% distribution of tunings - anterior bulb
% When we look at the range of pol angles represented in the ant bulb, we
% find a roughly uniform distribution on each side of the brain:

loadR34H10_Bu
printpath = fig5path;
superUseMSP(x,1)
superPolThreshold(x,-1);

useTuningCols = 1;
usePSIthreshold = 1;

prefix = 'polTuning';

suffix = 'hist_R';
superPolAngHist(x(inclR),useTuningCols,usePSIthreshold)
printAVP

suffix = 'hist_L';
superPolAngHist(x(inclL),useTuningCols,usePSIthreshold)
printAVP
