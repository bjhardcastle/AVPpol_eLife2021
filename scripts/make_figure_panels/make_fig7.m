% Generate plots for Fig7 (TuBu/R4m population resultant tunings)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig7path),'dir')
    try rmdir(fullfile(fig7path),'s'),end
end
% % % add left and right resultant plots and print info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TuBu anterior bulb
% When we look at the representation of pol tunings in a TuBu_a population,
% in an individual fly, we see a range of angles with a roughly uniform
% distribution and a short resultant:

loadR34H10_Bu
printpath = fig7path;

superUseMSP(x,1)
superPolThreshold(x,-1)

usePSIweighting = 1;
usePSIthreshold = 0;

disp('[Left+Right recordings]')
% % print plot (pixel-based to compare with R4m EB):
superPolAngResultant(x,usePSIweighting, usePSIthreshold)
suffix = 'pixel';
prefix = 'polResultant_';
printAVP

% % print stats (ROI-based to compare with pixel-based method):
superPolAngResultantROI(x,usePSIweighting, usePSIthreshold)
close(gcf)

disp('[Left recordings]')
% % print stats (pixel-based to compare with R4m EB):
superPolAngResultant(x(inclL),usePSIweighting, usePSIthreshold)
close(gcf)
% % print stats (ROI-based to compare with pixel-based method):
superPolAngResultantROI(x(inclL),usePSIweighting, usePSIthreshold)
close(gcf)

disp('[Right recordings]')
% % print stats (pixel-based to compare with R4m EB):
superPolAngResultant(x(inclR),usePSIweighting, usePSIthreshold)
close(gcf)
% % print stats (ROI-based to compare with pixel-based method):
superPolAngResultantROI(x(inclR),usePSIweighting, usePSIthreshold)
close(gcf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% R4m anterior bulb
% When we look at the representation of pol tunings in an R4m population,
% in an individual fly, we see a range of angles with a less uniform
% distribution and a longer resultant:

loadR34D03_Bu
printpath = fig7path;

superUseMSP(x,1)
superPolThreshold(x,-1)

usePSIweighting = 1;
usePSIthreshold = 0;

disp('[Left+Right recordings]')
% % print plot (pixel-based to compare with R4m EB):
superPolAngResultant(x,usePSIweighting, usePSIthreshold)
suffix = 'pixel';
prefix = 'polResultant_';
printAVP

% % print stats (ROI-based to compare with pixel-based method):
superPolAngResultantROI(x,usePSIweighting, usePSIthreshold)
close(gcf)

disp('[Left recordings]')
% % print stats (pixel-based to compare with R4m EB):
superPolAngResultant(x(inclL),usePSIweighting, usePSIthreshold)
close(gcf)
% % print stats (ROI-based to compare with pixel-based method):
superPolAngResultantROI(x(inclL),usePSIweighting, usePSIthreshold)
close(gcf)

disp('[Right recordings]')
% % print stats (pixel-based to compare with R4m EB):
superPolAngResultant(x(inclR),usePSIweighting, usePSIthreshold)
close(gcf)
% % print stats (ROI-based to compare with pixel-based method):
superPolAngResultantROI(x(inclR),usePSIweighting, usePSIthreshold)
close(gcf) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% R4m EB
% When we look at the representation of pol tunings in an R4m population,
% in an individual fly, we see a range of angles with a less uniform
% distribution and a longer resultant:
loadR34D03_EB
printpath = fig7path;

superUseMSP(x,1)
superPolThreshold(x,-1)

usePSIweighting = 1;
usePSIthreshold = 0;

% Only a single ROI was drawn around the entire EB (which would always give
% a resultant length of 1 in the ROI-based analysis), so we only plot the
% pixel-based approach (and there are no L/R recordings for the EB):
superPolAngResultant(x,usePSIweighting, usePSIthreshold)
suffix = 'pixel';
prefix = 'polResultant_';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% R4m EB selectivity, tuning maps, polar histograms
loadR34D03_EB
printpath = fig7path;
superUseMSP(x,1)
superPolThreshold(x,-1)

usePSIweighting = 1;
usePSIthreshold = 0;
useLayerMask = 1;

for sIdx = [1,3,6,7] %[1,3,6,7,11,13]% [1,3,6,7]
    
    suffix = num2str(sIdx);
    
    prefix = 'polSel';
    plotPolSelImg( x(selectObj(sIdx)).MIP ,1,-1)
    printAVP
    
  
    prefix = 'polHist';
    superPolAngHist(x(selectObj(sIdx)),1, usePSIthreshold);
    printAVP 
    
    if sIdx == 7
         prefix = 'avgActivity';
         plotCombPolImgManual( x(selectObj(sIdx)).MIP ,0,0)
         printAVP
	else
		prefix = 'polTuning';
		plotCombPolImgManual(x(selectObj(sIdx)).MIP ,1,-1,~useLayerMask,1)
		printAVP
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example average pol response from whole EB/all R4m
loadR34D03_EB
printpath = fig7path;
superUseMSP(x,1)

for sIdx = [1,3,6]
    if sIdx == 6
        plotPolMapTrialManual( x(selectObj(sIdx)), 1) % single layer, no MSP backed up due to size: use manually saved ROI of layer_mask
    else
        plotPolMapTrialManual( x(selectObj(sIdx)), -1) % use layer mask
    end
    ylim([-0.5 0.5])
    scalebarF(gca)
    scalebarT(gca)
    prefix = 'polMapping_';
    suffix = ['layerMask_' num2str(sIdx)];
    printAVP
end