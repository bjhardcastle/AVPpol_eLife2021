function PSI = getPSIstruct(objarray,excludeSingleLayerRecs)
if nargin < 2 || isempty(excludeSingleLayerRecs)
    excludeSingleLayerRecs = 0;
end
normMethod = 'probability';

histBinEdges = objarray(1).MIP.polSelHistBins;
histBinEdgesRelative = sort(unique([-histBinEdges histBinEdges]));

uniqueFlies = unique([objarray.Fly]);
PSI = struct; % non-existent data will be nan values. inclFly numbers will be empty
for fIdx = 1:length(uniqueFlies)
    
    % For each fly, pool data from all of its recordings, then find counts
    % and median values
    bkgPSIvectorExp = [];
    cellPSIvectorExp =[];
    bkgPSIvectorCtrl = [];
    cellPSIvectorCtrl =[];    
    layerMaskPSIvectorExp = [];
    layerMaskPSIvectorCtrl = [];
    for n = find([objarray.Fly]==uniqueFlies(fIdx))
        
        % check mask exists
        loadLayerMasks(objarray(n).MIP)
        if isempty( objarray(n).MIP.layerMask )
            disp(['No layer mask: x(' num2str(n) ')'])
            maskLayers(objarray(n))
            [~] = input('waiting..');
            loadLayerMasks(objarray(n).MIP)
        end
        
        if objarray(n).containsPolMapExp >= 4
            if ~excludeSingleLayerRecs || objarray(n).MIP.numZPlanes > 1
                % Get PSI values within cells and background of MSP
                bkgPSIvectorExp = [bkgPSIvectorExp; getPolSelMaskedVector( objarray(n).MIP, 'bkg' )];
                cellPSIvectorExp = [cellPSIvectorExp; getPolSelMaskedVector( objarray(n).MIP, 'cell' )];
                layerMaskPSIvectorExp = [layerMaskPSIvectorExp; getPolSelMaskedVector( objarray(n).MIP, 'all' )];
            end
        elseif objarray(n).containsPolMapExp == 2
            if ~excludeSingleLayerRecs || objarray(n).MIP.numZPlanes > 1
                % Get PSI values within cells and background of MSP
                bkgPSIvectorCtrl = [bkgPSIvectorCtrl; getPolSelMaskedVector( objarray(n).MIP, 'bkg' )];
                cellPSIvectorCtrl = [cellPSIvectorCtrl; getPolSelMaskedVector( objarray(n).MIP, 'cell' )];
                layerMaskPSIvectorCtrl = [layerMaskPSIvectorCtrl; getPolSelMaskedVector( objarray(n).MIP, 'all' )];
            end
        end
    end
    
    % Keep track of which flies contributed data (if no data, fly numbers will
    % be empty)
    if ~isempty(bkgPSIvectorExp)
        PSI(fIdx).expInclFlies = [uniqueFlies(fIdx)];
    end
    if ~isempty(bkgPSIvectorCtrl)
        PSI(fIdx).ctrlInclFlies = [uniqueFlies(fIdx)];
    end
    
    % Now collect data for this fly
    
    
    % Store the raw data in case we need it later
    PSI(fIdx).expBkgData = bkgPSIvectorExp;
    PSI(fIdx).expCellData = cellPSIvectorExp;
    PSI(fIdx).expLayerMaskData = layerMaskPSIvectorExp;
    PSI(fIdx).ctrlBkgData = bkgPSIvectorCtrl;
    PSI(fIdx).ctrlCellData = cellPSIvectorCtrl;
    PSI(fIdx).ctrlLayerMaskData = layerMaskPSIvectorCtrl;

    % Store the median of each distribution and their difference:
    PSI(fIdx).expBkgMedian = nanmedian(bkgPSIvectorExp);
    PSI(fIdx).expCellMedian = nanmedian(cellPSIvectorExp);
    PSI(fIdx).expLayerMaskMedian = nanmedian(layerMaskPSIvectorExp);
    PSI(fIdx).ctrlBkgMedian = nanmedian(bkgPSIvectorCtrl);
    PSI(fIdx).ctrlCellMedian = nanmedian(cellPSIvectorCtrl);
    PSI(fIdx).ctrlLayerMaskMedian = nanmedian(layerMaskPSIvectorCtrl);

    
    % Collect counts for histograms:
    % absolute values (bins 0 to 1)
    PSI(fIdx).expBkgCts = histcounts( bkgPSIvectorExp , histBinEdges,'Normalization',normMethod)';
    PSI(fIdx).expCellCts = histcounts( cellPSIvectorExp , histBinEdges,'Normalization',normMethod)';
    PSI(fIdx).expLayerMaskCts = histcounts( layerMaskPSIvectorExp , histBinEdges,'Normalization',normMethod)';
    PSI(fIdx).ctrlBkgCts = histcounts( bkgPSIvectorCtrl , histBinEdges,'Normalization',normMethod)';
    PSI(fIdx).ctrlCellCts = histcounts( cellPSIvectorCtrl , histBinEdges,'Normalization',normMethod)';
    PSI(fIdx).ctrlLayerMaskCts = histcounts( layerMaskPSIvectorCtrl , histBinEdges,'Normalization',normMethod)';

    % relative values (median of background subtracted): (using bins -1 to 1)
    PSI(fIdx).expBkgCtsRelative = histcounts( bkgPSIvectorExp - PSI(fIdx).expBkgMedian , histBinEdgesRelative,'Normalization',normMethod)';
    PSI(fIdx).expCellCtsRelative = histcounts( cellPSIvectorExp - PSI(fIdx).expBkgMedian , histBinEdgesRelative,'Normalization',normMethod)';
    PSI(fIdx).expLayerMaskCtsRelative = histcounts( layerMaskPSIvectorExp - PSI(fIdx).expLayerMaskMedian , histBinEdgesRelative,'Normalization',normMethod)';
    PSI(fIdx).ctrlBkgCtsRelative = histcounts( bkgPSIvectorCtrl - PSI(fIdx).ctrlBkgMedian , histBinEdgesRelative,'Normalization',normMethod)';
    PSI(fIdx).ctrlCellCtsRelative = histcounts( cellPSIvectorCtrl - PSI(fIdx).ctrlBkgMedian , histBinEdgesRelative,'Normalization',normMethod)';
    PSI(fIdx).ctrlLayerMaskCtsRelative = histcounts( layerMaskPSIvectorCtrl - PSI(fIdx).ctrlLayerMaskMedian , histBinEdgesRelative,'Normalization',normMethod)';
  
    PSI(fIdx).histBinEdgesRelative = histBinEdgesRelative;
    PSI(fIdx).histBinEdges = histBinEdges;

end

end
