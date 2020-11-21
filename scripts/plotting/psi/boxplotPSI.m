function [notBoxStruct,ax] = boxplotPSI(varargin)
% (see end for input arguments)

% Initialize input arguments 
ax = [];
PSI = struct();
objNames = {};
plotDistribution = [];
polVals = '';
useCommonCellColor = [];
parseInputs


% Get common plotting values
params = getAVPplotParamsWrapper();
objCols = params.objCols;
lightGreyCol = params.lightGreyCol;
darkGreyCol = params.darkGreyCol;
defaultAxisHeight_cm = params.defaultAxisHeight_cm;
defaultMarkerSize = params.defaultMarkerSize;


objSel = ismember(objNames(:),fieldnames(PSI));
assert( sum(objSel)>0 ,'No valid object names entered. Check for missing "_BU" or typos')
nonNanCells = @(c) ~all(isnan(c(:))); % also captures empty cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot values for background regions

for bIdx = 1:size(objNames,1)
    if objSel(bIdx)

        % Setup groups (single number corresponding to x-position in plot,
        % one per data point - in this case one fly) and fetch data
        switch polVals
            case 'exp' % default. Exp4 values
                bkgXpos{bIdx} = repmat(bIdx,1,length([PSI.(objNames{bIdx}).expBkgMedian]));
                bkgData{bIdx} = [PSI.(objNames{bIdx}).expBkgMedian];
            case 'ctrl' % Exp2 values
                bkgXpos{bIdx} = repmat(bIdx,1,length([PSI.(objNames{bIdx}).ctrlBkgMedian]));
                bkgData{bIdx} = [PSI.(objNames{bIdx}).ctrlBkgMedian];
            case 'compare' % Exp4 - mean(Exp2)
                bkgXpos{bIdx} = repmat(bIdx,1,length([PSI.(objNames{bIdx}).expBkgMedian]));
                bkgData{bIdx} = [PSI.(objNames{bIdx}).expBkgMedian] - nanmean([PSI.(objNames{bIdx}).ctrlBkgMedian]);
        end
    end
end

% Get CI size with Bonferroni correction applied
grpLength = sum(cellfun(nonNanCells,bkgData));
alphaVal = 0.05/grpLength;
[bkgHandles,bkgStats] = notBoxPlot([bkgData{:}],[bkgXpos{:}],'style','sdline','markMedian',true,'jitter',0.7,'manualCI',alphaVal);

hold(ax,'on')

for bIdx = 1:size(objNames(objSel),1)
    
    if ~isempty(bkgHandles(bIdx).data)
        bkgHandles(bIdx).semPtch.FaceColor = darkGreyCol;
        bkgHandles(bIdx).semPtch.EdgeColor = 'none';
        bkgHandles(bIdx).semPtch.FaceAlpha = 0.8;
        bkgHandles(bIdx).semPtch.EdgeAlpha = 0.8;

        bkgHandles(bIdx).data.Marker = '.';
        bkgHandles(bIdx).data.Color = 'k';
        
        bkgHandles(bIdx).data.MarkerSize = defaultMarkerSize;
        
        bkgHandles(bIdx).mu.Color = 'k';
        bkgHandles(bIdx).med.Color ='k';
        
        bkgHandles(bIdx).sd.Visible = 'off';
        
        bkgHandles(bIdx).mu.LineWidth = 0.5;
        bkgHandles(bIdx).med.LineWidth = 0.5;
        bkgHandles(bIdx).med.LineStyle = ':';
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot values for cell regions

for cIdx = 1:size(objNames,1)
    if objSel(cIdx)
        % Setup groups (single number corresponding to x-position in plot,
        % one per data point - in this case one fly)
        cellXpos{cIdx} = repmat(cIdx,1,length([PSI.(objNames{cIdx}).expCellMedian]));
        
        % get distribution (subset or whole mask)
        if strcmp(plotDistribution,'cells')
            distStr = 'Cell';
        elseif strcmp(plotDistribution,'mask')
            distStr = 'LayerMask';
        end
        
        % fetch data
        switch polVals
            case 'exp' % default. Exp4 values
                cellData{cIdx} = [PSI.(objNames{cIdx}).(['exp' distStr 'Median'])];
            case 'ctrl' % Exp2 values
                cellData{cIdx} = [PSI.(objNames{cIdx}).(['ctrl' distStr 'Median'])];
            case 'compare' % Exp4 - mean(Exp2)
                cellData{cIdx} = [PSI.(objNames{cIdx}).(['exp' distStr 'Median'])] - nanmean([PSI.(objNames{cIdx}).(['ctrl' distStr 'Median'])]);
        end
    end
end

hold(ax,'on')

% Get CI size with Bonferroni correction applied
grpLength = sum(cellfun(nonNanCells,cellData));
alphaVal = 0.05/grpLength;
[cellHandles,cellStats] = notBoxPlot([cellData{:}],[cellXpos{:}],'style','sdline','markMedian',true,'jitter',0.7,'manualCI',alphaVal);
for cIdx = 1:size(objNames(objSel),1)
    if ~isempty(cellHandles(cIdx).data)
        
        usedNames = objNames(objSel);
        if ~useCommonCellColor && any( contains(fieldnames(objCols),usedNames{cIdx}) )
            cellCol = objCols.(usedNames{cIdx}) ;
        else
            cellCol = 'r';
        end
        
        cellHandles(cIdx).semPtch.FaceColor = lightGreyCol;
        cellHandles(cIdx).semPtch.EdgeColor = 'none';
        cellHandles(cIdx).semPtch.FaceAlpha = 0.8;
        cellHandles(cIdx).semPtch.EdgeAlpha = 0.8;
        
        cellHandles(cIdx).data.Marker = '.';
        cellHandles(cIdx).data.Color = darkGreyCol;
        
        cellHandles(cIdx).data.MarkerSize = defaultMarkerSize;
        
        cellHandles(cIdx).mu.Color = cellCol;
        cellHandles(cIdx).med.Color =cellCol;
        
        cellHandles(cIdx).sd.Visible = 'off';
        
        cellHandles(cIdx).mu.LineWidth = 0.5;
        cellHandles(cIdx).med.LineWidth = 0.5;
        cellHandles(cIdx).med.LineStyle = ':';
    end
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Adjust plot 

uistack([bkgHandles(:).data],'top')
uistack([cellHandles(:).data],'top')
uistack([bkgHandles(:).mu],'top')
uistack([bkgHandles(:).med],'top')
uistack([cellHandles(:).mu],'top')
uistack([cellHandles(:).med],'top')

objNums = [1:size(objNames,1)];
objNameEmpty = ~strcmp('',objNames(:));
ax.XTick = objNums(objNameEmpty);
ax.XTickLabel = objNames(objNameEmpty);
ax.XTickLabelRotation = 45;
ax.XAxis.TickLabelInterpreter = 'none';

switch polVals
    
    case{'exp','ctrl'}
        ax.YLim = [0 1];
        ax.YTick = [0:0.2:1];
        ax.YLabel.String = 'PSI';

    case 'compare'
        ax.YLim = [-0.2 0.8];
        ax.YTick = [-0.2:0.2:0.8];
        ax.YLabel.String = ['\Delta' 'PSI'];
end

ax.YLabel.Rotation = 90;

ax.XLim = [0.3 size(objNames,1)+0.7];
ax.YGrid = 'on';
daspect([6,1,1])
setAVPaxes(ax,defaultAxisHeight_cm)
tightfig(gcf)
offsetAxesXTickONLY(ax)
ax.Layer = 'bottom';

addExportFigToolbar(gcf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup structure for output

notBoxStruct = struct;
assert(size(objNames(objSel),1) == size(cellStats,2),'Mismatch between plot and data: cannot assign outputs')
nonEmptyBkgEntries = cellfun(nonNanCells,bkgData);
nonEmptyCellEntries = cellfun(nonNanCells,cellData);
fN=0;
for fIdx= 1:size(objNames,1)
    if objSel(fIdx)
        fN=fN+1;
        if nonEmptyCellEntries(fIdx) % Only return non-empty entries
            notBoxStruct.(objNames{fIdx}).cellStats = cellStats(fN);
            notBoxStruct.(objNames{fIdx}).cellHandles = cellHandles(fN);
        end
        if nonEmptyCellEntries(fIdx)
            notBoxStruct.(objNames{fIdx}).bkgStats = bkgStats(fN);
            notBoxStruct.(objNames{fIdx}).bkgHandles = bkgHandles(fN);
        end
        
        if nonEmptyCellEntries(fIdx)
        if ~strcmp(polVals,'compare')
%             pval = ranksum(cellStats(fN).vals,bkgStats(fN).vals);
            [~,pval] = ttest(cellStats(fN).vals,bkgStats(fN).vals);

            fprintf([objNames{fIdx} ' [%s]: cells (%1.3f CI%1.3f) v bkg (%1.3f CI%1.3f), p=%1.6f t-test, N=%d\n'],...
                polVals,cellStats(fN).mu,cellStats(fN).sd,bkgStats(fN).mu,bkgStats(fN).sd,pval,cellStats(fN).N)
            if pval<alphaVal
                text(fIdx,ax.YLim(2),'*','FontSize',8,'HorizontalAlignment','Center','BackGroundColor','none', 'Tag','asterisk','VerticalAlignment','middle');
            end
        else
%             pvalcell = signrank(cellStats(fN).vals);
%             pvalbkg = signrank(bkgStats(fN).vals);
            [~,pvalcell] = ttest(cellStats(fN).vals);
            [~,pvalbkg] = ttest(bkgStats(fN).vals);

            fprintf([objNames{fIdx} ' [' char(916) 'cells]: %1.3f CI%1.3f p=%1.6f t-test, N=%d\n'],...
                cellStats(fN).mu,cellStats(fN).sd,pvalcell,cellStats(fN).N) 
            if pvalcell<alphaVal
                text(fIdx,ax.YLim(2),'*','FontSize',8,'HorizontalAlignment','Center','BackGroundColor','none', 'Tag','asterisk','VerticalAlignment','middle');
            end
            fprintf([objNames{fIdx} ' [' char(916) 'bkg]: %1.3f CI%1.3f p=%1.6f t-test, N=%d\n'],...
                bkgStats(fN).mu,bkgStats(fN).sd,pvalbkg,bkgStats(fN).N) 
        end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sort input arguments
    function parseInputs
        
        ip = inputParser;
        
        % % Set axes to plot in. Make a new figure by default
        defaultAxes = [];
        addParameter(ip,'ax',defaultAxes)
        
        % % Set values to plot with 'pol' argument:
        % %  'pol','exp'      exp4 values, with pol filter (default)
        % %  'pol','ctrl'      exp2 values, no pol filter
        % %  'pol','compare'    exp4 values minus mean(exp2) value for that line
        validPolValues = {'exp','ctrl','compare'};
        checkPolValues = @(str) any(validatestring(str,validPolValues));
        addParameter(ip,'pol','exp',checkPolValues)
                
        % Get distribution to plot (For cells only): top 10% (default) or
        % all pixels within layerMask
        addParameter(ip,'dist','cells',@(str) any(validatestring(str,{'cells','mask'})) )

        % Get color to use for cell data (default is each object's assigned color:
        % alternative is all the same color). Bkg colors always dark grey
        addParameter(ip,'col','obj',@(str) any(validatestring(str,{'obj','common'})) )
                
        % Allow custom PSI structure to be passed to function. Default
        % action is to load structure on file
        addParameter(ip,'PSI',struct,@isstruct)

        % % Specify objects to plot. Accepts empty strings or strings which don't
        % match an object so groups can be spaced out in the plot or additional
        % labels added
        invalidObjNames = [validPolValues(:);ip.Parameters(:)];
        checkObjNames = @(str) (any(iscell(str) || ischar(str) || isstring(str) ) && all(~contains(str,invalidObjNames)) );
        addOptional(ip,'objNames',{},checkObjNames)
        

        parse(ip,varargin{:})
        
        if any(contains(ip.UsingDefaults,{'objNames', 'PSI'}))
            % One or more inputs not specified
            % Load existing data
            [storedPSI, storedObjNames] = loadPSIstruct;
            
            % If we're plotting everything, skip the pan-tubu driver
            storedObjNames( contains(storedObjNames,'48B06') ) = [];
            storedObjNames( contains(storedObjNames,'R7R8') ) = [];
            
        end
        
        % Use existing data where required
        if isempty(ip.Results.objNames)
            objNames = fullSetObjBrainAreaNames;
        else
            objNames = ip.Results.objNames;
        end
        if ~iscell(objNames)
            objNames = {objNames};
        end
        if size(objNames,2) > size(objNames,1)
            objNames = objNames';
        end
        
        if isempty(fieldnames(ip.Results.PSI))
            PSI = storedPSI;
        else
            PSI = ip.Results.PSI;
        end
                
        % Get plot distribution option
        plotDistribution = ip.Results.dist;

        % Get pol plot type:
        polVals = ip.Results.pol;
        
        % Sort color option
        if strcmp(ip.Results.col,'common')
            useCommonCellColor = 1;
        else
            useCommonCellColor = 0;
        end
        
        % Sort figure axes
        if isempty(ip.Results.ax) || any(contains(ip.UsingDefaults,'ax'))
            figure('color','w')
            ax = gca;
        else
            ax = ip.Results.ax;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Default fieldnames to plot (grouped by brain area)
        function names = fullSetObjBrainAreaNames()
            names = {
                'DmDRA';  ...
                ''; ...
                ''; ...
                'R56F07_AOTU';  ...
                'R73C04_AOTU';  ...
                ''; ...
                ''; ...
                'R17F12_AOTU';  ...
                ''; ...
                ''; ...
                'R49E09_AOTU';  ...
                'R88A06_AOTU';  ...
                'R34H10_AOTU';  ...
                ''; ...
                ''; ...
                'R49E09_Bu';  ...
                'R88A06_Bu';  ...
                'R88A06_Bu_ant';  ...
                'R34H10_Bu';  ...
                ''; ...
                ''; ...
                'R19C08_Bu';  ...
                'R34D03_Bu';  ...
                ''; ...
                ''; ...
                'R34D03_EB';  ...
                ''; ...
                ''; ...
                'SS00096_PB'; ...
                };
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function params = getAVPplotParamsWrapper()
% This just runs a script and returns some of the variables. Since we have
% nested funcs / scoped variables in the main function workspace, we'd have
% to initialize every one of the variables created in the script before
% running it, and if more were added in the future we'd have to come back
% and also add them here. This is a work around that avoids having to do
% that.
getAVPplotParams

params.objCols = objCols;
params.lightGreyCol = lightGreyCol;
params.darkGreyCol = darkGreyCol;
params.defaultLineWidth = defaultLineWidth;
params.axisLabelFontSize = axisLabelFontSize;
params.defaultAxisHeight_cm = defaultAxisHeight_cm;
params.defaultMarkerSize = defaultMarkerSize;

end

