function [msebStruct,ax] = pdfplotPSI(varargin)

% Initialize input arguments
ax = [];
PSI = struct();
objNames = {};
plotDistribution = [];
mainLineType = '';
additionalLineType = '';
plotStyleVal = '';
useCommonCellColor = [];
parseInputs


% Get common plotting values
params = getAVPplotParamsWrapper();
objCols = params.objCols;
lightGreyCol = params.lightGreyCol;
darkGreyCol = params.darkGreyCol;
defaultLineWidth = params.defaultLineWidth;
axisLabelFontSize = params.axisLabelFontSize;
defaultAxisHeight_cm = params.defaultAxisHeight_cm;


objSel = ismember(objNames(:),fieldnames(PSI));
assert( sum(objSel)>0 ,'No valid object names entered. Check for missing "_BU" or typos')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot main lines for each object specified

for mIdx = 1:size(objNames,1)
    if objSel(mIdx)
        
        % get distribution (subset or whole mask)
        if strcmp(plotDistribution,'cells')
            distStr = 'Cell';
        elseif strcmp(plotDistribution,'mask')
            distStr = 'LayerMask';
        end
        
        % fetch data
        switch mainLineType
            case 'exp' % default. Exp4 values
                mData{mIdx} = [PSI.(objNames{mIdx}).(['exp' distStr 'Cts'])];
            case 'ctrl' % Exp2 values
                mData{mIdx} = [PSI.(objNames{mIdx}).(['ctrl' distStr 'Cts'])];
            case 'compare' % Exp4 - mean(Exp2)
                mData{mIdx} = [PSI.(objNames{mIdx}).(['exp' distStr 'CtsRelative'])];
        end
        
        
        % calculate parameters then plot line for this object:
        mainData(mIdx).mean = nanmean(mData{mIdx},2)';
        mainData(mIdx).N = sum(~isnan(mean(mData{mIdx},1)));
        mainData(mIdx).SEM = nanstd(mData{mIdx},[],2)'./sqrt( mainData(mIdx).N);
         if isfield(PSI.(objNames{mIdx}),'HistBinEdges')
            if strcmp('additionalLineType','compare')
                histBinEdges = PSI.(objNames{mIdx}).HistBinEdgesRelative;
            else
                histBinEdges = PSI.(objNames{mIdx}).HistBinEdges;
            end
        else
             if strcmp('additionalLineType','compareBkg')
                histBinEdges = linspace(-1,1,size(mainData(mIdx).mean,2)+1);
            else
                histBinEdges = linspace(0,1,size(mainData(mIdx).mean,2)+1);
            end
         end
         histBinCenters = histBinEdges(1:end-1) + 0.5*mode(diff(histBinEdges));
         mainData(mIdx).bincenters = histBinCenters;
         mainData(mIdx).binedges = histBinEdges;
         mainData(mIdx).counts = mData{mIdx}';
        
        if ~useCommonCellColor && any( contains(fieldnames(objCols),(objNames{mIdx})) )
            mainlineprops.col = { objCols.(objNames{mIdx}) };
        else
            mainlineprops.col = { ROIcolor(2) };
        end
        
        switch plotStyleVal
            case 'line'
                mainlineprops.width = defaultLineWidth;
                mainHandles(mIdx) = mseb(mainData(mIdx).bincenters,mainData(mIdx).mean,mainData(mIdx).SEM,mainlineprops,1);
                mainHandles(mIdx).patch.FaceColor = darkGreyCol;
                mainHandles(mIdx).patch.FaceAlpha = 0.15;
            case 'hist'
                mainlineprops.width =0.1;
                h = hseb(mainData(mIdx).bincenters,mainData(mIdx).mean,mainData(mIdx).SEM,mainlineprops,1);
                mainHandles(mIdx).mainLine = h.mainLine;
                mainHandles(mIdx).histPatch = h.mainPatch;
                mainHandles(mIdx).errPatch = h.errPatch;
                mainHandles(mIdx).errPatch.FaceColor = lightGreyCol;
        end
        
    end
end

hold(ax,'on')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot additional lines for each object specified

for aIdx = 1:size(objNames,1)
    if objSel(aIdx)
        
        if strcmp(plotDistribution,'cells')
            distStr = 'Cell';
        elseif strcmp(plotDistribution,'mask')
            distStr = 'LayerMask';
        end
        
        % fetch data
        switch additionalLineType
            case 'none' % default. Exp4 values
                addData{aIdx} = nan(2,50);
            case 'expBkg' % Exp4 background values
                addData{aIdx} = [PSI.(objNames{aIdx}).expBkgCts];
            case 'ctrlBkg' % Exp2 background values
                addData{aIdx} = [PSI.(objNames{aIdx}).ctrlBkgCts];
            case 'compareBkg' % Exp4 - mean(Exp2)
                addData{aIdx} = [PSI.(objNames{aIdx}).expBkgCtsRelative];
            case 'ctrlCell' % Exp2 cell values
                addData{aIdx} = [PSI.(objNames{aIdx}).(['ctrl' distStr 'Cts'])];
                
        end
        
        
        % calculate parameters then plot line for this object:
        additionalData(aIdx).mean = nanmean(addData{aIdx},2)';
        additionalData(aIdx).N = sum(~isnan(mean(addData{aIdx},1)));
        additionalData(aIdx).SEM = nanstd(addData{aIdx},[],2)'./sqrt( additionalData(aIdx).N);
        if isfield(PSI.(objNames{aIdx}),'HistBinEdges')
            if strcmp('additionalLineType','compareBkg')
                histBinEdges = PSI.(objNames{aIdx}).HistBinEdgesRelative;
            else
                histBinEdges = PSI.(objNames{aIdx}).HistBinEdges;
            end
        else
             if strcmp('additionalLineType','compareBkg')
                histBinEdges = linspace(-1,1,size(additionalData(aIdx).mean,2)+1);
            else
                histBinEdges = linspace(0,1,size(additionalData(aIdx).mean,2)+1);
            end
        end
        histBinCenters = histBinEdges(1:end-1) + 0.5*mode(diff(histBinEdges));
        additionalData(aIdx).bincenters = histBinCenters;
        additionalData(aIdx).binedges = histBinEdges;
        additionalData(aIdx).counts = addData{aIdx}';
        
        
        addlineprops.col = {darkGreyCol};
        
        switch plotStyleVal
            case 'line'
                addlineprops.width =defaultLineWidth;
                additionalHandles(aIdx) = mseb(additionalData(aIdx).bincenters, additionalData(aIdx).mean,additionalData(aIdx).SEM,addlineprops,1);
                additionalHandles(aIdx).patch.FaceColor = darkGreyCol;
                additionalHandles(aIdx).patch.FaceAlpha = 0.15;
            case 'hist'
                addlineprops.width =0.1;
                h = hseb(additionalData(aIdx).bincenters, additionalData(aIdx).mean,additionalData(aIdx).SEM,addlineprops,1);
                additionalHandles(aIdx).mainLine = h.mainLine;
                additionalHandles(aIdx).histPatch = h.mainPatch;
                additionalHandles(aIdx).errPatch = h.errPatch;
                additionalHandles(aIdx).errPatch.FaceColor = lightGreyCol*0.8;
                
        end
        
        
    end
end

hold(ax,'on')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Adjust plot

switch plotStyleVal
    case 'line'
        uistack([additionalHandles(:).mainLine],'top')
        uistack([mainHandles(:).mainLine],'top')
    case 'hist'
        uistack([mainHandles(:).histPatch],'top')
        uistack([mainHandles(:).mainLine],'top')
        uistack([additionalHandles(:).errPatch],'bottom')
        uistack([additionalHandles(:).histPatch],'bottom')
end

objNamesNotEmpty = objNames(~strcmp('',objNames(:)));
[L] = legend([mainHandles.mainLine],objNamesNotEmpty);
L.Box = 'off';
L.Interpreter = 'tex';
counter = 0;
for idx = 1:length(mainHandles)
    if ~isempty(mainHandles(idx).mainLine)
        counter = counter+1;
        L.String{counter} = ['\fontsize{' num2str(axisLabelFontSize) '}\color[rgb]{' num2str(mainHandles(idx).mainLine.Color) '}\bf{' strrep(objNamesNotEmpty{counter},'_',' ') '}'];
    end
end
L.Location = 'bestoutside';



switch mainLineType
    
    case{'exp','ctrl'}
        ax.XLim = [0 1];
        ax.XTick = [0:0.2:1];
        ax.XTickLabels ={'0','','','','','1'};
        
        ax.XLabel.String = 'PSI';
        
    case 'compare'
        ax.XLim = [-0.2 0.8];
        ax.XTick = [-0.2:0.2:0.8];
        ax.XLabel.String = ['PSI'];
end

ax.YLabel.String = 'probability';
ax.YLabel.Rotation = 90;
ax.YTick = [0:0.02:0.1];
ax.YTickLabels ={'0','','','','','0.1','','','','','0.2'};
ax.YLim = [0 0.136];
ax.YLabel.Position(2) = 0.05;
trimAxesToLims(ax)
ax.Clipping = 'off';

daspect([range(ax.XTick)/range(ax.YTick),1,1])
% pbaspect([1,1,1])
dratio = daspect(ax);
% for boxplotpsi we use defaultAxisHeight_cm for the print size. Here we'd
% like the y-axes to be the same size on paper as they'll often sit next to each other. 
% (although data may exceed this range since clipping is off)
setAVPaxes(ax,defaultAxisHeight_cm*range(ax.YLim)*dratio(1))

tightfig(gcf)
addExportFigToolbar(gcf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup structure for output

msebStruct = struct;
assert(size(objNames,1) == size(mainData,2),'Mismatch between plot and data: cannot assign outputs')
fN=0;
for fIdx= 1:size(objNames,1)
    if objSel(fIdx)
        fN=fN+1;
        msebStruct.(objNames{fIdx}).mainData = mainData(fN);
        msebStruct.(objNames{fIdx}).mainHandles = mainHandles(fN);
        msebStruct.(objNames{fIdx}).additionalData = additionalData(fN);
        msebStruct.(objNames{fIdx}).additionalHandles = additionalHandles(fN);
    end
end


% Only nested functions beyond here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sort input arguments
    function parseInputs
        
        ip = inputParser;
        
        % % Set axes to plot in. Make a new figure by default
        defaultAxes = [];
        addParameter(ip,'ax',defaultAxes)
        
        % % Set values to plot with 'pol' argument:
        % %  'pol','exp'        (default) exp4 values, with pol filter
        % %  'pol','ctrl'       exp2 values, no pol filter
        % %  'pol','compare'    exp4 values minus mean(exp2) value for that line
        validPolValues = {'exp','ctrl','compare'};
        checkPolValues = @(str) any(validatestring(str,validPolValues));
        addParameter(ip,'pol','exp',checkPolValues)
        
        % % In addition, a second line can be plot for each object:
        % % 'add','none'          (default)
        % % 'add','expBkg'        exp4 background values, with pol filter
        % % 'add','ctrlCell'      exp2 cell values
        % % 'add','ctrlBkg'       exp2 background values, with pol filter
        % % 'add','compareBkg'    bkg corresponding to 'pol','compare'
        validAdditionalValues = {'none','expBkg','ctrlCell','ctrlBkg','compareBkg'};
        checkAdditionalValues = @(str) any(validatestring(str,validAdditionalValues));
        addParameter(ip,'add','expBkg',checkAdditionalValues)
        
        % Get distribution to plot (For cells only): top 10% (default) or
        % all pixels within layerMask
        addParameter(ip,'dist','cells',@(str) any(validatestring(str,{'cells','mask'})) )

        % Get color to use for cell data (default is each object's assigned color:
        % alternative is all the same color). Bkg colors always dark grey
        addParameter(ip,'col','obj',@(str) any(validatestring(str,{'obj','common'})) )
        
        % Get plotting style: line plot (default) or blocky histogram style with
        % filled patch under line
        addParameter(ip,'style','line',@(str) any(validatestring(str,{'line','hist'})) )
     
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

        % Get pol line type:
        mainLineType = ip.Results.pol;
        
        % Get additional line type:
        additionalLineType = ip.Results.add;
        
        % Get plotting style option
        plotStyleVal = ip.Results.style;
        
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
                'R34H10_Bu';  ...
                ''; ...
                ''; ...
                'R19C08_Bu';  ...
                'R78B06_Bu';  ...
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
% running it, and ff more were added in the future we'd have to come back
% and also add them here. This is a work around that avoids having to do
% that.
getAVPplotParams

params.objCols = objCols;
params.lightGreyCol = lightGreyCol;
params.darkGreyCol = darkGreyCol;
params.defaultLineWidth = defaultLineWidth;
params.axisLabelFontSize = axisLabelFontSize;
params.defaultAxisHeight_cm = defaultAxisHeight_cm;
end

%{
function STD = manualSTD(x,nanmeanx,dim)
% built-in matlab nanvar function modified for manually providing mean (not used)  
n = sum(~isnan(x), dim);

% The unbiased estimator: divide by (n-1).  Can't do this when
% n == 0 or 1
denom = n-1;
denom(n == 1) = 1;
denom(n == 0) = 0;

% calc variance
xs = abs(x - nanmeanx).^2;
y = sum(xs, dim, 'omitnan') ./ denom; % abs guarantees a real result
ind = sum(~isnan(xs), dim) < n; % did computation of xs add NaNs
y(ind) = NaN;

STD = sqrt(y);

end
%}