function setAVPaxes(ax,axHeight_cm,axWidth_cm)

getAVPplotParams

if nargin <2   ||  nargin==3 && ( isempty(axHeight_cm) && isempty(axWidth_cm) )
    axHeight_cm = defaultAxisHeight_cm;
elseif nargin == 3 && ( isempty(axHeight_cm) && ~isempty(axWidth_cm) )
    axHeight_cm = [];
end
if nargin < 3
    axWidth_cm = [];
end



axUnitStr = ax.Units;
origDim = [ax.Position(3) ax.Position(4)];
[maxVal,maxDim] = max(origDim);
ax.ActivePositionProperty = 'position';
ax.Units = 'centimeters';
% Change axis width and height according to specified input
% if ~isempty(axWidth_cm)
%     rescale_factorX = ax.Position(4)/ax.Position(3);
% end
% if ~isempty(axHeight_cm)
%     rescale_factorY = ax.Position(3)/ax.Position(4);
% end
% ax.TickLength(1) = axisTickLength./ax.Position(2+maxDim);
if isprop(ax,'PlotBoxAspectRatio')
    pbar = ax.PlotBoxAspectRatio(2)/ax.PlotBoxAspectRatio(1);
else
    pbar = ax.Position(4)/ax.Position(3);
end
if ~isempty(axWidth_cm)
    rescale_factorX = axWidth_cm*pbar;%/origDim(1);
    ax.Position(3) = axWidth_cm;
    if isempty(axHeight_cm)
        ax.Position(4) =rescale_factorX;
    end
end
if ~isempty(axHeight_cm)
    rescale_factorY = axHeight_cm/pbar;%*origDim(2);
    ax.Position(4) = axHeight_cm;
    if isempty(axWidth_cm)
        ax.Position(3) = rescale_factorY;
    end
end

% ax.TickLength(1) = ax.TickLength(1)*(origDim(maxDim)/ax.Position(2+maxDim));
% if strcmp(ax.DataAspectRatioMode,'manual') || strcmp(ax.PlotBoxAspectRatioMode,'manual')
%
%     %%% figure out how to get plotbox size in cm , then scale tick length
%     %%% appropriately..
%     [maxVal,maxDim] = max([ax.Position(3) ax.Position(4)]);
%     lims = xlim;
%     lims(2,:) = ylim;
% ax.TickLength = [1 1].*axisTickLength.*(maxVal/range(lims(maxDim,:)));
% else
%     ax.TickLength = [1 1].*axisTickLength./(max([ax.Position(3) ax.Position(4)]));
% end
% ax.TickLength(1) = axisTickLength./ax.Position(4);
if ~strcmp(ax.Type,'polaraxes') ...
    && (strcmp(ax.DataAspectRatioMode,'manual') || strcmp(ax.PlotBoxAspectRatioMode,'manual')) ...
        && ( ax.DataAspectRatio(1) == ax.DataAspectRatio(2) || ax.PlotBoxAspectRatio(1) == ax.PlotBoxAspectRatio(2) )
    % height and width must be the same, else tick length behavior is
    % controlled by the shortest dimension, not the longest as we assume
    % below
    [minVal,minDim] = min([axWidth_cm axHeight_cm]);
    
    % ax.Position(3:4) = minVal;
    ax.TickLength = [1 1].*axisTickLength./minVal;
    
else
    ax.TickLength = [1 1].*axisTickLength./(max([ax.Position(3) ax.Position(4)]));
end
% ax.Units = axUnitStr;
ax.TickDir = 'out';

ax.Box = 'off';

ax.FontSize = axisLabelFontSize;
ax.FontWeight = 'normal';

if ~strcmp(ax.Type,'polaraxes')
    ax.XAxis.Label.FontSize = axisLabelFontSize;
    ax.YAxis.Label.FontSize = axisLabelFontSize;
    ax.XAxis.Color = [0 0 0];
    ax.YAxis.Color = [0 0 0];
    ax.XAxis.FontName = 'arial';
    ax.YAxis.FontName = 'arial';
else
    ax.ThetaAxis.Label.FontSize = axisLabelFontSize;
    ax.RAxis.Label.FontSize = axisLabelFontSize;
    ax.ThetaAxis.Color = [0 0 0];
    ax.RAxis.Color = [0 0 0];
    ax.ThetaAxis.FontName = 'arial';
    ax.RAxis.FontName = 'arial';
end
ax.Title.FontSize = axisLabelFontSize;
ax.Title.FontWeight = 'normal';
% ax.Title.FontName = 'arial';

set(ax,'layer','top')

