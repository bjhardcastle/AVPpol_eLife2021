function scalebarT(ax)
if ~exist('ax', 'var'), ax = gca; end
% holdstate = ishold;
% hold(ax,'on')
% 

% only set offset axes once per axes
% if ~event.hasListener(ax, 'MarkedClean')
    addlistener (ax, 'MarkedClean', @(obj,event)xScalebar(ax));
% end
% thanks to Pierre Morel, undocumented Matlab
% and https://stackoverflow.com/questions/38255048/separating-axes-from-plot-area-in-matlab
%
% by Anne Urai, 2016
% if ~holdstate
%     hold(ax,'off')
% end
end

function xScalebar ( ax )
warning('off','MATLAB:callback:error')

% decide how big to make the time scalebar: try to have the left edge aligned
% with 0, then use either 1, 4 or 20.

if max(ax.XLim)>=20 && min(ax.XLim)<=0
    ax.XTick = [0,0.2,20];
elseif max(ax.XLim)>=4 && min(ax.XLim)<=0 
    ax.XTick = [0,0.1,4];
elseif max(ax.XLim)>=1 && min(ax.XLim)<=0 
    ax.XTick = [0,0.5,1];
else
    scaleLength = 0.5*round(range(ax.XLim),2,'significant');
    ax.XTick = [min(ax.XLim),min(ax.XLim)+0.1*scaleLength, min(ax.XLim)+scaleLength];
end

% % repeat for Y (set 2nd row)
ax.XRuler.LineWidth = ax.LineWidth*2;
ax.XRuler.TickLength = [0 0];
ax.XAxis.Label.String = '';
ax.XTickLabelRotation = 0;

ax.XAxis.TickLabel = {''; [num2str(range(get(ax, 'Xtick'))) ' s']; ''};


drawnow

ax.XRuler.Axle.VertexData(1,1) = min(get(ax, 'Xtick'));
ax.XRuler.Axle.VertexData(1,2) = max(get(ax, 'Xtick'));
% drawnow

warning('on','MATLAB:callback:error')
end