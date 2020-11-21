function axesAnatomyLabels(ax)
% thanks to Pierre Morel, undocumented Matlab
% and https://stackoverflow.com/questions/38255048/separating-axes-from-plot-area-in-matlab
%
% by Anne Urai, 2016

if ~exist('ax', 'var'), ax = gca; end


% only set offset axes once per axes
if ~event.hasListener(ax, 'MarkedClean')
%     
%     % modify the x and y limits to below the data (by a small amount)
%     ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/4;
%     ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/4;
%     
%     % ax.YLim(1)-(ax.YTick(2)-ax.YTick(end))/4;    
    
    addlistener (ax, 'MarkedClean', @(obj,event)setLabels(ax));
end

end

function setLabels ( ax )

% ax.XLim(2) = round(ax.XLim(2), -floor( log10( ax.XTick(2) - ax.XTick(1) ) ) );
ax.FontWeight = 'bold';
% extract the x axis vertext data
% X, Y and Z row of the start and end of the individual axle.
ax.XTickLabel =' ';
for x = 2:length(ax.XTick)
                ax.XTickLabel = [ax.XTickLabel; ' '];
end
ax.XTickLabel(1) = 'm';
ax.XTickLabel(end) = 'l';

ax.YTickLabel =' ';
for x = 2:length(ax.YTick)
                ax.YTickLabel = [ax.YTickLabel; ' '];
end
ax.YTickLabel(1) = 'p';
ax.YTickLabel(end) = 'a';

ax.ZTickLabel =' ';
for x = 2:length(ax.ZTick)
                ax.ZTickLabel = [ax.ZTickLabel; ' '];
end
ax.ZTickLabel(1) = 'v';
ax.ZTickLabel(end) = 'd';

ax.XRuler.Axle.VertexData = single(zeros(3,2));
ax.YRuler.Axle.VertexData= single(zeros(3,2));
ax.ZRuler.Axle.VertexData= single(zeros(3,2));
end
