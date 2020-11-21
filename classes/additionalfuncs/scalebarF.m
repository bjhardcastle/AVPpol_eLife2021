function scalebarF(ax)
if ~exist('ax', 'var'), ax = gca; end
holdstate = ishold;
hold(ax,'on')

% add a line at y = 0 since the y-axis will no longer convey that
% information
p = plot(ax,[ax.XLim],[0 0],':','Color',[0.5 0.5 0.5],'LineWidth',ax.LineWidth);
uistack(p,'bottom')

% only set offset axes once per axes
% if ~event.hasListener(ax, 'MarkedClean')
addlistener (ax, 'MarkedClean', @(obj,event)yScalebar(ax));
% end
% thanks to Pierre Morel, undocumented Matlab
% and https://stackoverflow.com/questions/38255048/separating-axes-from-plot-area-in-matlab
%
% by Anne Urai, 2016
if ~holdstate
    hold(ax,'off')
end
end

function yScalebar ( ax )
warning('off','MATLAB:callback:error')

% decide how big to make the dF/F scalebar: try to have the bottom aligned
% with 0, then use either 1 or 0.5.
if max(ax.YLim)>=1 && min(ax.YLim)<=0
    ax.YTick = [0,0.2,1];
elseif  max(ax.YLim)>=0.5 && min(ax.YLim)<=0
    ax.YTick = [0,0.25,0.5];
    % elseif max(ax.YLim)>=0 && min(ax.YLim)<=0 ...
    %         && max(ax.YTick)>=0 && min(ax.YTick)<=0
    %     ax.YTick = [min(ax.YLim), 0, max(ax.YLim)];
else
    scaleLength = 0.5*round(range(ax.YLim),2,'significant');
    ax.YTick = [min(ax.YLim),min(ax.YLim)+0.5*scaleLength, min(ax.YLim)+scaleLength];
end

% % repeat for Y (set 2nd row)
ax.YRuler.LineWidth = ax.LineWidth*2;
ax.YRuler.TickLength = [0 0];
ax.YAxis.Label.String = '';
ax.YTickLabelRotation = 0;



%%%% https://www.mathworks.com/matlabcentral/answers/101922-how-do-i-create-a-multi-line-tick-label-for-a-figure-using-matlab-7-10-r2010a#answer_445691
%%%% attributed to Adam Danz, June 2020

% Define each row of labels.
row1 = {'' num2str(range(get(ax, 'Ytick'))) ''};
row2 = {'' ['\Delta' 'F/F'] ''};
% Combine the rows of labels into a cell array; convert non-strings to strings/character vectors.
% labelArray is an nxm cell array with n-rows of m-tick-lables.
labelArray = [row1; row2];
% Combine the rows of labels into individual tick labels
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% tickLabels = strsplit(tickLabels);  % Optional
% Assign ticks and labels
ax.YAxis.TickLabel = tickLabels;
% ax.YAxis.TickLabels = {'';[num2str(max(get(ax, 'Ytick'))) ' ' char(916) 'F/F'];''};
% pause(0.01)
drawnow

ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
ax.YRuler.Axle.VertexData(2,2) = max(get(ax, 'Ytick'));
% drawnow

warning('on','MATLAB:callback:error')
end