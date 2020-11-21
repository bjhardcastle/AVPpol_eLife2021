function axesAnatomyLabelsR(ax)
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
ticks = {};
for x = 1:length(ax.XTick)
               ticks{x} = ' ';
end
ticks{1} = 'm';
ticks{end} = 'l';
xticklabels(ax,ticks)

ticks = {};
for x = 1:length(ax.YTick)
               ticks{x} = ' ';
end
ticks{1} = 'p';
ticks{end} = 'a';
yticklabels(ax,ticks)

ticks = {};
ax.ZTickLabel =' ';
for x = 1:length(ax.ZTick)
 ticks{x} = ' ';
end
ticks{1} = 'v';
ticks{end} = 'd';
zticklabels(ax,ticks)

   if strcmp(ax.XTickLabel(2),' ')
        ax.XTickLabel{2} = [num2str(ax.XTick(2)) ' \mum'];
    elseif  strcmp(ax.XTickLabel(2),'l')
          ax.XTickLabel{2} = [num2str(ax.XTick(2)) ' \mum      l'];
   end 
    
      if strcmp(ax.YTickLabel(2),' ')
        ax.YTickLabel{2} = [num2str(ax.YTick(2)) ' \mum'];
    elseif  strcmp(ax.YTickLabel(2),'p')
          ax.YTickLabel{2} = [num2str(ax.YTick(2)) ' \mum      p'];
      end 
      if strcmp(ax.ZTickLabel(2),' ')
        ax.ZTickLabel{2} = [num2str(ax.ZTick(2)) ' \mum'];
    elseif  strcmp(ax.ZTickLabel(2),'d')
          ax.ZTickLabel{2} = [num2str(ax.ZTick(2)) ' \mum      d'];
   end 
    

          ax.TickLength = [0 0];  


%     ax.XRuler.Axle.VertexData(1,1) = single(ax.XTick(1));
%     ax.XRuler.Axle.VertexData(1,2) = single(ax.XTick(2));
%         ax.XRuler.Axle.VertexData(2,1) = single(ax.XTick(1));
%     ax.XRuler.Axle.VertexData(2,2) = single(ax.XTick(2));
%      ax.XRuler.Axle.VertexData(3,1) = single(ax.XTick(1));
%     ax.XRuler.Axle.VertexData(3,2) = single(ax.XTick(2));
% % %    ax.XRuler.LineWidth = 2;
%  
%     
%          ax.YRuler.Axle.VertexData(1,1) = 0;
%     ax.YRuler.Axle.VertexData(1,2) = 0;
%         ax.YRuler.Axle.VertexData(2,1) = ax.YTick(1);
%     ax.YRuler.Axle.VertexData(1,2) = ax.YTick(2);
%   ax.YRuler.Axle.VertexData(3,1) = 0;
%     ax.YRuler.Axle.VertexData(3,2) = 0;
% 
%     
%     ax.ZRuler.Axle.VertexData(1,1) = 0;
%     ax.ZRuler.Axle.VertexData(1,2) = 0;
%   ax.ZRuler.Axle.VertexData(2,1) = 0;
%     ax.ZRuler.Axle.VertexData(2,2) = 0;
%       ax.ZRuler.Axle.VertexData(3,1) = ax.ZTick(1);
%     ax.ZRuler.Axle.VertexData(3,2) = ax.ZTick(2);
% 
   
ax.XRuler.Axle.VertexData = single(zeros(3,2));
ax.YRuler.Axle.VertexData= single(zeros(3,2));
ax.ZRuler.Axle.VertexData= single(zeros(3,2));


end
