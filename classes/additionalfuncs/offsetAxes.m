function offsetAxes(ax)
% thanks to Pierre Morel, undocumented Matlab
% and https://stackoverflow.com/questions/38255048/separating-axes-from-plot-area-in-matlab
%
% by Anne Urai, 2016
%
% % % Copyright (c) 2016, Anne Urai
% % % All rights reserved.
% % % 
% % % Redistribution and use in source and binary forms, with or without
% % % modification, are permitted provided that the following conditions are met:
% % % 
% % % * Redistributions of source code must retain the above copyright notice, this
% % %   list of conditions and the following disclaimer.
% % % 
% % % * Redistributions in binary form must reproduce the above copyright notice,
% % %   this list of conditions and the following disclaimer in the documentation
% % %   and/or other materials provided with the distribution
% % % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% % % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% % % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% % % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% % % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% % % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% % % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% % % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% % % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% % % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if ~exist('ax', 'var'), ax = gca; end


% only set offset axes once per axes
if ~event.hasListener(ax, 'MarkedClean')
    
    % modify the x and y limits to below the data (by a small amount)
    ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/4;
    ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/4;
    
    % ax.YLim(1)-(ax.YTick(2)-ax.YTick(end))/4;
    
    addlistener (ax, 'MarkedClean', @(obj,event)resetVertex(ax));
end

end

function resetVertex ( ax )
warning('off','MATLAB:callback:error')

% ax.XLim(2) = round(ax.XLim(2), -floor( log10( ax.XTick(2) - ax.XTick(1) ) ) );
% extract the x axis vertex data
% X, Y and Z row of the start and end of the individual axle.
if ~isprop(ax.XRuler.Axle,'VertexData')
    ax.XRuler.Axle.VertexData = [0 0;single(ax.YLim(1)) single(ax.YLim(1));1 1];
end
    ax.XRuler.Axle.VertexData(1,1) = single(min(ax.XTick(:)));
    ax.XRuler.Axle.VertexData(1,2) = single(max(ax.XTick(:)));

% repeat for Y (set 2nd row)
if ~isprop(ax.YRuler.Axle,'VertexData')
    ax.YRuler.Axle.VertexData = [single(ax.XLim(1)) single(ax.XLim(1));0 0;1 1];
end
ax.YRuler.Axle.VertexData(2,1) = single(min(ax.YTick(:)));
ax.YRuler.Axle.VertexData(2,2) = single(max(ax.YTick(:)));

warning('on','MATLAB:callback:error')

end

