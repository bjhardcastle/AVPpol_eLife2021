function Layersummary(obj,layerSelect)
% layerSelect an array defining subset of layers
if ~isempty(obj.Layers)
    if nargin < 2 || isempty(layerSelect)
    layerSelect = 1:length(obj.Layers);
    end
    numLayers = length(layerSelect);

    if isempty(obj.Layers(1).fullAngMaskImg)
        for oidx = 1:length(obj.Layers)
            getPolROIs(obj.Layers(oidx))
            obj.Layers(oidx).Frames = [];
            if oidx < length(obj.Layers)
                obj.Layers(oidx + 1).Daq = obj.Layers(oidx).Daq;
                obj.Layers(oidx).Daq = [];
            else
                obj.Layers(oidx).Daq = [];
            end
        end
    end
    
    % Plot individually
    %{
    for oidx = 1:numLayers 
        
        L = layerSelect(oidx);
        
        axR1(oidx) = subplot(3,numLayers,oidx);
        plotFullAngMaskImg(obj.Layers(L),axR1(oidx))
        colorbar off
        title(num2str(L))
        
        axR2(oidx) = subplot(3,numLayers,oidx + numLayers);
        plotPolSelImg(obj.Layers(L),axR2(oidx))
        colorbar off
        title([])
        
        axR3(oidx) = subplot(3,numLayers,oidx + 2*numLayers);
        plotPolMagImg(obj.Layers(L),axR3(oidx))
        colorbar off
        title([])
        
        
        
    end
    tightfig;
    %}
    
    % Plot as a montage
     s = floor(obj.Layers(1).micronsPerZStep./obj.Layers(1).micronsPerPixel);
    for oidx = 1:numLayers
        L = layerSelect(oidx);
        R1{oidx} = obj.Layers(L).exampleAngMaskImg;%(1:s:end,1:s:end,:);
        R2{oidx} = obj.Layers(L).polSelImg;
        R3{oidx} = obj.Layers(L).fftMagImg;
        R4(:,:,oidx) = obj.Layers(L).k;
    end
    
    figure

    if numLayers < 9
        axR1 = subplot(2,1,1);
        image(cat(2,R1{:}))
        axis(axR1,'off')
        axis image
        
        axR2 = subplot(2,1,2);
        imagesc(cat(2,R2{:}))
        colormap(axR2,flipud(hot))
        axis(axR2,'off')
        axR2.CLim(1) = obj.Layers(1).polSelThreshold;
        axis image
    
    else
       nRows = ceil(numLayers/10);
       for rIdx = 1:nRows
           subplot(2*nRows,1,rIdx);
           if rIdx < nRows
               image(cat(2,R1{(rIdx-1)*10+1:(rIdx-1)*10+10}))
           else
                image(cat(2,R1{(rIdx-1)*10+1:end}))
           end
           axis off
           axis image
           
           subplot(2*nRows,1,rIdx+nRows);
            if rIdx < nRows
                imagesc(cat(2,R2{(rIdx-1)*10+1:(rIdx-1)*10+10}))
            else
                 imagesc(cat(2,R2{(rIdx-1)*10+1:end}))
            end
           colormap(flipud(hot))
           axis off
           clim = get(gca,'CLim');
           set(gca,'CLim',[obj.Layers(1).polSelThreshold clim(2)]);
           axis image
       end
        
    end
        
else
    disp('No Layers objarray exists')
end
end

function plotFullAngMaskImg(obj,ax)
image(obj.fullAngMaskImg,'Parent',ax)
axis(ax,'image')
ax.XTick = [];
ax.YTick = [];
end

function plotPolSelImg(obj,ax)
imagesc(obj.polSelImg,'Parent',ax) % threshold at 10 is just to get rid of edges
colormap(ax,flipud(hot))
title(ax,['selectivity (Avg bins)[thresh>' num2str(obj.polSelThreshold) ']']);
axis(ax,'image');
ax.XTick = []; ax.XTickLabel = {};
ax.YTick = []; ax.YTickLabel = {};

ax.Color = [0 0 0];
colorbar(ax);
ax.CLim(1) = obj.polSelThreshold;
end

function plotPolMagImg(obj,ax)
imagesc(obj.fftMagImg,'Parent',ax)
colormap(ax,flipud(hot))
title(ax,['mag(FFT)'])
title(ax,['magnitude (FFT)[take top' num2str(obj.polMagThreshold) '%]']);
axis(ax,'image');
ax.XTick = []; ax.XTickLabel = {};
ax.YTick = []; ax.YTickLabel = {};

ax.Color = [0 0 0];
colorbar(ax);

end