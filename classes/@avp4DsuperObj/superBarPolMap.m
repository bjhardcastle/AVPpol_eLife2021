function superBarPolMap(obj,layerSpec,noMask)
if nargin<3 || isempty(noMask) || ~any(noMask==0 || noMask==1)
    noMask = 0;
end
mipFlag = 0;
if nargin<2 || isempty(layerSpec)
    numLayers = length(obj.Layers);
    layerVec = 1:numLayers;
    mipFlag = 0;
elseif layerSpec == 0 % special case for plotting MIP data
    numLayers = 1;
    layerVec = 1;
    mipFlag = 1;
else
    numLayers = 1;
    layerVec = layerSpec;
    mipFlag = 0;
end
darkMode = 0;
if numLayers >=5
    numxplots = 5;
else
    numxplots = numLayers;
end

if  ~any(obj.Layers(1).TrialSeqNum == 5)
    disp('No bar exp exists')
   return 
elseif  ~any(obj.Layers(1).TrialSeqNum == 4)
    disp('No pol exp exists')
   return 

end

f = figure('color','w');

nidx=0;
for oidx = layerVec
    if mipFlag
        objStr = 'MIP';
    else
        objStr = 'Layers';
    end
    if isempty(obj.(objStr)(oidx).barAvgImg)
        getBarMaps(obj.(objStr)(oidx))
    end
     if isempty(obj.(objStr)(oidx).polSelImg)
        getPolMaps(obj.(objStr)(oidx))
     end
    if isempty(obj.(objStr)(oidx).AverageFrame)
        if isempty(obj.(objStr)(oidx).Frames)
            getFrames(obj.(objStr)(oidx))
        end
        obj.(objStr)(oidx).AverageFrame = mean(obj.(objStr)(oidx).Frames,3);
        obj.(objStr)(oidx).Frames = [];
    end
    nidx = nidx+1;
    subplot(1+floor(numLayers/5),numxplots,nidx)
    
    avgFrame = sqrt(obj.(objStr)(oidx).avgPolImg./max(obj.(objStr)(oidx).avgPolImg(:)));
    
    barFrame = obj.(objStr)(oidx).barSelImg;
    polFrame = obj.(objStr)(oidx).polSelImg;
%       barFrame([barFrame==0]) = nan;
%        barFrame(barFrame==1) = nan;
       if ~noMask && ~isempty(obj.(objStr)(oidx).layerMask) 
           maskFrame = obj.(objStr)(oidx).layerMask.mask;
       else
%            maskFrame = obj.(objStr)(oidx).AverageFrame> ( median(obj.(objStr)(oidx).AverageFrame(:))  );
maskFrame = ones(size( obj.(objStr)(oidx).AverageFrame));
       end
       if ~noMask
           combFrame = zeros(size(avgFrame));
           combFrame(~maskFrame==0) = 90.*( polFrame(~maskFrame==0) - barFrame(~maskFrame==0)  );
           combFrame(maskFrame==0) =  -90-180.*avgFrame(maskFrame==0);
           
           imagesc(combFrame)
           
           ax  = gca;
           
           if darkMode
               cmap = [flipud(gray(180));[1 1 1];(polarmap(180))];
           else
               cmap = ([ gray(180);[1 1 1];(polarmap(180))]);
           end
           colormap(ax,cmap);
           set(ax,'CLim',[-270 90])
           
       else
           diffFrame = (polFrame-barFrame).*maskFrame;
           diffFrame(isnan(diffFrame)) = 0;
           imagesc(diffFrame)
           ax  = gca;
           cmap = polarmap;
           colormap(ax,cmap);
           set(ax,'CLim',[-1 1])
           
       end
       
       pbaspect([1,1,1])
       
%        if ~noMask
%                   indFrame = combFrame + 270;       
%                   indFrame = indFrame./360;
%        else
%        indFrame = (combFrame + 1)/2;
%        end
%        indFrame = round(size(cmap,1).*indFrame);
% %        obj.(objStr)(oidx).ActivityFrame12 = ind2rgb(indFrame,cmap);

       ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    if mipFlag
        title(ax,'MIP','FontSize',8)
    else
        title(ax,['Layer ' num2str(oidx)],'FontSize',8)
    end
end
bentitle([obj.Cell ' ' obj.Line ' ' obj.Area ' | ' obj.Name]);
addExportFigToolbar(gcf)
end