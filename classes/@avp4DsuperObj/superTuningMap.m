function superTuningMap(obj,layerSpec,discretize,darkMode)
if nargin<4 || isempty(darkMode) || ~any(darkMode==0 || darkMode==1)
    darkMode = 0;
end
if nargin<3 || isempty(discretize)
    discretize = 0;
    % if set to 1, angle tunings will be presented with closest of 6 colors, making it easier to read preferences
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
    layerVec = layerSpec;
    numLayers = length(layerVec);
    mipFlag = 0;
end

if numLayers >=5
    numxplots = 5;
else
    numxplots = numLayers;
end

f = figure('color','w');

nidx=0;
for oidx = layerVec
    if mipFlag
        objStr = 'MIP';
    else
        objStr = 'Layers';
    end
%     if isempty(obj.(objStr)(oidx).avgPolImg)
%         getPolMaps(obj.(objStr)(oidx))
%         obj.(objStr)(oidx).Frames =[];
%         obj.(objStr)(oidx).Daq = [];
%     end
    nidx = nidx+1;
    subplot(1+floor(numLayers/5),numxplots,nidx)
    
   if isempty(obj.(objStr)(oidx).AverageFrame)
    mspState = obj.(objStr)(oidx).UseMSP;
    unattendedState = obj.(objStr)(oidx).Unattended;
    
    obj.(objStr)(oidx).UseMSP = 0;
    obj.(objStr)(oidx).Unattended = 1;
    
    getFrames(obj.(objStr)(oidx))
    
    obj.(objStr)(oidx).UseMSP = mspState;
    obj.(objStr)(oidx).Unattended = unattendedState;
    
    obj.(objStr)(oidx).Frames = [];
end


    
lim_bright = max(obj.(objStr)(oidx).AverageFrame(obj.(objStr)(oidx).cellMask));
lim_dark = mean(obj.(objStr)(oidx).AverageFrame(obj.(objStr)(oidx).bkgMask));

avgFrame = ( (obj.(objStr)(oidx).AverageFrame - lim_dark) ./ max( obj.(objStr)(oidx).AverageFrame(:) - lim_dark )  );
avgFrame(avgFrame>1) = 1;
avgFrame(avgFrame<0) = 0;
% lim_bright = max(obj.(objStr)(oidx).InactivityFrame(obj.(objStr)(oidx).cellMask));
% lim_dark = mean(obj.(objStr)(oidx).InactivityFrame(obj.(objStr)(oidx).bkgMask));
% 
% avgFrame = ( (obj.(objStr)(oidx).InactivityFrame - lim_dark) ./ max( obj.(objStr)(oidx).InactivityFrame(:) - lim_dark )  );
% avgFrame(avgFrame>1) = 1;
% avgFrame(avgFrame<0) = 0;

%     avgFrame = sqrt(obj.(objStr)(oidx).avgPolImg./max(obj.(objStr)(oidx).avgPolImg(:)));
    
    tuningFrame = obj.(objStr)(oidx).polPix;
    
    
    
    combFrame = zeros(size(avgFrame));
    if discretize
        combFrame(~isnan(tuningFrame)) = obj.(objStr)(oidx).polPix(~isnan(tuningFrame))./180;
                combFrame(isnan(tuningFrame)) = -1.*avgFrame(isnan(tuningFrame));

    else
        combFrame(~isnan(tuningFrame)) = tuningFrame(~isnan(tuningFrame));  
        combFrame(isnan(tuningFrame)) = -180.*avgFrame(isnan(tuningFrame));
    end
  
    
    imagesc(combFrame)
    
    ax  = gca;
    pbaspect([1,1,1])
    
    % colormap([gray(180);[1 1 1];(colorcet('R2','N',180, 'shift', 0.25))])
    % colormap([gray(180);[1 1 1];equalisecolourmap('RGB',flipud(hsv(180)),'CIE76',[1 1 0],0,1)])
    if darkMode
            cmap = ([flipud(gray(180));[1 1 1];flipud(hsv(180))]);
             if discretize
            cmap = ([flipud(gray(180));[1 1 1];flipud(hsv(180))]);
    end
    else
    cmap = ([gray(180);[1 1 1];flipud(hsv(180))]);
    if discretize
            cmap = ([gray(180);[1 1 1];flipud(hsv(6))]);
    end
    end
    colormap(ax,cmap);
    set(ax,'CLim',[-180 180])
    if discretize
            set(ax,'CLim',[-1 1])

    end
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    if mipFlag
        title(ax,'MIP','FontSize',8)
    else
        title(ax,['Layer ' num2str(oidx)],'FontSize',8)
    end
    
       indFrame = combFrame + 180;
       indFrame = indFrame./360;
       indFrame = round(size(cmap,1).*indFrame);
%        obj.(objStr)(oidx).ActivityFrame9 = ind2rgb(indFrame,cmap);


end
bentitle([obj.Cell ' ' obj.Line ' ' obj.Area ' | ' obj.Name]);
addExportFigToolbar(gcf)
end