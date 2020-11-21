function superSelectMap(obj,layerSpec,darkMode)
if nargin<3 || isempty(darkMode) || ~any(darkMode==0 || darkMode==1)
    darkMode = 0;
end
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
    if isempty(obj.(objStr)(oidx).avgPolImg)
        getPolMaps(obj.(objStr)(oidx))
        obj.(objStr)(oidx).Frames =[];
        obj.(objStr)(oidx).Daq = [];
    end
    nidx = nidx+1;
    subplot(1+floor(numLayers/5),numxplots,nidx)
    
    fftMagFrame = obj.(objStr)(oidx).fftMagImg;
    
    selectivityFrame = obj.(objStr)(oidx).polSelImg;
    
%     combFrame = fftMagFrame.*selectivityFrame;
    combFrame = selectivityFrame;

    imagesc(combFrame)
    
    ax  = gca;
    pbaspect([1,1,1])
    
    % colormap([gray(180);[1 1 1];(colorcet('R2','N',180, 'shift', 0.25))])
    % colormap([gray(180);[1 1 1];equalisecolourmap('RGB',flipud(hsv(180)),'CIE76',[1 1 0],0,1)])
if darkMode 
    cmap = ((magma(256)));
else
   cmap = (flipud(magma(256)));
end
colormap(ax,cmap);
    set(ax,'CLim',[0 1])
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    if mipFlag
        title(ax,'MIP','FontSize',8)
    else
        title(ax,['Layer ' num2str(oidx)],'FontSize',8)
    end
    
    indFrame = combFrame;
       indFrame = round(size(cmap,1).*indFrame);
%        obj.(objStr)(oidx).ActivityFrame8 = ind2rgb(indFrame,cmap);

    
end
bentitle([obj.Cell ' ' obj.Line ' ' obj.Area ' | ' obj.Name]);
addExportFigToolbar(gcf)

end