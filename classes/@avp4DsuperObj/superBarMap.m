function superBarMap(obj,layerSpec,noMask)
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
   
     if ~noMask && ~isempty(obj.(objStr)(oidx).layerMask) 
           maskFrame = obj.(objStr)(oidx).layerMask.mask;
     elseif ~noMask
           maskFrame = obj.(objStr)(oidx).AverageFrame> ( median(obj.(objStr)(oidx).AverageFrame(:))  );
     elseif noMask
         maskFrame = ones(size(obj.(objStr)(oidx).AverageFrame));
       end

    nidx = nidx+1;
    subplot(1+floor(numLayers/5),numxplots,nidx)
    
    avgFrame = sqrt(obj.(objStr)(oidx).AverageFrame./max(obj.(objStr)(oidx).AverageFrame(:)));
      
    
    barSelframe = obj.(objStr)(oidx).barSelImg;
    if noMask
        if isempty(obj.(objStr)(oidx).barMaxImg)
                getBarMaps(obj.(objStr)(oidx));
            end
            [maxframe,idxframe] = max(obj.(objStr)(oidx).barMaxImg,[],3);
            % convert pref pos idx into approx bar azimuth position on screen
            pos = linspace(-90,90,5);
            posArray = pos(idxframe);
     
            tuningFrame = reshape( posArray , size(idxframe) );
    else
    tuningFrame = obj.(objStr)(oidx).barPix;
    end
    polFrame = obj.(objStr)(oidx).polPix;
    
    
    combFrame = zeros(size(avgFrame));
        combFrame(~isnan(tuningFrame)) = tuningFrame(~isnan(tuningFrame));
 
    combFrame(isnan(tuningFrame)) = -90-180.*avgFrame(isnan(tuningFrame));
        combFrame(maskFrame == 0) = -90-180.*avgFrame(maskFrame == 0);
        combFrame(barSelframe < 0.5) = -90-180.*avgFrame(barSelframe < 0.5);
%         combFrame(polFrame > 0.4) = -90-180.*avgFrame(polFrame > 0.4);

    imagesc(combFrame)
    
    ax  = gca;
    pbaspect([1,1,1])
    
    % colormap([gray(180);[1 1 1];(colorcet('R2','N',180, 'shift', 0.25))])
    % colormap([gray(180);[1 1 1];equalisecolourmap('RGB',flipud(hsv(180)),'CIE76',[1 1 0],0,1)])
    if darkMode
            cmap = ([flipud(gray(180));[1 1 1];flipud(hsv(180))]);
    else
    cmap=([gray(180);[1 1 1];flipud(hsv(180))]);
    end
    colormap(ax,cmap);
      set(ax,'CLim',[-270 90])
    
       indFrame = combFrame + 270;
       indFrame = indFrame./360;
       indFrame = round(size(cmap,1).*indFrame);
%        obj.(objStr)(oidx).ActivityFrame11 = ind2rgb(indFrame,cmap);% 
       
       
       
%     polFrame = obj.(objStr)(oidx).polSelImg;
% %       barFrame([barFrame==0]) = nan;
% %        barFrame(barFrame==1) = nan;
%        if ~noMask && ~isempty(obj.(objStr)(oidx).layerMask) 
%            maskFrame = obj.(objStr)(oidx).layerMask.mask;
%        else
% %            maskFrame = obj.(objStr)(oidx).fftMagImg>0.1;
%            maskFrame = obj.(objStr)(oidx).AverageFrame> ( median(obj.(objStr)(oidx).AverageFrame(:))  );
% 
%        end
%        if ~noMask
%            combFrame = zeros(size(avgFrame));
%            combFrame(~maskFrame==0) = 90.*( polFrame(~maskFrame==0) - barFrame(~maskFrame==0)  );
%            combFrame(maskFrame==0) =  -90-180.*avgFrame(maskFrame==0);
%            
%            imagesc(combFrame)
%            
%            ax  = gca;
%            
%            if darkMode
%                colormap([flipud(gray(180));(polarmap(180))])
%            else
%                colormap([ [1 1 1];gray(180);(polarmap(180))])
%            end
%            
%            set(ax,'CLim',[-270 90])
%            
%        else
%            diffFrame = (barFrame - polFrame).*maskFrame;
%            diffFrame(isnan(diffFrame)) = 0;
%            imagesc(diffFrame)
%            ax  = gca;
%            
%            colormap(ax,polarmap);
%            set(ax,'CLim',[-1 1])
%            
%        end
%        
%        pbaspect([1,1,1])
       
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