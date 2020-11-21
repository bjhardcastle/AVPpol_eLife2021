function superPBTMap(obj,layerSpec,noMask)
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

noThrustMap = 0;
noBarMap= 0;
if  ~any(obj.Layers(1).TrialSeqNum == 4)
    disp('No pol exp exists')
    return
elseif  ~any(obj.Layers(1).TrialSeqNum == 5) && ~any(obj.Layers(1).TrialSeqNum == 6)
    disp('No bar or thrust exps exists')
    return
elseif ~any(obj.Layers(1).TrialSeqNum == 6)
    noThrustMap = 1;
 elseif ~any(obj.Layers(1).TrialSeqNum == 5)
    noBarMap = 1;   
end
loadLayerMasks([obj.Layers,obj.MIP]);
f = figure('color','w');

nidx=0;
for oidx = layerVec
    if mipFlag
        objStr = 'MIP';
    else
        objStr = 'Layers';
    end
    if ~noBarMap && isempty(obj.(objStr)(oidx).barAvgImg)
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
      
    % Construct frames for selectivity for one stimulus over the other two
    % combined
    if noThrustMap
        polSelFrame = (obj.(objStr)(oidx).polSelImg - obj.(objStr)(oidx).barSelImg);
        barSelFrame = (obj.(objStr)(oidx).barSelImg - obj.(objStr)(oidx).polSelImg);
        thrustSelFrame =  zeros(size(polSelFrame));
    elseif noBarMap
        thrustSelFrame = (obj.(objStr)(oidx).thrustSelImg - obj.(objStr)(oidx).polSelImg);
        polSelFrame = (obj.(objStr)(oidx).polSelImg - obj.(objStr)(oidx).thrustSelImg);
        barSelFrame = zeros(size(polSelFrame));
    else
        thrustSelFrame = (obj.(objStr)(oidx).thrustSelImg - obj.(objStr)(oidx).barSelImg - obj.(objStr)(oidx).polSelImg);
        polSelFrame = (obj.(objStr)(oidx).polSelImg - obj.(objStr)(oidx).thrustSelImg - obj.(objStr)(oidx).barSelImg);
        barSelFrame = (obj.(objStr)(oidx).barSelImg - obj.(objStr)(oidx).thrustSelImg - obj.(objStr)(oidx).polSelImg);
    end
    % Assign initial colors by converting to 3Darrays ( can remap colors
    % later)    
    % magenta - bar
    % green - pol
    % thrust - blue
    p = cat(3,thrustSelFrame, 0.*thrustSelFrame, 0.*thrustSelFrame) +cat(3,0*polSelFrame, polSelFrame, 0.*polSelFrame) +cat(3, 0.*barSelFrame, 0.*barSelFrame,barSelFrame)  ;

    % We're only interested in regions where selectivity ends up positive:
    combFrame = double(p>0);
%Remap red to magenta
redpix = combFrame(:,:,1)==1;
combFrame(:,:,3) = combFrame(:,:,3) + redpix;

% Fill blank pixels with average intensity image as a background to show
% anatomy
avgbkg = 1-avgFrame;
M = cat(3,avgbkg,avgbkg,avgbkg);
bkg = cat(3,sum(combFrame,3)==0,sum(combFrame,3)==0,sum(combFrame,3)==0);
combFrame(bkg) =M(bkg);

% Apply mask
mask = cat(3,maskFrame,maskFrame,maskFrame);
combFrame(~mask) = M(~mask);

    

    image(imrotate(combFrame,90))
    
    ax  = gca;
    pbaspect([1,1,1])
    
    obj.(objStr)(oidx).ActivityFrame10 = combFrame;
    
    
    
    % colormap([gray(180);[1 1 1];(colorcet('R2','N',180, 'shift', 0.25))])
    % colormap([gray(180);[1 1 1];equalisecolourmap('RGB',flipud(hsv(180)),'CIE76',[1 1 0],0,1)])
%     if darkMode
%             colormap([flipud(gray(180));[1 1 1];flipud(hsv(180))])
%     else
%     colormap([gray(180);[1 1 1];flipud(hsv(180))])
%     end
%       set(ax,'CLim',[-270 90])
    
    
%             tuningFrame = obj.(objStr)(oidx).barPix;
% 
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
addScalebar(obj.Layers(1))
bentitle([obj.Cell ' ' obj.Line ' ' obj.Area ' | ' obj.Name]);
addExportFigToolbar(gcf)
end