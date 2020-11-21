function plotCombPolImgManual(obj,showPolCols,ROIangs,noMask,noiseFilter)
%PLOTCOMBPOLIMG Plot avg frame with tuning colors and/or ROIs overlaid
% [] = plotCombPolImg(obj,showPolCols,ROIangs)


% run on superobj
% default behavior is to plot b&w avg img with pol tuning cols overlaid.
% if ROIangs are specified, they will be drawn instead.
if nargin < 3 || isempty(ROIangs)
    ROIangs = [];
    manualROIs=0;
elseif all(ROIangs<30) && all(ROIangs>0)
    manualROIs = 1;
elseif ROIangs == -1
    % use layermask
    disp('Assigning layer mask as ROI(1)')
    loadLayerMasks(obj)
    if ~isempty(obj.layerMask)
    obj.ROI(1).mask = obj.layerMask.mask;
        obj.ROI(1).position = obj.layerMask.position;
        obj.ROI(1).color = ROIcolor(8);

        manualROIs = 1;
        ROIangs = 1;
    else
            ROIangs = [];

    end
else
        manualROIs = 0;
end
if nargin < 2
    showPolCols = 1;
    ROIidx = [];
elseif isempty(showPolCols)
    showPolCols = 1;
end
if nargin<4 || isempty(noMask)
    noMask = 0;
end
discretize = 0;
if nargin<5 || isempty(noiseFilter)
    noiseFilter = 1;
end

if ~manualROIs
    if isempty(obj.polROI)
        loadPolROIs(obj)
    end
    ROI = obj.polROI;
elseif  manualROIs
    if isempty(obj.ROI)
        loadROIs(obj)
    end
    ROI = obj.ROI;
end



titleStr = [obj.File];
f = figure('Name',titleStr,'color','w');
ax(1) = subplot(1,1,1);

loadLayerMasks(obj)
if isempty(obj.AverageFrame)
    mspState = obj.UseMSP;
    unattendedState = obj.Unattended;
    
    obj.UseMSP = 0;
    obj.Unattended = 1;
    
    getFrames(obj)
    
    obj.UseMSP = mspState;
    obj.Unattended = unattendedState;
    
    obj.Frames = [];
end
if isempty(obj.InactivityFrame)
    getInactivityFrame(obj)
end
    im = obj.InactivityFrame;
lim_bright = max(im(obj.cellMask));
lim_dark = mean(im(obj.bkgMask));
% normalize:
avgFrame = im - lim_dark;
avgFrame = avgFrame./(mean(avgFrame(obj.layerMask.mask(:))) + 3*std(avgFrame(obj.layerMask.mask(:))) );
avgFrame = 1-avgFrame;
% avgFrame = 1-( (im - lim_dark) ./ max(im(:) - lim_bright - lim_dark )  );
avgFrame(avgFrame>1) = 1;
avgFrame(avgFrame<0) = 0;
% lim_bright = max(obj.InactivityFrame(obj.cellMask));
% lim_dark = mean(obj.InactivityFrame(obj.bkgMask));
% 
% avgFrame = 1-( (obj.InactivityFrame - lim_dark) ./ max( obj.InactivityFrame(:) - lim_dark )  );
% avgFrame(avgFrame>1) = 1;
% avgFrame(avgFrame<0) = 0;

if showPolCols
    tuningFrame = obj.polPix;
    polAngImg = obj.polTuningImg;
    if noMask
        tuningFrame = polAngImg;
    end
    tuningFrame(obj.polSelImg<obj.polSelThreshold) = nan;
    
    if noiseFilter

        s{2} = [0 1 0; 0 1 0; 0 0 0];
        s{3} = [0 0 0; 0 1 0; 0 1 0];
        s{4} = [0 0 0; 0 1 1; 0 0 0];
        s{5} = [0 0 0; 1 1 0; 0 0 0];
        s{6} = [1 0 0; 0 1 0; 0 0 0];
        s{7} = [0 0 1; 0 1 0; 0 0 0];
        s{8} = [0 0 0; 0 1 0; 1 0 0];
        s{9} = [0 0 0; 0 1 0; 0 0 1];
        s{end+1} = [0 0 0; 0 1 0; 0 0 0];
        
        for n = 2
            for sIdx = 2:length(s)
                
                % Hit or miss removal of isolated single pixels in center of 3x3 sq
                notNan = ~isnan(tuningFrame);
                a = imerode(notNan,s{sIdx});
                b = imerode(~notNan,~s{sIdx});
                
                no_islands = (a&b);
                tuningFrame(no_islands) = nan;
                
           
            end
        end
    end
    
    combFrame = zeros(size(avgFrame)    );
    if discretize
        combFrame(~isnan(tuningFrame)) = polAngImg(~isnan(tuningFrame));
    else
        combFrame(~isnan(tuningFrame)) = obj.polTuningImg(~isnan(tuningFrame));
    end
    combFrame(isnan(tuningFrame)) = -180.*avgFrame(isnan(tuningFrame));
    layer_im = imagesc(imrotate(combFrame,90), 'Parent', ax(1));
    colormap(ax(1),[flipud(gray(181));[0 0 0];flipud(hsv(180))])
    set(ax(1),'CLim',[-180 180])
    
else
    layer_im = imagesc(imrotate(-avgFrame,90), 'Parent', ax(1));

    colormap(ax(1),flipud(gray))
    set(ax(1),'CLim',[-1 0])
end


if ~isempty(ROIangs) && ~isempty(ROI)
ylims = ylim;
greycols = linspace(0,1,length(ROIangs));
    for ridx = 1:length(ROIangs)
        ROIang = ROIangs(ridx);

        if ~manualROIs && ~isempty(find([ROI.angle]==ROIang,1,'first'))
            ROIidx = find([ROI.angle]==ROIang,1,'first');
        elseif manualROIs
            ROIidx = ROIangs(ridx);
        else
            continue
        end

            hold(ax(1),'on')
            
            if showPolCols
                %roicol = greycols(ridx).*[1 1 1];
                roicol = 0.3.*[1 1 1];
            else
                roicol = ROI(ROIidx).color;
            end
            
            if isequal( size(ROI(ROIidx).position) , [1,4] )  % Ellipse ROI
                % [ xmin ymin width height ]
                xmin = ROI(ROIidx).position(1);
                ymin = ROI(ROIidx).position(2);
                width = ROI(ROIidx).position(3);
                height = ROI(ROIidx).position(4);
                
                % ROIs were drawn as close to circles as possible, so
                % approximate as a circular path
                r = 0.5*max([width height]);
                xcenter = xmin + 0.5*width;
                ycenter = ymin + 0.5*width;
                
                pos1 = r*cos(linspace(0,2*pi,20)) + xcenter;
                pos2 = r*sin(linspace(0,2*pi,20)) + ycenter;
                
                
                
            else % Polygon
                % make line ends connect:
                pos1 = [ROI(ROIidx).position(:,1); ROI(ROIidx).position(1,1)];
                pos2 = [ROI(ROIidx).position(:,2); ROI(ROIidx).position(1,2)];
            
        end
            % plot rotated 90 deg
            plot(pos2,ylims(2)-pos1,'Color',roicol ,'LineWidth',0.5)

    end
end
addScalebar(obj,ax(1),10)
axis(ax(1),'image')
axis(ax(1),'off')
getAVPplotParams
setAVPaxes(ax(1),defaultImageHeight_cm)
tightfig(f)
addExportFigToolbar(f)