function plotCombPolImg(obj,showPolCols,ROIangs,noMask)
%PLOTCOMBPOLIMG Plot avg frame with tuning colors and/or ROIs overlaid
% [] = plotCombPolImg(obj,showPolCols,ROIangs)


% run on superobj
% default behavior is to plot b&w avg img with pol tuning cols overlaid.
% if ROIangs are specified, they will be drawn instead.
if nargin < 3 || isempty(ROIangs)
    ROIangs = [];
% else
%     showPolCols = 0;
end
if nargin < 2
    showPolCols = 1;
    ROIangs = [];
elseif isempty(showPolCols)
    showPolCols = 1;
end
if nargin<4 || isempty(noMask)
    noMask = 0;
end
discretize = 0;

if isempty(obj.polROI)
    loadPolROIs(obj)
end

titleStr = [obj.File];
f = figure('Name',titleStr);
ax(1) = subplot(1,1,1);

loadLayerMasks(obj)

avgFrame = 1-(obj.avgPolImg./max(obj.avgPolImg(:)));

if showPolCols
    tuningFrame = obj.polPix;
    polAngImg = obj.polTuningImg;
    if noMask
        tuningFrame = polAngImg;
    end
    tuningFrame(obj.polSelImg<obj.polSelThreshold) = nan;
    combFrame = zeros(size(avgFrame)    );
    if discretize
        combFrame(~isnan(tuningFrame)) = polAngImg(~isnan(tuningFrame));
    else
        combFrame(~isnan(tuningFrame)) = obj.polTuningImg(~isnan(tuningFrame));
    end
    combFrame(isnan(tuningFrame)) = -180.*avgFrame(isnan(tuningFrame));
    layer_im = imagesc(imrotate(combFrame,90), 'Parent', ax(1));
    colormap(ax(1),[flipud(gray(180));[1 1 1];flipud(hsv(180))])
    set(ax(1),'CLim',[-180 180])
    
else
    layer_im = imagesc(imrotate(-avgFrame,90), 'Parent', ax(1));

    colormap(ax(1),flipud(gray))
    set(ax(1),'CLim',[-1 0])
end


if ~isempty(ROIangs) && ~isempty(obj.polROI)
ylims = ylim;
greycols = linspace(0,1,length(ROIangs));
    for ridx = 1:length(ROIangs)
        ROIang = ROIangs(ridx);
        if ~isempty(find([obj.polROI.angle]==ROIang,1,'first'))
            ROIidx = find([obj.polROI.angle]==ROIang,1,'first');
            hold(ax(1),'on')
            % make line ends connect:
            pos1 = [obj.polROI(ROIidx).position(:,1); obj.polROI(ROIidx).position(1,1)];
            pos2 = [obj.polROI(ROIidx).position(:,2); obj.polROI(ROIidx).position(1,2)];
            if showPolCols 
                %roicol = greycols(ridx).*[1 1 1];
                roicol = 0.3.*[1 1 1];
            else
               roicol = obj.polROI(ROIidx).color;
            end
            plot(pos2,ylims(2)-pos1,'Color',roicol ,'LineWidth',0.5)
        end
    end
end
addScalebar(obj,ax(1),10)
axis(ax(1),'image')
axis(ax(1),'off')
getAVPplotParams
setAVPaxes(ax(1),defaultImageHeight_cm)
tightfig(f)
addExportFigToolbar(f)