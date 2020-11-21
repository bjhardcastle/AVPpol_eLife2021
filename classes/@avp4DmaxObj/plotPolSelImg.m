function plotPolSelImg(obj,invertCols,ROIangs)
%PLOTSELPOLIMG Plot pol selecitivity image with option to add ROIs overlaid
% [] = plotPolSelImg(obj,invertCols,ROIangs)

getAVPplotParams

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
        obj.ROI(1).color = lightGreyCol;
        
        manualROIs = 1;
        ROIangs = 1;
    else
        ROIangs = [];
    end
else
    manualROIs = 0;    
end
if nargin < 2 || isempty(invertCols)
    invertCols = 0;
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

layer_im = imagesc(imrotate(obj.polSelImg,90), 'Parent', ax(1));
set(ax(1),'CLim',[0 1])

if invertCols
    colormap(ax(1),flipud(magma))
else
    colormap(ax(1),magma)
end


if ~isempty(ROIangs) && ~isempty(ROI)
    ylims = ylim;
    greycols = linspace(0.2,0.8,length(ROIangs));
    for ridx = 1:length(ROIangs)
        
        ROIang = ROIangs(ridx);
        
        roicol = greycols(ridx).*[1 1 1];

        
        if ~manualROIs && ~isempty(find([ROI.angle]==ROIang,1,'first'))
            ROIidx = find([ROI.angle]==ROIang,1,'first');
        elseif manualROIs
            ROIidx = ROIangs(ridx);
        else
            continue
        end
        hold(ax(1),'on')
            if isequal( size(obj.ROI(ROIidx).position) , [1,4] )  % Ellipse ROI
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
            plot(pos2,ylims(2)-pos1,'Color',roicol ,'LineWidth',1)
        
    end
end
addScalebar(obj,ax(1),10,~invertCols)
axis(ax(1),'image')
axis(ax(1),'off')
getAVPplotParams
setAVPaxes(ax(1),defaultImageHeight_cm)
tightfig(f)
addExportFigToolbar(f)