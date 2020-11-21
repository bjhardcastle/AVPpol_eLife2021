function varargout = plotSinglePolarResp(objarray,roiAngs,applyTuningColorsToManualROIs,addMeanROIresultants,plotCtrl)
if nargin<5 || isempty(plotCtrl)
    polExp=4;
else
    assert((plotCtrl==1|plotCtrl==0),'choose to plot ctrls only: 0 or 1')
    if plotCtrl == 1
        polExp=2;
    else
        polExp=4;
    end
end
if nargin < 2 || isempty(roiAngs)
    roiAngs = [90]; % angle as recorded
    manualROIs = 0;
    
elseif all(roiAngs<15) && all(roiAngs>0)
    manualROIs = 1;
else
    manualROIs = 0;
    
end
if nargin < 3 || isempty(roiAngs)
    applyTuningColorsToManualROIs = 0;
end
if nargin < 4 || isempty(addMeanROIresultants)
    addMeanROIresultants = 1;
end
plotIndividROIresultants = 0;
getAVPplotParams
% f = figure;
f = figure('position', [ 300   502   150  150],'menubar','none','toolbar','none');


if ~manualROIs
    superUseMSP(objarray,1)
end
angIdx = 1;
plotangs = 1;
ax(angIdx) = subplot(1,length(plotangs),angIdx);
hold on
%         ang = plotangs(angIdx);
meanAnglesOut = [];
for roiIdx = 1:length(roiAngs)
    alldata = [];
    nidx=0;
    inclObjIdx = [];
    ROImeanResultant=[];
    ROImeanTuningAng=[];
    
    ROIang = roiAngs(roiIdx);
    
    for oidx = 1:length(objarray)
        
        switch polExp
            % First check each object contains pol flash exp (3) and a pol map exp
            case 4
                % (4,8,10)
                if (objarray(oidx).containsPolMapExp < 4) ||  ...
                        ( isfield(objarray(oidx).pSet(4),'polOffBetweenTrials') && objarray(oidx).pSet(4).polOffBetweenTrials == 1  ) ||  ...
                        ( objarray(oidx).pSet(4).polAngleStep ~= 30  )
                    continue
                end
            case 2
                % (2)
                if isempty(objarray(oidx).containsPolMapExp) || objarray(oidx).containsPolMapExp > 2 ||  ...
                        ( isfield(objarray(oidx).pSet(2),'polOffBetweenTrials') && objarray(oidx).pSet(4).polOffBetweenTrials == 1  ) ||  ...
                        ( objarray(oidx).pSet(2).polAngleStep ~= 30  )
                    continue
                end
        end
        if isempty(objarray(oidx).MIP.polSelImg)
            continue
        end

        if ~manualROIs
            
            if isempty(objarray(oidx).MIP.polROI)
                getPolROIs(objarray(oidx).MIP)
                if isempty(objarray(oidx).MIP.polROI) % no rois made
                    continue
                end
            end
        else
            if isempty( objarray(oidx).MIP.ROI )
                loadROIs(objarray(oidx).MIP)
            end
        end
        if ~manualROIs && isempty(find([objarray(oidx).MIP.polROI.angle]==ROIang,1,'first'))
            continue
        else
            inclObjIdx(end+1) = oidx;
            nidx = nidx+1;
            objarray(oidx).MIP.Fly = objarray(oidx).Fly;
            
            if ~manualROIs
                polROIidx = find([objarray(oidx).MIP.polROI.angle]==ROIang,1,'first');
                data = mean(objarray(oidx).MIP.polAngResp(objarray(oidx).MIP.polROI(polROIidx).mask(:),:),1);
                ROIresp = mean(objarray(oidx).MIP.polAngResp(objarray(oidx).MIP.polROI(polROIidx).mask,:),1);
                
            else
                data = mean(objarray(oidx).MIP.polAngResp(objarray(oidx).MIP.ROI(ROIang).mask(:),:),1);
                ROIresp = mean(objarray(oidx).MIP.polAngResp(objarray(oidx).MIP.ROI(ROIang).mask,:),1);
            end
            
            
            % Find circular mean for each pixel:
            polAngsOrig = [30:30:360];
            [actualROIVec,axtheta] = circ_rangle(circ_axial(deg2rad(repmat(polAngsOrig(1:length(polAngsOrig)),size(ROIresp,1),1)),2), ROIresp, deg2rad(mode(diff(polAngsOrig))), 2);
            theta = axtheta/2;
            Ang =0.5.*(wrapTo360((2.*(rad2deg(theta)))));
            ROImeanResultant(end+1) = actualROIVec;
            ROImeanTuningAng(end+1) = Ang; % Expressed in original angles
            
            
            normdata = (data-min(data(:)))./(max(data(:))-min(data(:)));
            alldata(:,1,nidx) = circshift(normdata',-5); % first data to plot is 180(rec) / -90(weir) / hoiz right
           
            if manualROIs
                roicol = objarray(oidx).MIP.ROI(ROIang).color;
            end
        end
      
    end
    clear opt_area opt_lines opt_axes
    
    polcols = flipud(hsv(180));
    if ~manualROIs
        roicol = polcols(ROIang,:);
    end
    
    opt_area.Color = lightGreyCol;
    opt_area.err = 'sem';
    opt_area.FaceAlpha = 0.8;
    
    
    opt_lines.LineWidth = 1;
    opt_lines.Color = roicol; %polcols(ROIang,:);
    opt_axes.Ticks = [0:0.25:1];
    opt_axes.Ticks = [];

    
    %opt_axes.Labels     = {'-90','-120','-150','180','150','120','90','60','30','0','-30','-60'};
    opt_axes.Labels     = {['  -90' char(176)],'','',['180' char(176)],'','',['90' char(176) ' '],'','',['0' char(176)],'',''};
    
    if ~isempty(alldata)
        polygonplot(alldata,opt_axes,opt_lines,opt_area);
        % Add mean tuning pref line in terms of new, published angles:
        pixPolWeight = mean(squeeze(alldata),2)';
        polAngs = wrapTo180(-circshift([30:30:360],-5)-270);
        [r,axtheta] = circ_rangle(circ_axial(deg2rad(polAngs),2), pixPolWeight, deg2rad(mode(diff(polAngs))), 2);
        theta = axtheta/2;
        Ang =0.5.*(wrapTo180((2.*(rad2deg(theta)))));
        hold on
        R  = [0; 1];
        TH = [deg2rad(Ang-90) deg2rad(Ang-270)].*ones(2,2);
        
        [X,Y] = pol2cart(TH, R);
        
        % plot mean resultant line and tuning color
        if addMeanROIresultants
            
            if ~manualROIs
            for n = 1:2,  plot(X(2,:),Y(2,:),'color',polcols(ROIang,:),'LineWidth',2,'LineStyle',':'); end
        elseif ~applyTuningColorsToManualROIs
            for n = 1:2,  plot(X(2,:),Y(2,:),'color',objarray(oidx).MIP.ROI(ROIang).color,'LineWidth',2,'LineStyle',':'); end
        else
            tuningAngle = round(atand(abs(Y(2,2))/X(2,1)));
            if tuningAngle<0
                tuningAngle = tuningAngle + 180;
            elseif tuningAngle == 180
                tuningAngle = 0;
            end
            col = polcols( tuningAngle + 1 ,:) ;
            for n = 1:2,  plot(X(2,:),Y(2,:),'color',col,'LineWidth',2,'LineStyle',':'); end
            
        end
        end
        
        % plot individual recording resultant lines and tuning colors
        if plotIndividROIresultants
            for ridx = 1:length(ROImeanResultant)
                tuningAngle = ROImeanTuningAng(ridx);
                ROItuning = wrapTo180(-ROImeanTuningAng(ridx)-270);
                TH = [deg2rad(ROItuning-90) deg2rad(ROItuning-270)].*ones(2,2);
                R  = [0; ROImeanResultant(ridx)];
                [X,Y] = pol2cart(TH, R);
                if ~manualROIs
                    for n = 1:2,  plot(X(2,:),Y(2,:),'color',polcols(ROIang,:),'LineWidth',0.5,'LineStyle','-'); end
                elseif ~applyTuningColorsToManualROIs
                    for n = 1:2,  plot(X(2,:),Y(2,:),'color',objarray(oidx).MIP.ROI(ROIang).color,'LineWidth',0.5,'LineStyle','-'); end
                else
                    col = polcols( floor(tuningAngle) + 1 ,:) ;
                    for n = 1:2,  plot(X(2,:),Y(2,:),'color',col,'LineWidth',0.5,'LineStyle','-'); end
                    
                end
            end
        end
        
        meanAnglesOut(roiIdx,:) = ROImeanTuningAng;
        resultantVectorOut(roiIdx,:) = ROImeanResultant;
        
    end
    
    
    
    
end

if nargout >0
    varargout{1} = meanAnglesOut;
end
if nargout >1
    varargout{2} = resultantVectorOut;
end

for tidx = 1:length(ax)
    ax(tidx).YAxis.Label.Visible = 'off';
        ax(tidx).XAxis.Label.Visible = 'off';

    for cidx = 1:length(ax(tidx).Children)
        if strcmp(ax(tidx).Children(cidx).Type,'text')
            ax(tidx).Children(cidx).FontSize = axisLabelFontSize;
        end
    end
    setAVPaxes(ax(tidx),2)
    tightfigadv(f);
    
end
addExportFigToolbar(f)
end