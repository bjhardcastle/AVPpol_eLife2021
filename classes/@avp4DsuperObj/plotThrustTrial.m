function plotThrustTrial(objarray,roiAngs,errorbar,numCycles)
manualROIs = 0;
if nargin < 2 || isempty(roiAngs)
    roiAngs = [30 90 120]; % angle as recorded
    
elseif all(roiAngs<15) && all(roiAngs>0)
    manualROIs = 1;
end
if nargin < 3  || isempty(errorbar)
    errorbar = 1;
end
if nargin<4 || isempty(numCycles)
    numCycles = [1,4];
elseif length(numCycles) == 1
    numCycles = [1,numCycles];
end
plotangs = [90];

addMotorPositionIndicator = 0;
addMotorVelocityIndicator = 0;
numAxes = 1+addMotorPositionIndicator;
f=figure;
% for  oidx = 1:length(objarray)
%     if objarray(oidx).containsPolMapExp > 2
%         getPolMaps(objarray(oidx).MIP)
%         getPolROIs(objarray(oidx).MIP)
%     end
% end
superUseMSP(objarray)

ax(1) = subplot(numAxes,1,1);
hold on
superInclObjIdx = [];
for angIdx = 1:length(plotangs)
    ang = plotangs(angIdx);
    
    for roiIdx = 1:length(roiAngs)
        ROIang = roiAngs(roiIdx);
        
        inclObjIdx = [];
        
        for oidx = 1:length(objarray)
            
            
            % First check each object contains pol flash exp (3) and a pol map exp
            % (4,8,10)
            if ~any(objarray(oidx).Exps == 6)
                continue
            end
            
            if ~manualROIs
                loadPolROIs(objarray(oidx).MIP)
                
                if isempty(objarray(oidx).MIP.polROI)
                    getPolROIs(objarray(oidx).MIP)
                    if isempty(objarray(oidx).MIP.polROI) % no rois made
                        continue
                    end
                end
                
                if isempty(find([objarray(oidx).MIP.polROI.angle]==ROIang,1,'first'))
                    continue
                else
                    inclObjIdx(end+1) = oidx;
                    
                    % assign this polROI to ROI(1)
                    objarray(oidx).MIP.ROI = objarray(oidx).MIP.polROI(find([objarray(oidx).MIP.polROI.angle]==ROIang,1,'first'));
                    objarray(oidx).MIP.Fly = objarray(oidx).Fly;
                end
                
            else
                if isempty(objarray(oidx).MIP.ROI(ROIang))
                    continue
                else
                    inclObjIdx(end+1) = oidx;
                    objarray(oidx).MIP.Fly = objarray(oidx).Fly;
                    
                end
            end
        end
        if ~isempty(inclObjIdx)
                if ~manualROIs
                    [~,~,stimTimes]=plotTrialsExp6([objarray([inclObjIdx]).MIP],1,[],errorbar,[],f);
                else
                    [~,~,stimTimes]=plotTrialsExp6([objarray([inclObjIdx]).MIP],ROIang,[],errorbar,[],f);
                end
                superInclObjIdx = [superInclObjIdx inclObjIdx];
            end
        
    end
        if isempty(superInclObjIdx)
            close(f)
            disp('Wrong experiment, or no matching ROIs')
            return
        end
        % add stim lines in response plot
        getAVPplotParams
        xstim=[mean(stimTimes,1); mean(stimTimes,1)];
        ylims = ylim;
        ymax=ylims(2);
        ymin=ylims(1);
        ystim= [ymax;ymin].*ones(2,size(xstim,2));
        % Lhandle = line(ax(1), xstim,ystim,'color',lightGreyCol,'LineWidth',0.5);
        ylim([ylims(1) ylims(2)]);
        % uistack(Lhandle,'bottom')
        addStimShading(ax(1),xstim(1,:))
        %  addAngleLabels(ax(1),xstim(1,:))
        
        
        if addMotorPositionIndicator
            % add aligned stim trace in 2nd axes underneath
            ax(2) = subplot(numAxes,1,2);
            addMotorPositionPlot(ax(2),xstim)
            addStimShading(ax,xstim(1,:))
            pbaspect(ax,[3 1 1])
            setAVPaxes(ax,2.5)
            addAngleIndicators(ax,xstim(1,:))
        end
        pbaspect(ax(1),[3 1 1])
        if addMotorVelocityIndicator
            addMotorVelocityPlot(ax(1),xstim(1,:));
        end
        linkaxes(ax,'x')
        
%         addAngleIndicators(ax(1),xstim(1,:))
        if ~isequal(numCycles,[1,4])
            ax(1).XLim(1) =  xstim(1,floor(numCycles(1)*size(xstim,2)/4));
            
            ax(1).XLim(2) =  xstim(1,floor(numCycles(2)*size(xstim,2)/4));
        end
        
        setAVPaxes(ax(1),3,5)
        tightfig(f)
    end
    
end


    function addStimShading(ax,xstim)
        clims = [0.85,0.95];
        cols = linspace(clims(2),clims(1),4)'; % grey cols for 4 angles
        cols = [cols;cols;cols;cols];
        ylims = ylim(ax);
        ymax=ylims(2);
        ymin=ylims(1);
        for n = 1:length(xstim)-1
            
            if mod(n,2) % large patch for angle presentation
                c0 = cols(floor(n/2)+1);
                p = patch(ax,[xstim(n) xstim(n+1) xstim(n+1) xstim(n)], [ymin ymin ymax ymax], c0.*ones(1,3));
                
            else % small gradient patch for changing angle
                %         c1 = cols(floor((n+1)/2)+1);
                %         p = patch(ax,[xstim(n) xstim(n+1) xstim(n+1) xstim(n)], [ymin ymin ymax ymax], [c0 c1 c1 c0]);
                
            end
            colormap(ax, gray)
            p.EdgeAlpha = 0;
            ax.CLim = [0 1];
            uistack(p,'bottom')
        end
        ylim(ax,[ylims(1) ylims(2)]);
        
    end

    function addAngleLabels(ax,xstim)
        ylims = ylim(ax);
        ymax=ylims(2);
        ymin=ylims(1);
        angs = wrapTo180(-[30:30:720]-270);
        for n = 1:length(angs)
            
            if abs(angs(n)) == 0 || abs(angs(n)) ==90
                x0 = xstim(2*n-1);
                x1 = xstim(2*n);
                x = mean([x0,x1]);
                text(x,ymax,[num2str(abs(angs(n))) char(176)],'FontSize',7);
                
            end
        end
        ylim(ax,[ylims(1) ylims(2)]);
        ax.Clipping = 'off';
        
        
    end

    function addAngleIndicators(ax,xstim)
        ylims = ylim(ax);
        ymax=ylims(2);
        ymin=ylims(1);
        xlims = xlim(ax);
        xmax=xlims(2);
        xmin=xlims(1);
        angs = [30:30:720];
        col = flipud(hsv(6));
        for n = 1:length(angs)
            
            %  if abs(angs(n)) == 0 || abs(angs(n)) ==90
            % First find the centre of the line:
            % x-position for marker
            x0 = xstim(2*n-1);
            x1 = xstim(2*n);
            x = mean([x0,x1]);
            
            % y-position for maker
            y = ymax*1.08;
            
            % Now find the coordinates of the ends of the line
            
            % line is expressed in data units, so work out x/y scaling to give
            % appearance of equal lengths
            
            xyr = (xmax-xmin)/(ymax-ymin); %ratio x:y units
            pb = pbaspect(ax);
            pbxyr = pb(1)/pb(2); %ratio of x:y aspect
            
            Ay=(ymax-ymin)*0.04; % found after some trial and erro
            Ax = Ay*xyr/(pbxyr); % scale x accordingly
            
            Y = [Ay*sind(angs(n)) Ay*sind(angs(n)+180)];
            X = [Ax*cosd(angs(n)) Ax*cosd(angs(n)+180)];
            
            % add to line center coords
            X = X+x;
            Y = Y+y;
            
            L = line(ax,X,Y,'LineWidth',1);
            L.Color = col(mod(n-1,6)+1,:);
            %  end
            
            ylim(ax,[ylims(1) ylims(2)]);
            
        end
        ax.Clipping = 'off';
        
        
    end


    function addMotorPositionPlot(ax,xstim)
        
        angs = wrapTo180([60:-30:-630]);
        angs(angs==-180)=180;
        yangs = reshape([angs;angs],1,size(xstim,2));
        splitIdx = find(yangs==180);
        
        yAngs1 = [yangs(1:splitIdx(1)-1) -180];
        xStim1 = xstim(1,1:splitIdx(1));
        
        yAngs2 = [yangs(splitIdx(1):splitIdx(3)-1) -180];
        xStim2 = xstim(1,splitIdx(1):splitIdx(3));
        
        yAngs3 = yangs(splitIdx(3):end);
        xStim3 = xstim(1,splitIdx(3):end);
        
        line(ax,xStim1',yAngs1,'color','k');
        line(ax,xStim2',yAngs2,'color','k');
        line(ax,xStim3',yAngs3,'color','k');
        ax.XAxis.Visible = 'off';
        % pbaspect(ax(2),[10 1 1])
        ax.YLim = [-180 180];
        ax.YTick = [-180:90:180];
        sstim= [180;-180].*ones(2,size(xstim,2));
        % Lhandle = line(ax(2), xstim,sstim,'color',darkGreyCol,'LineWidth',0.1);
        % uistack(Lhandle,'bottom')
        
        % ax(2).XAxis.MinorTickValues = xstim(1,:);
        % ax(2).XMinorGrid = 'on';
        % ax(2).GridAlpha = 0;
        xlims = xlim(ax);
        xmax=xlims(2);
        xmin=xlims(1);
        astim = [-180:30:180; -180:30:180];
        tstim= [xmax;xmin].*ones(2,size(astim,2));
        % Lhandle = line(ax(2), tstim,astim,'color',darkGreyCol,'LineWidth',0.1);
        
        % ax(2).YAxis.MinorTickValues = [-180:30:180];
        % ax(2).YMinorGrid = 'on';
        % uistack(Lhandle,'bottom')
        ax.XLim(2) = xstim(1,end);
    end



    function addMotorVelocityPlot(ax,xstim)
        % assumes axes already exist and other data are plotted
        xlims = xlim(ax);
        xmax=xlims(2);
        xmin=xlims(1);
        ylims = ylim(ax);
        ymax=ylims(2);
        ymin=ylims(1);
        
        angs = wrapTo180([60:-30:-630]);
        angs(angs==-180)=180;
        yangs = reshape([angs;angs],1,size(xstim,2));
        splitIdx = find(yangs==180);
        xyscale = 1/min(diff(xstim)); % short time interval where motor moves
        
        % get ang vel from position: Take abs value to avoid dividing by zero later (and we don't care if
        % the actualy plot shows upward or downward deflection). Apply scale to get
        % angles per second
        yAngVel1 = abs(xyscale*diff([yangs(1:splitIdx(1)-1) -180]));
        yVel1 = reshape([yAngVel1;yAngVel1],1,[]);
        x1 = xstim(1,1:splitIdx(1));
        xStim1 = reshape([x1;x1],1,[]);
        xStim1(1) = []; xStim1(end) = [];
        
        yAngVel2 = abs(xyscale*diff([yangs(splitIdx(1):splitIdx(3)-1) -180]));
        yVel2 = reshape([yAngVel2;yAngVel2],1,[]);
        x2 = xstim(1,splitIdx(1):splitIdx(3));
        xStim2 = reshape([x2;x2],1,[]);
        xStim2(1) = []; xStim2(end) = [];
        
        yAngVel3 = abs(xyscale*diff(yangs(splitIdx(3):end)));
        yVel3 = reshape([yAngVel3;yAngVel3],1,[]);
        x3 = xstim(1,splitIdx(3):end);
        xStim3 = reshape([x3;x3],1,[]);
        xStim3(1) = []; xStim3(end) = [];
        
        yscale = (ymax-ymin)/(35*max(yVel3(:)));
        
        % ax.XAxis.Visible = 'off';
        % pbaspect(ax(2),[10 1 1])
        % ax.YLim = [-180 180];
        % ax.YTick = [-180:90:180];
        % sstim= [180;-180].*ones(2,size(xstim,2));
        % Lhandle = line(ax(2), xstim,sstim,'color',darkGreyCol,'LineWidth',0.1);
        % uistack(Lhandle,'bottom')
        
        % ax(2).XAxis.MinorTickValues = xstim(1,:);
        % ax(2).XMinorGrid = 'on';
        % ax(2).GridAlpha = 0;
        
        yposition = 1.01*ymax;
        line(ax,[xStim1'; xStim2'; xStim3'],yposition + yscale*[yVel1 yVel2 yVel3],'color','k','LineWidth',0.2);
        
        
        t0 = text(ax,0,yposition,['0' char(176) '/s']);
        t1 = text(ax,0,yposition+yscale*max(yVel3(:)),['30' char(176) '/s']);
        %  t0 = text(ax,0,yposition,[num2str(floor(min(yVel3(:)))) char(176) '/s ']);
        %  t1 = text(ax,0,yposition+yscale*max(yVel3(:)),[num2str(round(max(yVel3(:)))) char(176) '/s ']);
        % t1.VerticalAlignment = 'bottom';
        % t0.VerticalAlignment = 'top';
        t0.HorizontalAlignment = 'right';
        t1.HorizontalAlignment = 'right';
        t0.FontSize = 2;
        t1.FontSize = 2;
        % astim = [-180:30:180; -180:30:180];
        % tstim= [xmax;xmin].*ones(2,size(astim,2));
        % Lhandle = line(ax(2), tstim,astim,'color',darkGreyCol,'LineWidth',0.1);
        
        % ax(2).YAxis.MinorTickValues = [-180:30:180];
        % ax(2).YMinorGrid = 'on';
        % uistack(Lhandle,'bottom')
        ax.XLim(2) = xstim(1,end);
        ax.YLim(1) = ymin;
        ax.YLim(2) = ymax;
        ax.Clipping = 'off';
    end