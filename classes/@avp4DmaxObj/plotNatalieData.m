function plotNatalieData(obj)
if isempty(obj.fullAngMaskImg)
   getPolROIs(obj); 
end
pathname = 'C:\Users\sapiens\Desktop\DmDRA split data\';
% ax1a = subplot(6,8,[1:2,9:10,17:18]);
% ax1b = subplot(6,8,[3:4,11:12,19:20]);
% ax1c = subplot(6,8,[5:6,13:14,21:22]);
% ax1d = subplot(6,8,[7:8,15:16,23:24]);
% % ax1 = subplot(6,8,[1:4,9:12,17:20]);
% % ax2 = subplot(6,8,[5:8,13:16,21:24]);
% 
% ax3 = subplot(6,8,[26:28,34:36]);
% hold(ax3,'on');
% % hold(ax4,'on');
% 
expNum = 4;

for oidx = 1:length(obj)
    
 figure('color','w')
   plotFullAngMaskImg(obj,gca)
          set(gcf, 'InvertHardCopy', 'off');
    print(gcf, fullfile(pathname,'pol_tuning_map') , '-dpng','-r450');
imgPos = get(gcf,'position');
axPos = get(gca,'position');
   figure('color','w')
   plotPolSelImg(obj,gca)
             set(gcf, 'InvertHardCopy', 'off');
    print(gcf, fullfile(pathname,'pol_selectivity_map') , '-dpng','-r450');

     figure('color','w')
      plotPolAvgImg(obj,gca)
     set(gcf,'position',imgPos);
          set(gca,'position',axPos);

          set(gcf, 'InvertHardCopy', 'off');
    print(gcf, fullfile(pathname,'activity_image') , '-dpng','-r450');

      figure('color','w','position',[500 500 200 200])
pAx = polaraxes;

    polAngs = obj.expPolAngs(expNum);
plotPolarKey(obj,pAx)
          set(gcf, 'InvertHardCopy', 'off');
    print(gcf, fullfile(pathname,'pol_tuning_legend') , '-dpng','-r450');

      pF = figure('color','w','position',[50 500 1800 160]);
      
      p3 = figure('color','w');
ax3 = gca;
hold on
    exPolResp = nan(length(obj.ROI),0.5*length(polAngs));
    fullPolResp = nan(length(obj.ROI),0.5*length(polAngs));
    
    % For each ROI:
    for rIdx = 1:length(obj.ROI)            
        
        % Find F0 for this ROI
        F0FrameArray = findF0LogicArray(obj,expNum);
        exF0(rIdx) = mean(obj.ROI(rIdx).response(F0FrameArray));
        fullF0(rIdx) = mean(obj.fullROI(rIdx).response(F0FrameArray));
        
        % For each angle:
        for aIdx = 1:length(polAngs)
            
            % Get logic array that indicates every frame that was recorded at the
            % angle (NOT including angle+180)
            angleFrameArray = getAngleLogicArrayNot180(obj,polAngs(aIdx));
                       
            
            % Extract mean response for this angle:
%             exPolResp(rIdx,aIdx) = (mean(obj.ROI(rIdx).response(angleFrameArray))/exF0(rIdx)) - 1;
%             plot(ax4,polAngs(aIdx),exPolResp(rIdx,aIdx),'o','color',obj.ROI(rIdx).color,'MarkerFaceColor',obj.ROI(rIdx).color)
            
            fullPolResp(rIdx,aIdx) = (mean(obj.fullROI(rIdx).response(angleFrameArray))/fullF0(rIdx)) - 1;
            plot(ax3,polAngs(aIdx),fullPolResp(rIdx,aIdx),'o','color',obj.ROI(rIdx).color,'MarkerFaceColor',obj.ROI(rIdx).color)
        end
        plot(ax3,polAngs,fullPolResp(rIdx,:),'color',obj.ROI(rIdx).color,'LineWidth',2)
%         plot(ax4,polAngs,exPolResp(rIdx,:),'color',obj.ROI(rIdx).color,'LineWidth',2)
                  
    
% figure('color','w')

   pAxF(rIdx) = subplot(1,7,rIdx,polaraxes,'parent',pF);
    makePolarROIplot(obj,pAxF(rIdx),rIdx,polAngs,fullPolResp)


    end
   
    for pAxIdx = 1:length(pAxF)
    pAxF(pAxIdx).RLim(2) = max([pAxF.RLim]);
    end
    
              set(gcf, 'InvertHardCopy', 'off');
    print(pF, fullfile(pathname,'activity_polar_plots') , '-dpng','-r450');
        print(pF, fullfile(pathname,'activity_polar_plots') , '-svg');


    ax3.XTick = [30:30:360];   

    ax3.XTickLabels = {'';'';'90';'';'';'180';'';'';'270';'';'';'360';};
%     ax4.XTick = [ 90 180 270 360];                
    ax3.YLabel.String = '\DeltaF/F';
%     ax4.YLabel.String = 'dF/F';
 ax3.XLabel.String = 'angle (\circ)';
%  ax4.XLabel.String = 'angle (\circ)';
 ax3.XLim(1) = 0;
 
%   ax4.XLim(1) = 30;
ax3.YLim = [-1 1];
ax3.YTick = [-1:0.25:1];
%  axis(ax3,'tight')
ax3.YTickLabels = {'-1';'';'-0.5';'';'0';'';'0.5';'';'1'};
%  axis(ax4,'tight')
%     offsetAxesXTickONLY(ax4)
    offsetAxesXTickONLY(ax3)
title('GCaMP activity')
ax3.FontSize = 13;
ax3.Clipping = 'off';

recordStr = erase(obj.Folder,'Y:\ben\avp_pol\18_data\');
bslashIdx = strfind(recordStr,'\');
% bentitle(recordStr(1:bslashIdx(4)-1));
          set(gcf, 'InvertHardCopy', 'off');
    print(gcf, fullfile(pathname,'activity_time_series') , '-dpng','-r450');

end
end



function angleFrameArray = getAngleLogicArrayNot180(obj,polAng)

% Extract the start:stop frame numbers for the pol tuning experiment:

if any(obj.TrialSeqNum==4)
    expNum = 4;
elseif any(obj.TrialSeqNum==2)
    expNum = 2;
else
    return
    error('No pol tuning exp (2 or 4) available')
end

polStart = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first'));
polStop = obj.TrialEndFrame(find(obj.TrialSeqNum == expNum,1,'last'));


% Find the frames recorded with polarizer at polAng:
angStarts = obj.TrialStartFrame(obj.TrialPatNum == polAng);
angStops = obj.TrialEndFrame(obj.TrialPatNum == polAng);

% Add the equivalent +180deg angle
% angStarts = [angStarts, obj.TrialStartFrame(obj.TrialPatNum == polAng + 180)];
% angStops = [angStops, obj.TrialEndFrame(obj.TrialPatNum == polAng + 180)];

% Construct logical array for frameAngles == polAng
angleFrameArray = zeros(1,size(obj.Frames,3));
for n = 1:length(angStarts)
    if (angStarts(n) >= polStart) && (angStops(n) <= polStop)
        angleFrameArray(angStarts(n):angStops(n)) = 1;
    end
end

angleFrameArray = logical(angleFrameArray);
end


function angleFrameArray = getAngleLogicArray(obj,polAng)

% Extract the start:stop frame numbers for the pol tuning experiment:

if any(obj.TrialSeqNum==4)
    expNum = 4;
elseif any(obj.TrialSeqNum==2)
    expNum = 2;
else
    return
    error('No pol tuning exp (2 or 4) available')
end

polStart = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first'));
polStop = obj.TrialEndFrame(find(obj.TrialSeqNum == expNum,1,'last'));


% Find the frames recorded with polarizer at polAng:
angStarts = obj.TrialStartFrame(obj.TrialPatNum == polAng);
angStops = obj.TrialEndFrame(obj.TrialPatNum == polAng);

% Add the equivalent +180deg angle
angStarts = [angStarts, obj.TrialStartFrame(obj.TrialPatNum == polAng + 180)];
angStops = [angStops, obj.TrialEndFrame(obj.TrialPatNum == polAng + 180)];

% Construct logical array for frameAngles == polAng
angleFrameArray = zeros(1,size(obj.Frames,3));
for n = 1:length(angStarts)
    if (angStarts(n) >= polStart) && (angStops(n) <= polStop)
        angleFrameArray(angStarts(n):angStops(n)) = 1;
    end
end

angleFrameArray = logical(angleFrameArray);
end

function F0FrameArray = findF0LogicArray(obj,expNum)
% 5sec recorded after each exp, 5sec recorded before each. Use 5s before. 
% We don't have a marker for where the 'set' starts, only the first trial.
% However, before each trial there is at least 1sec where the LED is on
% before the trial start marker, which we can use

lastFrame = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first')) - 1;
firstFrame = lastFrame - (floor(10*obj.IFI/obj.AIrate));
if firstFrame < 1
    firstFrame = 1;
end

% lastFrame = obj.TrialStartFrame(1);
% 
% firstFrame = lastFrame - (ceil(1*obj.IFI/obj.AIrate));%- (floor(1*obj.IFI/obj.AIrate));
% if firstFrame < 1
%     firstFrame = 1;
% end
% Construct logical array for F0 == 1
F0FrameArray = zeros(1,size(obj.Frames,3));
F0FrameArray(firstFrame:lastFrame)=1;

F0FrameArray = logical(F0FrameArray);
end

function makePolarROIplot(obj,pAx,rIdx,polAngs,fullPolResp,exPolResp)
% if rIdx > 4
%    axNum = axNum + 4; 
% end
%       axEx = subplot(6,8,axNum+4,polaraxes);

set(pAx,'ThetaDir','clockwise','ThetaZeroLocation','right')
% pAx.RAxis.Visible = 'off';
% set(axEx,'ThetaDir','clockwise','ThetaZeroLocation','right')
% axEx.RAxis.Visible = 'off';
hold(pAx,'on');
% hold(axEx,'on');
    pAx.ThetaTick = [0:30:360];
%     axEx.ThetaTick = [0:30:360];                
    pAx.ThetaTickLabel = {'0\circ';'';'';'90\circ';'';'';'180\circ';'';'';'270\circ';'';''};
%     axEx.ThetaTickLabel = {'0';'';'';'90';'';'';'180';'';'';'270';'';''};

theta = [deg2rad(polAngs)];

% add +ve responses to polar axes
rhoFull = fullPolResp(rIdx,:);
rhoFull(rhoFull<0) = 0;
polarplot(pAx,[theta theta(1)],[rhoFull rhoFull(1)],'-','color',obj.ROI(rIdx).color,'LineWidth',2);
% add -ve responses to polar axes
rhoFull = -fullPolResp(rIdx,:);
rhoFull(rhoFull<0) = 0;
polarplot(pAx,[theta theta(1)],[rhoFull rhoFull(1)],':','color', obj.ROI(rIdx).color-(obj.ROI(rIdx).color)*0.5,'LineWidth',2);
pAx.RAxis.Label.String = '    \DeltaF/F';
pAx.RAxis.Label.Visible = 'on';
pAx.RAxis.Label.Rotation = 0;
pAx.RAxis.Label.HorizontalAlignment = 'left';
pAx.RAxis.Limits = [0 1];
pAx.RAxis.TickValues  = [0:0.25:1];
pAx.RAxis.TickLabels = {'';'';'0.5';'';'1'};
pAx.RAxisLocation = 45;
pAx.FontSize = 13;
% % add +ve responses to polar axes
% rhoEx = exPolResp(rIdx,:);
% rhoEx(rhoEx<0) = 0;
% polarplot(axEx,[theta theta(1)],[rhoEx rhoEx(1)],'-','color',obj.ROI(rIdx).color,'LineWidth',2);
% % add +ve responses to polar axes
% rhoEx = -exPolResp(rIdx,:);
% rhoEx(rhoEx<0) = 0;
% polarplot(axEx,[theta theta(1)],[rhoEx rhoEx(1)],'-','color', obj.ROI(rIdx).color-(obj.ROI(rIdx).color)*0.25,'LineWidth',2);

end

function plotPolSelImg(obj,ax)
imagesc(imrotate(obj.polSelImg,90),'Parent',ax) % threshold at 10 is just to get rid of edges
colormap(ax,flipud(hot))

title(ax,['polarization selectivity']);
axis(ax,'image');
    ax.XTick = []; ax.XTickLabel = {};
        ax.YTick = []; ax.YTickLabel = {};

ax.Color = [0 0 0];
    c= colorbar(ax);
    c.Ticks = [0:0.25:1];
    c.TickLabels = {'0';'';'';'';'1'};
    c.TickDirection = 'out';

% ax.CLim(1) = obj.polSelThreshold;
ax.CLim = [0 1];

ax.FontSize = 13;

end
function plotPolMagImg(obj,ax)
    imagesc(obj.fftMagImg,'Parent',ax)
             colormap(ax,flipud(hot))

    title(ax,['magnitude (FFT)[take top' num2str(obj.polMagThreshold) '%]']);
    axis(ax,'image');
        ax.XTick = []; ax.XTickLabel = {};
        ax.YTick = []; ax.YTickLabel = {};

    ax.Color = [0 0 0];
    colorbar(ax);
ax.FontSize = 13;

end
function plotPolAvgImg(obj,ax)    
    imshow(imrotate(6.*obj.avgPolImg./max(obj.avgPolImg(:)),90),'Parent',ax)
    title(ax,'average GCaMP activity');
    gmap = zeros(255,3);
    gmap(:,2) = linspace(0,1,255);
    
       colormap(ax,gmap);
    axis(ax,'image');
    ax.XTick = []; ax.XTickLabel = {};
        ax.YTick = []; ax.YTickLabel = {};

%     c= colorbar(ax);
%     c.Ticks = [0 1];
%     c.TickLabels = {'min';'max'};
%         c.TickDirection = 'out';

    ax.FontSize = 13;

end
function plotPolarKey(obj,pAx)
    hold(pAx,'on');
    set(pAx,'ThetaDir','clockwise','ThetaZeroLocation','right')
    pAx.RLim = [0 1];    
    pAx.RAxis.Visible = 'off';
   pAx.RAxis.TickValues = [];
   pAx.ThetaAxis.TickLabels = {'180\circ';'30\circ';'60\circ';'90\circ';'120\circ';'150\circ'};
for rIdx = 1:length(obj.fullROI)-1
            % add line to polar axes
            theta = deg2rad(obj.fullROI(rIdx).angle);
            if isnan(theta)
                theta = 0;
            end
            rho = 1;
            polarplot(pAx,[theta theta],[0.1 rho],'-','color',obj.fullROI(rIdx).color,'LineWidth',5)
            polarplot(pAx,[theta-pi theta-pi],[0.1 rho],'-','color',obj.fullROI(rIdx).color,'LineWidth',5)

            
end
pAx.FontSize = 11;
pAx.FontWeight = 'bold';

end
function plotFullAngMaskImg(obj,ax)
imagesc(imrotate(obj.fullAngMaskImg,90),'Parent',ax)
axis(ax,'image')
ax.XTick = [];
ax.YTick = [];
colormap(flipud(hsv(6)))
    c= colorbar(ax);
    c.Ticks = [1/12:1/6:1-1/12];
    c.TickDirection = 'out';
    c.TickLabels = {'30\circ';'60\circ';'90\circ';'120\circ';'150\circ';'180\circ'};
    c.FontWeight = 'bold';
    title('polarization tuning map')
    ax.FontSize = 13;
end
