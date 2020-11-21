function plotPolData(obj)
if isempty(obj.fullROI)
   getPolROIsOrig(obj); 
end
figure('color','w')
ax1a = subplot(6,8,[1:2,9:10,17:18]);
ax1b = subplot(6,8,[3:4,11:12,19:20]);
ax1c = subplot(6,8,[5:6,13:14,21:22]);
ax1d = subplot(6,8,[7:8,15:16,23:24]);
% ax1 = subplot(6,8,[1:4,9:12,17:20]);
% ax2 = subplot(6,8,[5:8,13:16,21:24]);

ax3 = subplot(6,8,[26:28,34:36]);
pAx = subplot(6,8,[29:31,37:39],polaraxes);
hold(ax3,'on');
% hold(ax4,'on');

expNum = 4;

for oidx = 1:length(obj)
    
 
   plotFullAngMaskImg(obj,ax1a)
   plotPolSelImg(obj,ax1b)
     plotPolMagImg(obj,ax1c)
      plotPolAvgImg(obj,ax1d)

    polAngs = obj.expPolAngs(expNum);
plotPolarKey(obj,pAx)
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
%             plot(ax3,polAngs(aIdx),fullPolResp(rIdx,aIdx),'o','color',obj.ROI(rIdx).color,'MarkerFaceColor',obj.ROI(rIdx).color)
        end
%         plot(ax3,polAngs,fullPolResp(rIdx,:),'color',obj.ROI(rIdx).color,'LineWidth',2)
normresp = fullPolResp(rIdx,:)./max(fullPolResp(rIdx,:));
normresp = normresp-mean(normresp);
  plot(ax3,polAngs,normresp,'color',obj.ROI(rIdx).color,'LineWidth',2)

%   plot(ax3,polAngs,normresp,'.','color',obj.ROI(rIdx).color,'LineWidth',50)
  scatter(ax3,polAngs,normresp,20,obj.ROI(rIdx).color)

%         plot(ax4,polAngs,exPolResp(rIdx,:),'color',obj.ROI(rIdx).color,'LineWidth',2)
        
     
   pAxF(rIdx) = subplot(6,8,40+rIdx,polaraxes);
    makePolarROIplot(obj,pAxF(rIdx),rIdx,polAngs,fullPolResp)


    end
   
    for pAxIdx = 1:length(pAxF)
    pAxF(pAxIdx).RLim(2) = max([pAxF.RLim]);
    end
    
    ax3.XTick = [ 90 180 270 360];
%     ax4.XTick = [ 90 180 270 360];                
    ax3.YLabel.String = 'dF/F';
%     ax4.YLabel.String = 'dF/F';
 ax3.XLabel.String = 'angle (\circ)';
%  ax4.XLabel.String = 'angle (\circ)';
 ax3.XLim(1) = 30;
%   ax4.XLim(1) = 30;

%  axis(ax3,'tight')
%  axis(ax4,'tight')
    offsetAxesXTickONLY(ax3)
%     offsetAxesXTickONLY(ax4)

recordStr = erase(obj.Folder,'Y:\ben\avp_pol\18_data\');
bslashIdx = strfind(recordStr,'\');
bentitle(recordStr(1:bslashIdx(4)-1));
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

lastFrame = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first')) - 1 - (floor(1*obj.IFI/obj.AIrate));
firstFrame = lastFrame - (floor(5*obj.IFI/obj.AIrate));
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
pAx.RAxis.Visible = 'off';
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

% % add +ve responses to polar axes
% rhoEx = exPolResp(rIdx,:);
% rhoEx(rhoEx<0) = 0;
% polarplot(axEx,[theta theta(1)],[rhoEx rhoEx(1)],'-','color',obj.ROI(rIdx).color,'LineWidth',2);
% % add +ve responses to polar axes
% rhoEx = -exPolResp(rIdx,:);
% rhoEx(rhoEx<0) = 0;
% polarplot(axEx,[theta theta(1)],[rhoEx rhoEx(1)],'-','color', obj.ROI(rIdx).color-(obj.ROI(rIdx).color)*0.25,'LineWidth',2);
pAx.RAxis.Label.String = '    dF/F';
pAx.RAxis.Label.Visible = 'on';
pAx.RAxis.Label.Rotation = 0;
pAx.RAxis.Label.HorizontalAlignment = 'left';
% pAx.RAxis.Limits = [0 1];
pAx.RAxis.TickValues  = [0:0.25:1];
pAx.RAxis.TickLabels = {'';'';'0.5';'';'1'};
pAx.RAxisLocation = 45;
pAx.FontSize = 13;

end

function plotPolSelImg(obj,ax)
imagesc(obj.polSelImg,'Parent',ax) % threshold at 10 is just to get rid of edges
colormap(ax,flipud(hot))

title(ax,['selectivity (Avg bins)[thresh>' num2str(obj.polSelThreshold) ']']);
axis(ax,'image');
    ax.XTick = []; ax.XTickLabel = {};
        ax.YTick = []; ax.YTickLabel = {};

ax.Color = [0 0 0];
colorbar(ax);
ax.CLim(1) = obj.polSelThreshold;
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
    clims = ax.CLim;
   
%     ax.CLim = [(100-obj.polMagThreshold)*clims(2)/100 clims(2)];
end
function plotPolAvgImg(obj,ax)    
    imagesc(obj.avgPolImg,'Parent',ax)
    title(ax,'median activity');
       colormap(ax,flipud(hot))
    axis(ax,'image');
    ax.XTick = []; ax.XTickLabel = {};
        ax.YTick = []; ax.YTickLabel = {};

    colorbar(ax);
end
function plotPolarKey(obj,pAx)
    hold(pAx,'on');
    set(pAx,'ThetaDir','clockwise','ThetaZeroLocation','right')
    pAx.RAxis.Visible = 'off';
    pAx.RLim = [0 1];
for rIdx = 1:length(obj.fullROI)
            % add line to polar axes
            theta = deg2rad(obj.fullROI(rIdx).angle);
            if isnan(theta)
                theta = 0;
            end
            rho = mean(obj.polSelImg(obj.fullROI(rIdx).mask));
            polarplot(pAx,[theta-pi theta],[rho rho],'-','color',obj.fullROI(rIdx).color,'LineWidth',3)
            
end
end
function plotFullAngMaskImg(obj,ax)
image(obj.fullAngMaskImg,'Parent',ax)
axis(ax,'image')
ax.XTick = [];
ax.YTick = [];
end
