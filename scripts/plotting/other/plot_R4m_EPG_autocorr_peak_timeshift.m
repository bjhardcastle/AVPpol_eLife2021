% Plot first peak of autocorrelation function of each ROI's response during
% pol mapping exp in terms of the timeshift of the peak relative to the
% period of the stimulus.
% For TuBu_a (currently not shown) and R4m in anterior bulb, and E-PGs in PB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % load stored mat file if it exists, or try to refresh

pathsAVP
if exist(autocorrStruct_path,'file')
    load(autocorrStruct_path)
else
    try
        refresh_TuBu_a_R4m_EPG_autocorrStruct
    catch
        warning(['Expected to find ' autocorrStruct_path])
        disp('If raw data are available, run ''refresh_TuBu_a_R4m_EPG_autocorrStruct.m''')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % make boxplot

figure,hold on
grps = [autocorrStruct.GrpIdx(:)];
pks = [autocorrStruct.Peakshift(:)];
% deal with nans ( shouldn't be any )
grps(isnan(pks)) = [];
pks(isnan(pks)) = [];
posMod = [3,4,6,7,9,10]; % x-axis positions
grpLength = length(posMod);
alphaVal = 0.05/grpLength;
[bkgHandles,bkgStats] = notBoxPlot(pks,grps,'style','sdline','markMedian',true,'jitter',0.4,'manualCI',alphaVal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % adjust aesthetics

getAVPplotParams
for bIdx = 1:length(bkgHandles)
    
    if ~isempty(bkgHandles(bIdx).data)
        bkgHandles(bIdx).semPtch.FaceColor = 'none';
        bkgHandles(bIdx).semPtch.EdgeColor = 'none';
        bkgHandles(bIdx).semPtch.FaceAlpha = 0.8;
        bkgHandles(bIdx).semPtch.EdgeAlpha = 0.8;
        bkgHandles(bIdx).semPtch.Visible = 'on';
        
        bkgHandles(bIdx).data.Marker = '.';
        bkgHandles(bIdx).data.Color = darkGreyCol;
        if ~mod(bIdx,2)
            bkgHandles(bIdx).data.Color = 'k';
        end
        bkgHandles(bIdx).data.MarkerSize = defaultMarkerSize;
        bkgHandles(bIdx).data.Visible = 'on';
        
        bkgHandles(bIdx).mu.Color = darkGreyCol;
        bkgHandles(bIdx).med.Color = 'k';
        bkgHandles(bIdx).med.Visible = 'on';
        bkgHandles(bIdx).sd.Visible = 'off';
        bkgHandles(bIdx).sd.Color = darkGreyCol;
        bkgHandles(bIdx).sd.LineWidth = 0.5;
        bkgHandles(bIdx).mu.LineWidth = 0.5;
        bkgHandles(bIdx).med.LineWidth = 0.5;
        bkgHandles(bIdx).med.LineStyle = ':';
        
        % set sd to top layer
        bkgHandles(bIdx).sd.ZData  = abs(bkgHandles(bIdx).sd.ZData);
    end
end
uistack([bkgHandles(:).sd],'top')
uistack([bkgHandles(:).mu],'top')
uistack([bkgHandles(:).med],'top')
uistack([bkgHandles(:).data],'top')

ax = gca;

% add limit lines at +/- 2 sec (half the static angle presentation time)
hold(ax,'on')
plot([0,12],[2,2],'r','LineWidth',0.5)
plot([0,12],[-2,-2],'r','LineWidth',0.5)
% add limit lines at +/- 4 sec (static angle presentation time)
% plot([0,12],[4,4],'r:','LineWidth',0.5)
% plot([0,12],[-4,-4],'r:','LineWidth',0.5)



xlim([min(posMod)-0.4, max(posMod)+0.5])
ylim([-5, mean([bkgStats(:).mu])+3*std([bkgStats(:).mu]) ] )


ax.YLabel.String = 'peak shift (-s)';
ax.XLabel.String = '';
% ax.XTick = posMod;
ax.XTick = [3.5 6.5 9.5];
ax.XTickLabel = autocorrStruct.lineCell;
ax.TickLabelInterpreter = 'none';
ax.XTickLabelRotation = 45;

offsetAxesXTickONLY(ax)
setAVPaxes(ax)
ax.YGrid = 'on';
pbaspect([2.5,1,1])
uistack(gca,'bottom');
ax.Layer = 'bottom';

set(gcf,'color','w')
tightfig(gcf)
addExportFigToolbar(gcf)

groupIdx = unique(autocorrStruct.GrpIdx_AllExamined);
if length(bkgStats)<length(groupIdx)
    warning('Some data not plot: likely that no R34H10 ctrl data passed peak-detection threshold.')
end
boxIdx = 0;
cStr = {'ctrl';'exp'};
for sidx = 1:3 % each line
    for cidx = 1:2 % exp2 (ctrl) and exp4
        
        xpos = groupIdx((sidx-1)*2+cidx); % x-axis position
        inclValIdx = find(autocorrStruct.GrpIdx(:)==xpos);
        
        if ~isempty(inclValIdx)
            boxIdx = boxIdx + 1;
            
            fprintf(['[%s %s]\n median=%1.2fs, SD=%1.2f, p=%1.5f signrank\n stim-locked=%1.2f%%, ',...
                '\n n=%d ROIs, N=%d flies included (N=%d examined)\n\n'] ,...
                autocorrStruct.lineCell{sidx} ,...
                cStr{cidx} ,...
                bkgStats(boxIdx).median , ...
                bkgStats(boxIdx).sd , ...
                signrank(autocorrStruct.Peakshift(inclValIdx)) ,...
                100*sum(autocorrStruct.StimLocked(inclValIdx))./length(inclValIdx) ,...
                length(inclValIdx) ,...
                length(unique(autocorrStruct.FlyNumPeriodic(inclValIdx))) ,...
                length(unique(autocorrStruct.FlyNum_AllExamined(autocorrStruct.GrpIdx_AllExamined==xpos))) ...
                );
        end
    end
end