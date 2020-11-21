function varargout = plotFlashResp(objarray,roiAngs,showErrorBar)
if nargin < 2 || isempty(roiAngs)
    roiAngs = [30 90 120]; % angle as recorded
        manualROIs = 0;
elseif all(roiAngs<15) && all(roiAngs>0)
    manualROIs = 1;
else
    manualROIs = 0;
end
if nargin < 3 || isempty(showErrorBar)
    showErrorBar = 1;
end
% f = figure; % changed to make new figure for each flash resp

padStart = 2;
getAVPplotParams

plotangs = [90 180];
if ~manualROIs
superUseMSP(objarray,1)
end
% for  oidx = 1:length(objarray)
%     if objarray(oidx).containsPolMapExp > 2
%    getPolMaps(objarray(oidx).MIP)
%    getPolROIs(objarray(oidx).MIP)
%     end
% end
try
    for angIdx = 1:length(plotangs)
        f(angIdx) = figure;
        ax(angIdx) = subplot(1,length(plotangs),angIdx);
        hold on
        ang = plotangs(angIdx);
        for roiIdx = 1:length(roiAngs)
            inclObjIdx = [];
            
            ROIang = roiAngs(roiIdx);
            
            for oidx = 1:length(objarray)
                
                
                % First check each object contains pol flash exp (3) and a pol map exp
                % (4,8,10)
                if ~any(ismember(3,objarray(oidx).Exps)) ||(objarray(oidx).containsPolMapExp < 4)
                    continue
                end
                if isempty(objarray(oidx).MIP.polSelImg)
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
                    
                else
                    if isempty( objarray(oidx).MIP.ROI )
                        loadROIs(objarray(oidx).MIP)
                    end
                end
                
                 if ~manualROIs && isempty(find([objarray(oidx).MIP.polROI.angle]==ROIang,1,'first'))
                    continue
                else
                    inclObjIdx(end+1) = oidx;
                    
                    % we will temporarily discard all but the exp3 trials. Store the original
                    % vectors for trial start stop times etc so we can restore them later
                    keepExpTrials(objarray(oidx).MIP,3)
                    
                    
                    % Pad each trial by +4s before and +4s after trial
                    addToTrialStart(objarray(oidx).MIP,padStart)
                    addToTrialEnd(objarray(oidx).MIP,padStart+2)
                    if ~manualROIs
                        % assign this polROI to ROI(1)
                        objarray(oidx).MIP.ROI = objarray(oidx).MIP.polROI(find([objarray(oidx).MIP.polROI.angle]==ROIang,1,'first'));
                    end
                    objarray(oidx).MIP.Fly = objarray(oidx).Fly;
                end
            end
            if ~isempty(inclObjIdx)
                % Plot
                fields = [];
                fields.TrialPatNum = ang;
                fields.TrialSeqNum = 3;
                if ~manualROIs
                    outputData(roiIdx,angIdx) = plotTrialsExp3([objarray([inclObjIdx]).MIP],1,fields,showErrorBar,[],f(angIdx));
                else
                    outputData(roiIdx,angIdx) = plotTrialsExp3([objarray([inclObjIdx]).MIP],ROIang,fields,showErrorBar,[],f(angIdx));
                end
                % Restore original data when we're finished:
                for oidx = inclObjIdx
                    restoreAllTrials(objarray(oidx).MIP)
                end
            end
        end
        pbaspect(ax(angIdx),[1 1.5 1])
        
%         ax(angIdx).XAxis.MinorTickValues = [-99.5:1:99.5];
%         ax(angIdx).XMinorGrid = 'on';
%         ax(angIdx).YAxis.MinorTickValues = [-990:0.25:99];
%         ax(angIdx).YMinorGrid = 'on';
%         ax(angIdx).YGrid = 'off';
%         ax(angIdx).XGrid = 'off';
%         ax(angIdx).Clipping = 'off'; 


    end
catch ME
    % Restore original data when we're finished:
    for oidx = 1:length(objarray)
        restoreAllTrials(objarray(oidx).MIP)
    end
    rethrow(ME)
end
ylims = [min(cellfun(@min,ylim(ax))) max(cellfun(@max,ylim(ax)))];
% linkaxes( ax)
for aidx = 1:length(ax)
    ax(aidx).YLim = ylims;
    % 3 part patch to show start/stop times may vary by up to 1s:
%     patchH(aidx) = plotStimPatch(objarray(1).MIP,ax(aidx),[padStart-1 padStart],[], [.95 .95 .95]);
%     patchH(aidx) = plotStimPatch(objarray(1).MIP,ax(aidx),[padStart padStart+3],[], lightGreyCol);
%     patchH(aidx) = plotStimPatch(objarray(1).MIP,ax(aidx),[padStart+3 padStart+4],[], [.95 .95 .95]);
patchH(aidx) = plotStimPatch(objarray(1).MIP,ax(aidx),[padStart-0.5 padStart+3.5],[], lightGreyCol);

outputData(1,aidx).patch(aidx) = patchH(aidx);
ax(aidx).XTick = [padStart-0.5:1:padStart+3.5];
ax(aidx).XTickLabel = {'0';'';'';'';'4 s'};
ax(aidx).XLabel = [];   
ax(aidx).ActivePositionProperty = 'position';

if ~mod( plotangs(aidx), 180) % show 90 instead of -90
title(ax(aidx),[num2str(abs(wrapTo180(-plotangs(aidx)-270))) char(176)])   
else
title(ax(aidx),[num2str(wrapTo180(-plotangs(aidx)-270)) char(176)] )        
end
offsetAxes(ax(aidx))
% scalebarF(ax(aidx))
% scalebarT(ax(aidx))
setAVPaxes(ax(aidx),defaultAxisHeight_cm)
tightfig(f(aidx))

end

if nargout==1
    varargout{1} = outputData;
end
end