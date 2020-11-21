function varargout = plotFlashRespManual(objarray,ROIidx,angs,plotCtrl)
assert(nargin>1 && ~isempty(ROIidx),'specify an ROI index')
if nargin<4 || isempty(plotCtrl)
    flashExp=3;
else
    assert((plotCtrl==1|plotCtrl==0),'choose to plot ctrls only: 0 or 1')
    if plotCtrl == 1
        flashExp=1;
    else
        flashExp=3;
    end
end
if nargin >2
    plotangs = angs;
else
    if flashExp == 1
        plotangs = 0;
    else
        plotangs = [90 180];
    end
end
% f = figure; % changed to make new figure for each flash resp

padStart = 2;
getAVPplotParams

superUseMSP(objarray,1)


try
    for angIdx = 1:length(plotangs)
        f(angIdx) = figure;
        ax(angIdx) = subplot(1,length(plotangs),angIdx);
        hold on
        ang = plotangs(angIdx);
        for roiIdx = 1:length(ROIidx)
            inclObjIdx = [];
            
            
            for oidx = 1:length(objarray)
                
                
                % First check each object contains pol flash exp (1,3)
                if ~any(ismember(flashExp,objarray(oidx).Exps))
                    continue
                end
                if flashExp ==3 && ang == 0
                    ang = 180;
                end
                
                
                
                if isempty(objarray(oidx).MIP.ROI) % no rois made
                    continue
                end
                
                inclObjIdx(end+1) = oidx;
                
                % we will temporarily discard all but the exp3 trials. Store the original
                % vectors for trial start stop times etc so we can restore them later
                
                keepExpTrials(objarray(oidx).MIP,flashExp)
                
                
                % Pad each trial by +4s before and +4s after trial
                addToTrialStart(objarray(oidx).MIP,padStart)
                addToTrialEnd(objarray(oidx).MIP,padStart+2)
                
                objarray(oidx).MIP.Fly = objarray(oidx).Fly;
                
            end
            try
                if ~isempty(inclObjIdx)
                    % Plot
                    fields = [];
                    fields.TrialPatNum = ang;
                    fields.TrialSeqNum = flashExp;
                    
                    outputData(roiIdx,angIdx) = plotTrialsExp3([objarray([inclObjIdx]).MIP],ROIidx(roiIdx),fields,1,[],f(angIdx));
                    
                    % Restore original data when we're finished:
                    for oidx = inclObjIdx
                        restoreAllTrials(objarray(oidx).MIP)
                    end
                end
            catch
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
if size(ax,2)>1
    ylims = [min(cellfun(@min,ylim(ax))) max(cellfun(@max,ylim(ax)))];
else
    ylims = ylim(ax);
end
% linkaxes( ax)
for aidx = 1:length(ax)
    ax(aidx).YLim = ylims;
    % 3 part patch to show start/stop times may vary by up to 1s:
    %     patchH(aidx) = plotStimPatch(objarray(1).MIP,ax(aidx),[padStart-1 padStart],[], [.95 .95 .95]);
    %     patchH(aidx) = plotStimPatch(objarray(1).MIP,ax(aidx),[padStart padStart+3],[], lightGreyCol);
    %     patchH(aidx) = plotStimPatch(objarray(1).MIP,ax(aidx),[padStart+3 padStart+4],[], [.95 .95 .95]);
    patchH(aidx) = plotStimPatch(objarray(1).MIP,ax(aidx),[padStart-0.5 padStart+3.5],[], lightGreyCol);
    outputData(1,aidx).patch(aidx) = patchH(aidx);
    %  addlistener(ax(aidx), 'MarkedClean', @(object,event)setPatchY(ax(aidx),patchH(aidx)));
    ax(aidx).XTick = [padStart-0.5:1:padStart+3.5];
    ax(aidx).XTickLabel = {'0';'';'';'';'4 s'};
    ax(aidx).XLabel = [];
    ax(aidx).ActivePositionProperty = 'position';
    if flashExp~=1 && ~mod( plotangs(aidx), 180) % show 90 instead of -90
        title(ax(aidx),[num2str(abs(wrapTo180(-plotangs(aidx)-270))) char(176)])
    elseif flashExp==1
        title(ax(aidx),'unpolarized')
    else
        title(ax(aidx),[num2str(wrapTo180(-plotangs(aidx)-270)) char(176)] )
    end
    offsetAxes(ax(aidx))
%     scalebarF(ax(aidx))
%     scalebarT(ax(aidx))
    setAVPaxes(ax(aidx),defaultAxisHeight_cm)
    tightfig(f(aidx))
    

end


if nargout==1
    varargout{1} = outputData;
end

end

function setPatchY(ax,patch)

upper  = ax.YLim(2);
lower = ax.YLim(1);
patch.YData(1) = lower;
patch.YData(4) = lower;
patch.YData(2) = upper;
patch.YData(3) = upper;
ax.YLimMode = 'auto';
end