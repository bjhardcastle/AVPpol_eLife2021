% Generate plots of blue flash responses for TuBu drivers in AOTU & BU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % load stored mat file if it exists, or try to refresh
pathsAVP
if exist(TuBu_flash_data_path,'file')
    load(TuBu_flash_data_path,'flash3','flash9')
else
    try
        refresh_TuBu_flash_resp
    catch
        warning(['Expected to find ' TuBu_flash_data_path])
        disp('If raw data are available, run ''refresh_TuBu_flash_resp.m''')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % prepare for plotting

getAVPplotParams
printpath = fig4s1path;
prefix = 'blue_flash';

maskStr = {'cell'}; %originally included 'layer'
areaStr = {'AOTU';'Bu'};
lineStr = { ...
    'R34H10_AOTU';...
    'R88A06_AOTU';...
    'R49E09_AOTU';...
    };

lineStr(:,2) = { ...
    'R34H10_Bu';...
    'R88A06_Bu';...
    'R49E09_Bu';...
    };

% store data in struct (no longer necessary)
clear resp
resp.cell.flash9 = flash9; % blue flash responses
% resp.cell.flash3 = flash3; % polarized UV flash responses (no longer
% stored)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot blue flash timeseries
savename = 'TuBu_resp';

for aidx = 1:length(areaStr)
    
    figure('color','w')
    addExportFigToolbar
    ax = gca;
    hold(ax,'on')
    
    for sidx = 1:length(lineStr(:,aidx))
        if ~isempty(lineStr{sidx,aidx})
            t1 = 1.5;
            t2 = 5.5;
            clear this*
            for midx = 1:length(maskStr)
                for eidx = [9] % blue flash only
                    for tidx = 1% blue flash only %:length(resp.(maskStr{midx}).(['flash' num2str(eidx)]).(lineStr{sidx,aidx}))
                        
                        thisTimeVec = resp.(maskStr{midx}).(['flash' num2str(eidx)]).(lineStr{sidx,aidx})(tidx).timeVector;
                        thisRespArray = [resp.(maskStr{midx}).(['flash' num2str(eidx)]).(lineStr{sidx,aidx})(tidx).flyMeanTraces];
                        
                        linecolor = objCols.(lineStr{sidx,aidx});
                        
                        lineprops.col = {linecolor};
                        lineprops.edgestyle = ':';
                        lineprops.width = 1;
                        
                        % Plot shaded error bar
                        mseb(thisTimeVec , mean(thisRespArray,1) , std(thisRespArray,[],1)./sqrt(size(thisRespArray,1)) ,lineprops, 1);
                        
                    end
                end
            end
        end
    end
    
    % Adjust plot aesthetics
    
    padStart = 2;
    ax.YLim = [-0.5 1.25];
    ax.YLabel.String = (['\Delta' 'F/F']);
    ax.XLabel.String = ('time (s)');
    
    ax.XTick = [padStart-0.5:1:padStart+3.5];
    ax.XTickLabel = {'0';'';'';'';'4 s'};
    ax.XLabel = [];
                pbaspect(ax,[1 1.5 1])
%     daspect(ax,[10 1 1])
    if ~exist('x','var') % for stim patch
        loadR49E09_AOTU
    end
    patchH = plotStimPatch(x(1).MIP,ax,[padStart-0.5 padStart+3.5],[], lightGreyCol);
    ax.YLim = [-0.25 1.05];
    setAVPaxes(ax,defaultAxisHeight_cm)
        offsetAxes(ax)
tightfig(gcf)
    
    
    suffix = areaStr{aidx};
    printAVP
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot max response to blue flash
savename = 'TuBu_max';

for aidx = 1:length(areaStr)
    
    xpos = 1:sum(~cellfun(@isempty,lineStr(:,aidx)));
    alphaVal = 0.05;%/length(xpos); % Bonferroni correction removed
    
    figure('color','w')
    ax = gca;
    hold(ax,'on')
    ax.YLim = [-0.5 1.25];
    
    for sidx = 1:length(lineStr(:,aidx))
        if ~isempty(lineStr{sidx,aidx})
            
            % Find max response amplitude during the flash stimulus
            t1 = 1.5;
            t2 = 5.5;
            clear this*
            for midx = 1:length(maskStr)
                for eidx = [9] % blue flash only
                    for tidx = 1% blue flash only %:length(resp.(maskStr{midx}).(['flash' num2str(eidx)]).(lineStr{sidx,aidx}))
                        
                        thisTimeVec = resp.(maskStr{midx}).(['flash' num2str(eidx)]).(lineStr{sidx,aidx})(tidx).timeVector;
                        thisRespArray = [resp.(maskStr{midx}).(['flash' num2str(eidx)]).(lineStr{sidx,aidx})(tidx).flyMeanTraces];
                        
                        [~,s1] = min(abs(thisTimeVec-t1));
                        [~,s2] = min(abs(thisTimeVec-t2));
                        
                        [absmax,maxidx] = nanmax(abs(thisRespArray(:,s1:s2)),[],2);
                        clear flyrespmax
                        for mm = 1:length(maxidx)
                            flyrespmax(mm) = sign(thisRespArray(mm,s1-1+maxidx(mm)))*absmax(mm);
                        end
                        this(eidx).Data(:,tidx) =flyrespmax';%[resp.(maskStr{midx}).(['flash' num2str(eidx)]).(lineStr{sidx,aidx})(tidx).maxOfFlies];
                        
                    end
                    
                end
                totalData = [this(9).Data(:)];
                
                % Add a boxplot to the axes for this group of cells
                
                totalGroup = xpos(sidx).*ones(size([this(9).Data(:)]));
                
                [bkgHandles,bkgStats] = notBoxPlot(totalData,totalGroup,'style','sdline','markMedian',true,'jitter',0.7);
                groups = unique(totalGroup);
                for bIdx = 1:length(groups)
                    
                    % Adjust boxplot aesthetics
                    
                    if ~isempty(bkgHandles(bIdx).data)
                        bkgHandles(bIdx).semPtch.FaceColor = lightGreyCol;
                        bkgHandles(bIdx).data.Color = darkGreyCol;
                        bkgHandles(bIdx).mu.Color = objCols.(lineStr{sidx,aidx});
                        bkgHandles(bIdx).med.Color = 'none';
                        
                        
                        bkgHandles(bIdx).semPtch.EdgeColor = 'none';
                        bkgHandles(bIdx).semPtch.FaceAlpha = 0.8;
                        bkgHandles(bIdx).semPtch.EdgeAlpha = 0.8;
                        
                        bkgHandles(bIdx).data.Marker = '.';
                        
                        bkgHandles(bIdx).data.MarkerSize = defaultMarkerSize;
                        
                        
                        bkgHandles(bIdx).sd.Visible = 'off';
                        
                        bkgHandles(bIdx).mu.LineWidth = 0.5;
                        bkgHandles(bIdx).med.LineWidth = 0.5;
                        bkgHandles(bIdx).med.LineStyle = ':';
                        
                        % run stats, print in command window, add asterisks
                        [~,pvalbkg] = ttest(bkgStats(bIdx).vals);
                        %                         pvalbkg = signrank(bkgStats(bIdx).vals);
                        fprintf([lineStr{sidx,aidx} ' [' char(916) 'F/Fmax]: %1.3f CI95[%1.3f %1.3f] p=%1.6f t-test, N=%d\n'],...
                            bkgStats(bIdx).mu,bkgStats(bIdx).mu-bkgStats(bIdx).interval,bkgStats(bIdx).mu+bkgStats(bIdx).interval,pvalbkg,bkgStats(bIdx).N)
                        if pvalbkg<alphaVal/length(xpos)
                            text(xpos(sidx),ax.YLim(2),'*','FontSize',8,'HorizontalAlignment','Center','BackGroundColor','none', 'Tag','asterisk','VerticalAlignment','middle');
                        end
                        
                    end
                end
            end
        end
    end
    
    
    % Adjust plot aesthetics
    
    ax.Clipping = 'off';
    ax.XTickLabelRotation = 45;
    ax.XTick = xpos;
    ax.XTickLabels = lineStr(1:length(xpos),aidx);
    ax.XAxis.TickLabelInterpreter = 'none';
    
    ax.YLabel.String = (['\Delta' 'F/F']);
ax.XLim = [xpos(1)-0.7 xpos(end)+0.7];
    %             pbaspect(ax,[1 1.5 1])
    
    ax.YGrid = 'on';
    ax.YLabel.Rotation = 90;
       daspect(ax,[6 range(ylim(ax)) 1])
       pbaspect([1,6,1])

    setAVPaxes(ax,defaultAxisHeight_cm)
    offsetAxesXTickONLY(ax)
    trimYAxisToLims(ax)
    ax.Layer = 'bottom';
    tightfig(gcf)
    % setAVPaxes(ax,defaultAxisHeight_cm)
%     tightfig(gcf)
    
    suffix = areaStr{aidx};
    printAVP
    
end