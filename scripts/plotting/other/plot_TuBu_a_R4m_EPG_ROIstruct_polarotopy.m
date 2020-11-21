% plot ROI-based linear and circular position vs tuning plots in the bulb
% or PB from stored .mat files

% requires 'lineStr' variable to be already set to one of the following:
%     'R34H10_Bu'
%     'R34D03_Bu'
%     'SS00096_PB'

% Note:  This E-PG section has been superceded by another scatter plot
% Here, the PSI values used were calculated from first 2 cycles of each pol
% mapping experiment. We later understood that these weren't reliable
% (indicated by the low PSI) and sought single cycles where all glomeruli
% contained reliable responses, and then plot their tuning in a scatter
% plot. The results turned out to be qualitatively similar (no correlation
% between position and tuning) but given the inconsistencies in E-PG
% responses, the later version is the more correct version (see
% 'make_various_EPG_panels.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % load stored mat file if it exists, or try to refresh

assert(exist('lineStr')==1 && ismember(lineStr,{'R34H10_Bu';'R34D03_Bu';'SS00096_PB'}),...
    ['Make sure lineStr, savePlotStr and printpath are in workspace (see ' mfilename ')']);

pathsAVP
if exist(ROIstruct_path(lineStr),'file')
    load(ROIstruct_path(lineStr))
else
    try
        refresh_TuBu_a_R4m_EPG_ROIstruct
    catch
        warning(['Expected to find ' ROIstruct_path(lineStr)])
        disp('If raw data are available, run ''refresh_TuBu_a_R4m_EPG_ROIstruct.m''')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot all scatter plots and then save
getAVPplotParams

cols = [0 0 0; flipud(magma(255))];

switch lineStr
    
    case 'SS00096_PB'
        
        clear ax f
        
        % right and left data are available in ROIstruct, but almost
        % identical, so we just plot the joint data (L+R joint ROIs)
        
        f(1) = figure('color','w');
        ax(1) = axes(f(1));
        
        hold(ax(1),'on')
        
        ax(1).XLim = [0.5 8.5];
        ax(1).XTickLabelRotation = 45;
        
        % ax(1).XTick = [1:1:8];
        % ax(1).XTickLabels = {'1L/2R';'2L/3R';'3L/4R';'4L/5R';'5L/6R';'6L/7R';'7L/8R';'8L/1R'};
        ax(1).XTick = [1:1:8];
        ax(1).XTickLabels = {'1L/2R';'';'';'';'5L/6R';'';'';'8L/1R'};
        
        for sidx = 1:length(ROIstructJoint)
            scatter(ax(1),[1:8],ROIstructJoint(sidx).tuning,defaultMarkerSize,cols(round([ROIstructJoint(sidx).PSI]*255+1),:),'filled');
        end
        
        joint_pos = repmat([1:8], length([ROIstructJoint.tuning])/8,1);
        joint_pos = joint_pos.';
        
        [rho,pval] = circ_corrcc(deg2rad([ROIstructJoint.tuning]),deg2rad(0.5.*[joint_pos(:)]));
        title(ax(1),sprintf('pooled: r=%1.2f, p=%1.5f, N=%d',...
            rho,pval,size([ROIstructJoint.Fly],2)));
        
        for aidx = 1:length(ax)
            pbaspect(ax(aidx),[1 1 1]);
            ax(aidx).YLim = [-90 90];
            ax(aidx).YTick = [-90:30:90];
            ax(aidx).YTickLabel = {['-90' char(176)];'';'';['0' char(176)];'';'';['90' char(176)]};
            ylabel(ax(aidx),'preferred AoP')
            uistack(ax(aidx),'bottom')
            xlabel(ax(1),'normalized position')
            setAVPaxes(ax(aidx),2)
            tightfig(f(aidx))
            
        end
        
        prefix = 'scatter';
        savename = lineStr;
        suffix = 'joint';
        printAVP
        
    otherwise
        
        clear ax f
        
        % all possible polarotopic organizations to test
        allPlotStr = { ...
            'horiz_right';...
            'horiz_left';...
            'horiz_all';...
            'vert_right';...
            'vert_left';...
            'vert_all';...
            'circ_right';...
            'circ_left';...
            'circ_all';...
            };
        
        % selected polarotopic organizations to save
        if ~exist('savePlotStr','var') || isempty(savePlotStr) || ~all(contains(savePlotStr,allPlotStr))
            savePlotStr = {};
        elseif ~iscell(savePlotStr)
            savePlotStr = {savePlotStr};
        end
        
        % first make scatter plots for all organizations and print stats, then decide
        % whether to save or close figure
        for pidx = 1:length(allPlotStr)
            
            f(pidx) = figure('color','w');
            ax(pidx) = gca;
            hold(ax(pidx),'on')
            
            % First collect data for this plot type:
            plotStr = allPlotStr{pidx};
            
            if contains(plotStr,'right')
                ROIstruct = ROIstructRight;
                statsStr = 'Right';
            elseif contains(plotStr,'left')
                ROIstruct = ROIstructLeft;
                statsStr = 'Left';
            elseif contains(plotStr,'all')
                ROIstruct = ROIstructPooled;
                statsStr = 'Pooled';
            end
            
            if contains(plotStr,'horiz')
                posStr = 'posX';
            elseif contains(plotStr,'vert')
                posStr = 'posY';
            elseif contains(plotStr,'circ')
                posStr = 'posAng';
            end
            
            % store as cells so we can run stats on each individual
            % recording later
            posValCell = {ROIstruct.(posStr)};
            tuningValCell = {ROIstruct.tuning};
            psiValCell = {ROIstruct.PSI};
            
            % pool values for plotting/pooled stats
            posVals = [posValCell{:}];
            if ~contains(plotStr,'circ')
                % position data are in the range [0:1], normalized to the
                % bounds of the layer mask. Here we normalize the data to
                % their own range, so [0:1] are the min and max values
                % (linear position data only)
                posVals = ( [posValCell{:}] - min([posValCell{:}]) ) ./ ( max([posValCell{:}]) - min([posValCell{:}]) );
            end
            tuningVals = [tuningValCell{:}];
            psiVals = [psiValCell{:}];
            statsVals =  ROIstats.(statsStr);
            
            % plot
            scatter(ax(pidx),posVals,tuningVals,defaultMarkerSize,cols(round(psiVals*255+1),:),'filled');
            
            
            
            % correlation, all ROIs pooled across recordings:
            
            if contains(plotStr,'circ')
                % circular-circular correlation:
                [rho,pval] = circ_corrcc(deg2rad(tuningVals),deg2rad(0.5.*posVals));
                
                % get straight line fit to add to plot:
                % first rescale position data to [0:1]
                circPlotPosVals = 0.5*(1+(posVals/180));
                [~,~,phi_0,amax] = circlin_wcorr(deg2rad(2.*tuningVals),circPlotPosVals,[-2 2]);
                
            else
                % linear-circular correlation:
                [rho,pval,phi_0,amax] = circlin_wcorr(deg2rad(2.*tuningVals),posVals,[-2 2]);
            end
            
            
            [~,perm_pval] = circlin_wcorr_permute(deg2rad(2.*tuningVals), posVals, [-2 2]);
            %{
        % add fit line for pooled recordings:
        if perm_pval<0.05
            [lineHandle, phi_fit_returned] = plotPhi_VentDors(posStr,amax,[],phi_0,ax(pidx),1);
        end
            %}
            
            % print stats for pooled recordings:
            fprintf(['[%s| %s | pooled]:\n  r=%1.2f, p=%1.5f, N=%d\n'...
                'ROIs: count=%d, mean=%1.2f' char(177) '%1.2fstd per fly\n'...
                'PSI: mean=%1.2f' char(177) '%1.2f\n'],...
                lineStr,...
                allPlotStr{pidx},...
                rho,...
                perm_pval,...
                statsVals.numFlies,...
                statsVals.numROIsTotal,...
                statsVals.numROIsPerFlyAvg,...
                statsVals.numROIsPerFlyStd,...
                mean(psiVals),...
                std(psiVals)...
                );
            %}
            
            
            % correlation, individual recordings:
            rhoIndiv = []; pvalIndiv = [];
            for ridx = 1:length(ROIstruct)
                if ~isempty(ROIstruct(ridx).PSI)
                    if ~contains(plotStr,'circ')
                        [rho,pval] = circ_corrcc(deg2rad(tuningValCell{ridx}),deg2rad(0.5.*posValCell{ridx}));
                        circPlotPosVals = 0.5*(1+(posValCell{ridx}/180));
                        [~,~,phi_0,amax] = circlin_wcorr(deg2rad(2.*tuningValCell{ridx}),circPlotPosVals,[-2 2]);
                    else
                        [rho, pval,phi_0,amax] = circlin_wcorr(deg2rad(2.*tuningValCell{ridx}),posValCell{ridx},[-2 2]);
                    end
                    rhoIndiv(end+1) = rho;
                    pvalIndiv(end+1) = pval;
                    
                    % % Add fit lines for individual recordings:
                    %[lineHandle, phi_fit_returned] = plotPhi_VentDors(posStr,amax,[],phi_0,ax(pidx),1);
                end
            end
            % for recordings with a pval less than 0.05, apply fisher
            % z-transform and average
            signifRhoIndiv = rhoIndiv(pvalIndiv<0.05);
            %[ try including all ]:
            %                       signifRhoIndiv = rhoIndiv;
            
            meanRhoIndiv = tanh(mean(atanh(signifRhoIndiv)));
            semRhoIndiv = tanh(std(atanh(signifRhoIndiv)/sqrt(length(signifRhoIndiv))));
            
            fprintf('[%s| %s | indiv]:\n  r=%1.2f, SEM=%1.2f, N=%d/%d\n',...
                lineStr,allPlotStr{pidx},meanRhoIndiv,semRhoIndiv,length(signifRhoIndiv),length(rhoIndiv));
            
            
            
            
            
            % adjust common plot aesthetics
            pbaspect(ax(pidx),[1 1 1]);
            ax(pidx).YLim = [-90 90];
            ax(pidx).YTick = [-90:30:90];
            ax(pidx).YTickLabel = {['-90' char(176)];'';'';['0' char(176)];'';'';['90' char(176)]};
            ylabel(ax(pidx),'preferred AoP')
            uistack(ax(pidx),'bottom')
            xlabel(ax(pidx),'normalized position')
            setAVPaxes(ax(pidx),2)
            if contains(plotStr,'horiz')
                ax(pidx).XLim = [0 1];
                ax(pidx).XTick = [0:0.25:1];
                if contains(plotStr,'left')
                    ax(pidx).XTickLabel = {'lateral';'';'';'';'medial'};
                else
                    ax(pidx).XTickLabel = {'medial';'';'';'';'lateral'};
                end
            elseif contains(plotStr,'vert')
                ax(pidx).XLim = [0 1];
                ax(pidx).XTick = [0:0.25:1];
                ax(pidx).XTickLabel = {'ventral';'';'';'';'dorsal'};
            elseif contains(plotStr,'circ')
                ax(pidx).XLim = [-180 180];
                ax(pidx).XTick = [-180:90:180];
                ax(pidx).XTickLabel = {['-' char(960)];'';'';'';char(960)};
            end
            
            
            tightfig(f(pidx))
            [~,plotIdx] = ismember(allPlotStr(pidx),savePlotStr);
            if plotIdx
                prefix = 'scatter';
                suffix = savePlotStr{plotIdx};
                savename = [lineStr '_ROI'];
                printAVP
            elseif ~isempty(savePlotStr)
                close(f(pidx))
            else
                f(pidx).Name = [allPlotStr{pidx} ' - ' lineStr];
            end
            
        end
        
end

function [rho,perm_pval] = circlin_wcorr_permute(phi, x, bounds)
% with help from https://courses.washington.edu/matlab1/Bootstrap_examples.html#1
showPlot = 0;

nreps = 10000;

P = phi;
X = x;


B = [-max(abs(bounds)) max(abs(bounds))];
perm = zeros(nreps,1);
rng default
if nreps > 1000
    fprintf('permutation loop')
end
for i=1:nreps
    if nreps > 1000 && mod(i,1000)==0
        fprintf('.')
    end
    %shuffle the lsat scores and recalculate the correlation
    perm(i) = circlin_wcorr(P,X(randperm(length(X))),B);
    perm(i) = -perm(i); % to express in published angles (see rec2Weir)
end
fprintf('\n')

[rho,pval,~,~] = circlin_wcorr(P,X,bounds);
rho=-rho; % to express in published angles (see rec2Weir)

perm_pval = (sum(abs(perm)>abs(rho)) + 1)/ (nreps + 1); % upperbound on pval

if showPlot
    figure(919)
    
    histbins = [-1:0.05:1];
    h = histogram(perm,histbins);
    h.FaceColor = [0.7 0.7 0.7];
    h.EdgeAlpha = 0;
    
    ylim = get(gca,'YLim');
    hold on
    l1 = plot(rho*[1,1],ylim,'r-','LineWidth',2);
    title(sprintf('P = %5.5f (%d/%d)',perm_pval, sum(abs(perm)>abs(rho)),nreps));
    xlim([-1 1])
    pbaspect([2 1 1])
    
end

end

function [lineHandle, phi_fit_returned] = plotPhi_VentDors(posstr,amax,positionData,phi_0,ax,convertAngles4Pub)

colstr = 'k';
acolstr = 'w';

lowerLim = -90;
upperLim=90;

if strcmp(posstr,'posAng')
    positionData = [-180:1:180];
    amax = amax/360;
else
    positionData = [0:0.01:1];
end

if nargin < 5 || isempty(ax),ax=gca;end
n=0;

phi_fit_returned = 180*amax*positionData + rad2deg(phi_0/2);

phi_fit_plot = 180*amax*positionData + rad2deg(phi_0/2);%;

hold(ax,'on')
for phi_shift =  [-720:180:720]
    
    this_phi_fit = phi_fit_plot  + phi_shift;
    h = plot(ax,positionData(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),this_phi_fit(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),colstr);
    if ~isempty(h)
        n=n+1;
        lineHandle(n) = h;
    end
end


end


