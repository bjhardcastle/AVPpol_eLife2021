% Plots using single cycle responses of E-PGS, to show the relative change
% in tuning across PB columns and to test whether there is local, momentary
% polarotopic organization of reponses

assert(exist('printpath')==1,'Specify printpath in workspace for ''printAVP''')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % load stored mat file if it exists, or try to refresh
pathsAVP
if exist(snapshotStruct_HalfCycle_path('SS00096_PB'),'file')
    load(snapshotStruct_HalfCycle_path('SS00096_PB'))
else
    try
        refresh_EPG_cycle_snapshots
    catch
        warning(['Expected to find ' snapshotStruct_HalfCycle_path('SS00096_PB')])
        disp('If raw data are available, run ''refresh_EPG_cycle_snapshots.m''')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % initial processing
% across all stimulus cycles in one recording, find a cycle where the mean
% psi across all ROIs is above a minimum threshold. Then find the
% glomerulus-pair with the highest PSI (G0,which we consider the most
% reliable tuning). For all other glomeruli, find their distance (in terms
% of glom. index / PB columns) and tuning (angle in range +/-90) relative
% to G0. [m x n] stats array returned, where (m,n) is the p-value from a
% test of all relative tunings in glom-pair (m) vs all relative tunings in
% glom-pair (n)

loadSS00096_PB

superPolThreshold(x,-2.5) % set threshold to 2.5 std above mean psi within cells in ctrl exps (no polarizer)
% This deviates from the threshold used elsewhere (see superPolThreshold)
% since the 2-cycle PSI in E-PGs was often lower than background
mdir =[];
tuning = [];
psi = [];
flyid = [];
unshiftedpsi=[];
unshiftedtuning=[];
% shuffle the order of cycles so we're not always picking the first cycle
% in time. Each exp has diferent numbers of cycles so we're also not always
% picking the first index in the shuffled order.
rng('default')
cycs = randperm(8);
cycs(cycs==1) = []; % remove first cycle, which could contain leftover 'light on' responses
for c = 1:length(cycs)
    
    cycIdx = cycs(c);
    
    for n = length(snapStruct):-1:1 % go through roi cycle data
        
        if ~isempty(snapStruct(n).individual) && size(snapStruct(n).individual.ROImeanTuningAng,1) >= cycIdx
            
            for xidx = 1:length(x)
                
                % Find the matching recording for this roi cycle data
                % Don't include ctrl exps
                % Don't include more than one cycle per fly
                if strcmp( [x(xidx).DateStr x(xidx).TimeStr] ,num2str(snapStruct(n).individual.ID)) ...
                        && x(xidx).MIP.polExp > 2 ...
                        && ~any(flyid==(snapStruct(n).individual.ID))
                    
                    thispsi = snapStruct(n).individual.ROIpsi(cycIdx,:);
                    [maxpsi,psiIdx] = max(thispsi);
                    
                    % Only consider a cycle where the mean psi across all ROIs is above a minimum threshold
                    if mean(thispsi) >=  x(xidx).MIP.polSelThreshold
                        %                     if maxpsi >=  x(xidx).MIP.polSelThreshold ...
                        %                         && all(thispsi>=0.75*maxpsi)
                        
                        abstuning = snapStruct(n).individual.ROImeanTuningAng(cycIdx,:);
                        
                        % rearrange psi and tuning arrays so the max psi index is
                        % in last position
                        shiftpsi = circshift(thispsi,-psiIdx);
                        shifttuning = circshift(abstuning,-psiIdx);
                        
                        % find the tunings relative to the max psi
                        % glomerulus and add to pooled array. (without the
                        % max psi glom, which would always have a relative tuning
                        % of 0deg)
                        tuning(end+1,:) = shifttuning - shifttuning(end) ; % range [-180:180]
                        
                        % collect some other data too
                        psi(end+1,:) = shiftpsi; % don't actually need this
                        mdir(end+1) = x(xidx).MIP.pSet(x(xidx).MIP.polExp).StepDIR;
                        flyid(end+1,:)=(snapStruct(n).individual.ID);
                        
                        unshiftedpsi(end+1,:) = thispsi;
                        unshiftedtuning(end+1,:) = abstuning;
                    end
                end
            end
        end
    end
end

% relative glomerulus positions
dist = [1 2 3 4 3 2 1 0];

% Express all tunings as the magnitude of relative shift [0:90deg] from the tuning of
% the max psi glomerlus, G0
pubVecDirArray = (deg2rad(tuning) );
q = (pubVecDirArray < -0.5*pi) | (0.5*pi < pubVecDirArray);
pubVecDirArray(q) = wrapToPi(pubVecDirArray(q) + pi) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Bar graph of probabilities of each relative tuning at each relative distance 

% bin relative tunings and get probability of each bin
hbins = [0:15:90];
N = [];
du = sort(unique(dist));
for dIdx = 1:length(du)
    [N(dIdx,:),edges]= histcounts(abs( rad2deg( pubVecDirArray(:,dist==du(dIdx)) )),hbins,'Normalization','probability');
end

% plot
f = figure('color','w');
ax = gca;
hold(ax,'on')
xpos = [1,2.5,4,5.5,7];
bH = bar(xpos,N,'stacked');
statsArray = [];
% run stats and add asterisks for reference
for dIdx = 1:length(du) % increasing distance
    statsArray(dIdx,:) = nan(1,4);
    
    % get the relative tuning angles for this glomerulus
    dd{dIdx} = abs( rad2deg( pubVecDirArray(:,dist==du(dIdx))));
    
    if dIdx > 1 % compare distribtuion to all previous glomeruli individually
        for sIdx = 1:dIdx-1
            statsArray(dIdx,sIdx) = ranksum(dd{dIdx}(:),dd{dIdx-sIdx}(:));
            if statsArray(dIdx,sIdx) >= 0.05
            sh = sigstar({[xpos(dIdx),xpos(dIdx-sIdx)]},statsArray(dIdx,sIdx) );
            end
        end
    end
end

% tweak plot
getAVPplotParams
cols = flipud(inferno(length(bH)+1));

% brown/yellow color scheme
cols = flipud(copper(length(bH)+1));
cols(:,1) = 0.9.*cols(:,1);
cols(2,:) = [0.85 0.68 0.45 ];
cols = [0.98 0.95 0.92;cols];

% gray color scheme
cols = flipud(gray(length(bH)+1));
m=length(bH)+1;
a= [255 255 255;189,189,189;99,99,99]./255;
a= [255 255 255;158,202,225;49,130,189]./255;

n = size(a,1);
% interpolate the originals to have m values as requested
% note that there are n-1 steps between 1 and n
cols = interp1(1:n,a,1:(n-1)/(m-1):n);

for bidx = 1:length(bH)
    bH(bidx).FaceColor = cols(bidx,:);
    bH(bidx).EdgeColor = darkGreyCol;
    
    bH(bidx).LineWidth = 0.5;
    
    yp = sum(N(end,1:bidx))- 0.5*N(end,bidx);
    xp = xpos(end)+1*bH(bidx).BarWidth;
    tH(bidx) = text(xp,yp,[char(8592) ' ' char(177) num2str( hbins((bidx))) ':' num2str(hbins((bidx+1))) char(176)]);
    tH(bidx).FontSize = axisLabelFontSize;
end

ax.XTick = xpos;
ax.XLim = [xpos(1)-1.2*bH(bidx).BarWidth xpos(end)+1.2*bH(bidx).BarWidth];
ax.XTickLabel = strcat(char(177),num2str(du'));
ax.XLabel.String = {'glomerulus position';'(relative to max)'};
ax.YLabel.String = 'probability';
ax.YLabel.Position(2) = 0.5;
ax.YTick = [0:0.25:1];
ax.YTickLabel = {'0','','','','1'};
% ax.YLim = [0 1];
pbaspect([1 1.75 1])
% daspect([10 1 1])
trimAxesToLims(ax)
setAVPaxes(ax,defaultAxisHeight_cm+1,defaultAxisHeight_cm+2)
set(ax,'layer','top')
uistack(ax,'top')
addlistener (ax, 'MarkedClean', @(obj,event)resetVertex(ax));

tightfig(f)
addExportFigToolbar(gcf)

prefix = 'relative_';
savename = 'EPG_PB';
suffix = 'pdf_bar';
printpath = fig8path;
printAVP
 
disp(['n=' num2str(length(tuning(:))) ' ROIs, N='  num2str(length(unique(flyid))) ' flies']);

statsArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Mean shift per column of grid plot (magnitude, unsigned)
shiftMaxPos=-4; % put G0 in the center
plottuning = circshift(rad2deg( pubVecDirArray),shiftMaxPos,2);

figure('color','w')
lineprops.width = 0.5;
lineprops.col= {darkGreyCol};
mseb([],mean(abs(plottuning),1),std(abs(plottuning),[],1),lineprops,1);
ax = gca;
ax.XTick = 1:length(dist);
ax.XTickLabel = circshift(num2str(dist'),shiftMaxPos);
ax.XLabel.String = {'glomerulus position';'(relative to max)'};
ax.YLabel.String = '||relative tuning||';
ax.YTick = [0:30:90];
ax.YTickLabel = {['0' char(176)],'','',['90' char(176)]};
ax.YLim = [0 90];

setAVPaxes(ax,[],defaultAxisHeight_cm-0.15)
offsetAxes(ax)
pbaspect([1 0.5 1])

set(ax,'layer','top')
uistack(ax,'top')

tightfig(gcf)
addExportFigToolbar(gcf)

prefix = 'mean_relative';
savename = 'PBgrid';
suffix = 'dist_unsigned';
printpath = fig8s6path;
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Grid plot with tuning vectors, colormap: signed shift distance +/- (blue/red)

shiftMaxPos = -4; % put G0 in the center
plottuning = circshift(rad2deg( pubVecDirArray),shiftMaxPos,2);
plotpsi = circshift(psi,shiftMaxPos,2);

ax = pbgridplot_distsigned(plottuning,plotpsi,num2str(mdir'),1,1);
ax.XTickLabel = circshift(num2str(dist'),shiftMaxPos); 
ax.XTickLabelRotation = 0;
setAVPaxes(gca,[],defaultAxisHeight_cm)
tightfig(gcf)

prefix = 'relative';
savename = 'PBgrid';
suffix = 'dist_signed';
printpath = fig8s6path;
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Grid plot with tuning vectors, colormap: unsigned magnitude of shift distance (greyscale)
%{
shiftMaxPos = -4; % put G0 in the center
plottuning = circshift(rad2deg( pubVecDirArray),shiftMaxPos,2);
plotpsi = circshift(psi,shiftMaxPos,2);

ax = pbgridplot_dist(plottuning,plotpsi,num2str(flyid),1,1);
ax.XTickLabel = circshift(num2str(dist'),shiftMaxPos); 
ax.XTickLabelRotation = 0;

prefix = 'relative';
savename = 'PBgrid';
suffix = 'dist_unsigned';
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Grid plot with tuning vectors, colormap: PSI 
%{
shiftMaxPos = -4; % put G0 in the center
plottuning = circshift(rad2deg( pubVecDirArray),shiftMaxPos,2);
plotpsi = circshift(psi,shiftMaxPos,2);

ax = pbgridplot(plottuning,plotpsi,num2str(flyid),1,1);
ax.XTickLabel = circshift(num2str(dist'),shiftMaxPos); 
ax.XTickLabelRotation = 0;

prefix = 'relative';
savename = 'PBgrid';
suffix = 'dist_PSI';
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Grid plot with tuning vectors, colormap: signed shift distance +/- (blue/red)
%{
shiftMaxPos = 0; % no longer using G0-relative tunings, but abs PB column index
plottuning =unshiftedtuning;
plotpsi = unshiftedpsi;

ax = pbgridplot(plottuning,plotpsi,num2str(flyid),1,1);
% ax.XTickLabels = {'8L/2R';'';'';'5L/5R';'';'';'';'1L/1R'};
% ax.XTickLabelRotation = 0;

prefix = 'regular';
savename = 'PBgrid';
suffix = 'PSI';
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % scatter plot using these half cycle data
figure('color','w')
ax = gca;
cols = flipud(magma(256));
xposArray = [1:8].*ones(size(unshiftedtuning));
scatter(xposArray(:),unshiftedtuning(:),defaultMarkerSize,cols(round(unshiftedpsi(:)*255+1),:),'filled');

ax.XLim = [0.5 8.5];
ax.XTick = [1:1:8];
ax.XTickLabelRotation = 45;
ax.XTickLabels = {'8L/2R';'';'';'5L/5R';'';'';'';'1L/1R'};
pbaspect(ax,[1 1 1]);
ax.YLim = [-90 90];
ax.YTick = [-90:30:90];
ax.YTickLabel = {['-90' char(176)];'';'';['0' char(176)];'';'';['90' char(176)]};
ax.XLabel.String = 'glomerulus pair';
ax.YLabel.String = 'preferred AoP';


% correlation, pooled recordings:
tun = deg2rad(2*unshiftedtuning);
pos = linspace(0,2*pi,17)-pi;
pos = pos(2:2:end);
pos = pos.*ones(size(tun));

[rho,pval] = circ_corrcc(tun(:),pos(:));
 [rho,perm_pval] = circ_corrcc_permute(tun(:),pos(:));
fprintf('[pooled]: r=%1.2f, p=%1.5f(circ-circ),p=%1.5f(perm), N=%d, mean PSI=%1.2f\n', ...
    rho,pval,perm_pval,size(unshiftedtuning,1),mean(unshiftedpsi(:)),std(unshiftedpsi(:)))


% correlation, individual recordings:
rhoIndiv = []; pvalIndiv = [];
for tidx = 1:size(tun,1)
    
            [rho,pval] = circ_corrcc(tun(tidx,:),deg2rad(pos(tidx,:)));
        rhoIndiv(end+1) = rho;
        pvalIndiv(end+1) = pval;
end
% for recordings with a pval less than 0.05, apply fisher
% z-transform and average
signifRhoIndiv = rhoIndiv(pvalIndiv<0.05);
% %[ try including all ]:
% signifRhoIndiv = rhoIndiv;

meanRhoIndiv = tanh(mean(atanh(signifRhoIndiv)));
semRhoIndiv = tanh(std(atanh(signifRhoIndiv)/sqrt(length(signifRhoIndiv))));

fprintf('[individual]:\n  r=%1.2f, N=%d/%d SEM=%1.2f (significant only)\n',...
    meanRhoIndiv,length(signifRhoIndiv),length(rhoIndiv),semRhoIndiv);


setAVPaxes(ax,[],defaultAxisHeight_cm)
tightfig(gcf)

suffix = 'scatter';
savename = 'EPG';
prefix = 'half_cyc_sel';
printpath = fig8path;
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resetVertex ( ax )
ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
ax.YRuler.Axle.VertexData(2,2) = max(get(ax, 'Ytick'));
end

function [rho,perm_pval] = circ_corrcc_permute(phi, x)
nreps = 10000;
showPlot = 0;

P = phi;
X = x;

perm = zeros(nreps,1);
rng default
fprintf('permutation loop')
for i=1:nreps
    if mod(i,1000)==0
        fprintf('.')
    end
    perm(i) = circ_corrcc(P,X);
    perm(i) = -perm(i); % to express in published angles (see rec2Weir)
end
fprintf('\n')

[rho,pval] =circ_corrcc(P,X);
rho=-rho; % to express in published angles (see rec2Weir)

perm_pval = (sum(abs(perm)>abs(rho)) + 1)/ (nreps + 1); % upperbound on pval

if showPlot
    figure(919)
    subplot(2,1,1)
    
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
