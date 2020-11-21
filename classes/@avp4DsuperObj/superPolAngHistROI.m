function varargout = superPolAngHistROI(objarray,useTuningColors,usePSIthreshold)
%SUPERPOLANGHIST(objarray,useTuningColors) Simple tuning angle histogram
% Plots half a polar histogram showing the probability of pixels within the
% MIP layer mask having a certain tuning. Adds a p-value for a Rayleigh
% test.
%
% useTuningColors       0 (default) or 1 to enable color coding of bins

if nargin <2 || isempty(useTuningColors)
    useTuningColors = 0;
end
if useTuningColors
    vers = 'color'; % Apply hsv color map to each histogram bin
else
    vers = 'mono';
end
if nargin <3 || isempty(usePSIthreshold)
    usePSIthreshold = 0;
end

getAVPplotParams
superUseMSP(objarray,1) % Use MSP instead of MIP

% histogram/histcounts can't wrap 0/180, and since we want bin centers to
% fall on presented angles (of which 0/180 is one), we instead have to
% shift all angle data between [0:0.5*binWidth] to [180:180 + 0.5*binWidth]
binWidth =15;
binCenters = [0:binWidth:180-binWidth];  % length n
shiftedBinEdges = [0.5*binWidth:binWidth:180+0.5*binWidth];  % length n+1, first/last are equivalent
% all angle data must be in the range [shiftedBinEdges(1) : shiftedBinEdges(end)]

% Initialize array for storing histogram counts
ctArray = nan(length(objarray),length(binCenters));
ROIcount = 0;

flyArray = nan(length(objarray),1);
for oidx = 1:length(objarray)
    if objarray(oidx).containsPolMapExp>=4 ...
            && ( objarray(oidx).MIP.UseMSP || objarray(oidx).numZPlanes==1 )
        
       tun = nan;
res = nan;
ct=0;
loadROIs(objarray(oidx).MIP)
for ridx = 1:length(objarray(oidx).MIP.ROI)
    [thisTh,thisR] = tuningROI(objarray(oidx).MIP,objarray(oidx).MIP.ROI(ridx).mask);
    thisPSI = mean(objarray(oidx).MIP.polSelImg(objarray(oidx).MIP.ROI(ridx).mask));
    if ~usePSIthreshold || ( usePSIthreshold && thisPSI > objarray(oidx).MIP.polSelThreshold )
        ct = ct+1;
        tun(ct) = thisTh;
        res(ct) = thisPSI;
     end
end

if length(res)>2 % if we use recordings with only one or two ROIs their data will skew the overall probabilities
            shiftAng = 0.5.*(wrapTo360((2.*( tun - shiftedBinEdges(1)))))  + shiftedBinEdges(1);
            ctArray(oidx,:) = histcounts(shiftAng(~isnan(shiftAng)),shiftedBinEdges, 'Normalization', 'probability');
            flyArray(oidx) = objarray(oidx).Fly;
            ROIcount = ROIcount + length(res);
        end
    end
end
% Get number of flies across all returned trials
included_flies = unique(flyArray(~isnan(flyArray)));
populationSize = length(included_flies);
if ~populationSize
    return
end
% Get individual fly mean-traces
fly_mean_traces = nan(populationSize,size(ctArray,2));
for ifidx = 1:populationSize
    fly_mean_traces(ifidx,:) = nanmean( ctArray( (flyArray == included_flies(ifidx)), : ) , 1);
end

% Work on these mean-traces to get SEM and the population mean trace
varTrials = nanstd(fly_mean_traces , [] , 1) ./ sqrt(populationSize);
meanTrials = nanmean(fly_mean_traces,1);

plotTrials = sqrt(meanTrials); % this way the area of each bin wedge will be proportional to the data.. can no longer read off radial values from plot though

% plot
f = figure('color','w');
ax = polaraxes;
hold(ax,'on')
plotBins = unique( [shiftedBinEdges,shiftedBinEdges+180] );
switch vers
    case 'mono'
        % regular version
        histHandles = polarhistogram('BinEdges',deg2rad(plotBins),'BinCounts',([meanTrials'; meanTrials']));
        histHandles.FaceColor = darkGreyCol;
        histHandles.FaceAlpha = 0.8;
        histHandles.EdgeColor = 'w';
        histHandles.EdgeAlpha = 0.8;
        histHandles.LineWidth = 0.5;
        
    case 'color'
        % fancy color version
        cols = flipud(hsv(length(meanTrials)));
        % plot each bin individually in the same axes (couldn't figure out
        % how to change bin face colors individually)
        for pidx = 1:2*length(meanTrials)
            histHandles(pidx) = polarhistogram('BinEdges',[deg2rad([plotBins(pidx),plotBins(pidx+1)])],'BinCounts',(meanTrials(mod(pidx-1,length(meanTrials))+1)));
            histHandles(pidx).FaceColor = cols(mod(pidx-1,length(meanTrials))+1,:);
            histHandles(pidx).FaceAlpha = 1;
            histHandles(pidx).EdgeColor = 'w';
            histHandles(pidx).EdgeAlpha = 0.8;
            histHandles(pidx).LineWidth = 0.5;
        end
end

% adjust angle representation to match Weir,Henze 2016
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'left';
ax.ThetaTickLabel = {['  90' char(176)],'','',['  0' char(176)],'','',['-90' char(176)]};
ax.ThetaLim = [0 180];
ax.ThetaTick = [0:30:180];

% ax.RLim = [0 0.2];
if max(meanTrials) > max(ax.RLim)
    disp('Data may be clipped: reduce RLim(2)')
end
% ax.RTick = [0.1:0.1:max(ax.RLim)]; % turned off since these numbers are no longer
% useful - they indicate probability^2
ax.RTick = [];

setAVPaxes(ax,[],defaultAxisHeight_cm)
set(ax,'layer','bottom')

% % Rayleigh test: h0 is uniform distribution around circle
% % Since we have axial data we double all the angles before testing
doubledBinCenters = deg2rad( 2.*(shiftedBinEdges + shiftedBinEdges(1)) );
ctVec = length(meanTrials).*reshape(fly_mean_traces',[],1); %so that N samples = numFlies*numAngleBins;
% % alternative using mean values replicated numFlies times:
%ctVec = length(meanTrials).*repmat(meanTrials,1,size(fly_mean_traces,1))'; %
angVec = repmat(doubledBinCenters(1:end-1),1,size(fly_mean_traces,1))';
[pVal] = circ_rtest(angVec,ctVec,deg2rad(2.*binWidth));

% add CI95 bar
[r,th] = circ_rangle(angVec, ctVec, deg2rad(2.*binWidth), 1);
th = th/2; % revert back to axial data
t = circ_confmean(angVec, 0.05, ctVec, deg2rad(2.*binWidth), 1);
t = t/2; % revert back to axial data
ciTh = [th-t:0.01:th+t];
ciR = repmat(0.95*ax.RLim(2),1,length(ciTh));

if ~isnan(t)
% draw an arc around the mean angle
arcHandle = polarplot(ciTh,ciR,'LineWidth',1.5,'Color','k');
% draw the diametrically oppoiste arc too, in case it wraps at 0/180
arcHandle(2) = polarplot(ciTh+pi,ciR,'LineWidth',1.5,'Color','k');
% keep it at the outside, in case we change the limits
addlistener (ax, 'MarkedClean', @(object,event)setCImarker(ax,arcHandle));
% % draw a stalk from the center
% polarplot([th th],[0 ax.RLim(2)],'LineWidth',1,'Color','k');
end

% % add stats info to plot for convenience:
% text(ax,90,-ax.RLim(2)+0.05,{['p=' num2str(pVal)];[num2str(r,'%.2f') ', ' num2str(rad2deg(th),'%.2f') char(176) ];[ 'CI:' num2str(rad2deg(t),'%.2f') char(176)]},'FontSize',axisLabelFontSize,'Color',[0.6 0.6 0.6]);
% % print instead:
fprintf(['\n[hist|ROI|PSIthresh=%d)]\n N=%d, %d ROIs, p=%1.5f Rayleigh\n resultant:%1.2f %1.2f' char(176) ', CI:%1.2f' char(176) '\n\n'],...
usePSIthreshold,...
populationSize,...
ROIcount,...
pVal,...
r,...
rad2deg(th),...
rad2deg(t) ...
);

tightfig(f)
addExportFigToolbar(gcf)


if nargout >0 
    varargout{1} = histHandles;
end
if nargout > 1
    varargout{2} = meanTrials;
end
%%
% old version with a mean line +/- sem using polygonplot
%{

figure
clear opt_axes opt_area opt_lines
   polcols = flipud(hsv(180));
    opt_area.Color = darkGreyCol;%polcols(ROIang,:);
    opt_area.err = 'sem';
    opt_area.FaceAlpha = 0.5;
    
    
    opt_lines.LineWidth = 1;
    opt_lines.Color = darkGreyCol;%polcols(ROIang,:);
%     opt_axes.Ticks = [0:0.25:1];
    opt_axes.NRadii = 12;
    opt_axes.Labels     = {['  -90' char(176)],'','','','','',['  90' char(176)],'','',['  0' char(176)],'',''};
%     opt_axes.Labels     = {['  -90' char(176)],'','',['180' char(176)],'','',['90' char(176) ' '],'','',['0' char(176)],'',''};
    clear alldata
    alldata(:,1,:) = [fly_mean_traces'; fly_mean_traces'] ;
%     alldata = alldata./(max(alldata,[],1)); % normalize each fly trace
   
    

    polygonplot(alldata,opt_axes,opt_lines,opt_area);
    hold on
    
    %% individual fly histograms
    if plotIndividFlies || plotIndividExps
    indivdata(:,:,1) = [fly_mean_traces'; fly_mean_traces'] ;
%     indivdata = indivdata./(max(indivdata,[],1)); % normalize each fly trace
    
    polygonplot(indivdata);
    end
        %%
    ylim([0 1])
        xlim([-1 1])
        pbaspect([ 1 0.5 1])
        %%

%         % Add mean tuning pref line:
%         pixPolWeight = mean(squeeze(alldata),2)';
%         polAngs = wrapTo180(-circshift([30:30:360],-5)-270);
%         [r,axtheta] = circ_rangle(circ_axial(deg2rad(polAngs),2), pixPolWeight, deg2rad(mode(diff(polAngs))), 2);
%         theta = axtheta/2;
%         Ang =0.5.*(wrapTo180((2.*(rad2deg(theta)))));
%         hold on
%         R  = [0; 1];
%         TH = [deg2rad(Ang-90) deg2rad(Ang-270)].*ones(2,2);
%
%         [X,Y] = pol2cart(TH, R);
%
%          for n = 1:2,  plot(X(2,:),Y(2,:),'color',polcols(ROIang,:),'LineWidth',2,'LineStyle',':'); end
%
setAVPaxes(gca,2.5)
tightfig(gcf)

%}

end
function setCImarker(pax,arcHandle)
upper  = pax.RLim(2);
for aIdx = 1:length(arcHandle)
    arcHandle(aIdx).RData = ones(size(arcHandle(aIdx).RData)).*upper;
end
end