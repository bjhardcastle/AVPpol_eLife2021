% Plots the cross-correlation coef. of L and R glomerulus responses, under
% different L+R pairing schemes 

if ~exist('x','var') || ~strcmp(x(1).Line,'SS00096-gal4')
    loadSS00096_PB
end

limitToPolExp = 0; % Use entirety of activity recorded or just the pol mapping
applyFisherZ = 1;

posMod = [-1:6]; % array of pairing schemes to test. 
% [-1, 0, 1] are the only ones of interest, but the others put the size of
% the differences in context. Also shows a sinusoidal relationship, which
% is reassuring
% [0] corresponds to currently used pairing scheme (1L/1R - as per Wolff et al 2015 )
% -1 corresponds to a leftwards shift of the rhs ROI for a given lhs ROI (1L/8R)
% +1 a rightwards shift of the rhs ROI (1L/2R)
pairingStr = {'1L/1R';'1L/2R';'1L/3R';'1L/4R';'1L/5R';'1L/6R';'1L/7R';'1L/8R'};
pairingVal = [1:8];

corrPeaks = {};
corrGrps = {};
for oidx = 1:length(x)
    if isempty( x(oidx).MIP.ROI )
        loadROIs( x(oidx).MIP)
    end
    for m = 1:length(posMod)
        
        pairingIdx = circshift([17:24],-posMod(m));
        
        corrPeaks{m,oidx} = nan(1,8);
        corrGrps{m,oidx} = nan(1,8);
        for n = 1:8
            
            % for each glomerulus position, we get the ROI time-series,
            % plus that of its partner: the position of the partner changes
            % with posMod
            a = x(oidx).MIP.ROI(n+8).response;
            
            b = x(oidx).MIP.ROI(pairingIdx(n)).response;
            
            if limitToPolExp
                if ~isempty(x(oidx).MIP.polExp)
                   a = a(x(oidx).MIP.expStart(x(oidx).MIP.polExp) : x(oidx).MIP.expStop(x(oidx).MIP.polExp));
                   b = b(x(oidx).MIP.expStart(x(oidx).MIP.polExp) : x(oidx).MIP.expStop(x(oidx).MIP.polExp));
                end
            end
            
            % cross-corr both time-series and find the normalized coeff at
            % zero lag. All coeffs are fairly high, since we correlate the
            % responses from the entire recording, with events like LEDs
            % turning on and off present in the responses and always being
            % correlated regardless of PB position
            [r,lags] = xcorr(a,b,'coeff');
            f = find(lags==0);
            % corrPeaks{m,oidx}(n) = max(abs(r));
            corrPeaks{m,oidx}(n) = r(f);
            corrGrps{m,oidx}(n) = posMod(m);
                        
        end
    end
end
%% boxplot
figure,hold on
% deal with nans ( shouldn't be any )
grps = [corrGrps{:}];
pks = [corrPeaks{:}];
grps(isnan(pks)) = [];
pks(isnan(pks)) = [];

grpLength = length(posMod).*8.*length(x);
alphaVal = 0.05/grpLength;
[bkgHandles,bkgStats] = notBoxPlot(pks,grps,'style','sdline','markMedian',true,'jitter',0.4,'manualCI',alphaVal);

getAVPplotParams
for bIdx = 1:length(posMod)
    
    if ~isempty(bkgHandles(bIdx).data)
        bkgHandles(bIdx).semPtch.FaceColor = lightGreyCol;
        bkgHandles(bIdx).semPtch.EdgeColor = 'none';
        bkgHandles(bIdx).semPtch.FaceAlpha = 0.8;
        bkgHandles(bIdx).semPtch.EdgeAlpha = 0.8;
        bkgHandles(bIdx).semPtch.Visible = 'off';
        
        bkgHandles(bIdx).data.Marker = '.';
        bkgHandles(bIdx).data.Color = lightGreyCol;
        bkgHandles(bIdx).data.MarkerSize = defaultMarkerSize;
        bkgHandles(bIdx).data.Visible = 'on';

        bkgHandles(bIdx).mu.Color = darkGreyCol;
        bkgHandles(bIdx).med.Color = darkGreyCol;
        bkgHandles(bIdx).med.Visible = 'on';
        bkgHandles(bIdx).sd.Visible = 'off';
        bkgHandles(bIdx).sd.Color = darkGreyCol;
        bkgHandles(bIdx).sd.LineWidth = 0.5;
        bkgHandles(bIdx).mu.LineWidth = 0.5;
        bkgHandles(bIdx).med.LineWidth = 0.5;
        bkgHandles(bIdx).med.LineStyle = ':';
        
        % set sd to top layer 
        bkgHandles(bIdx).sd.ZData  = abs(bkgHandles(bIdx).sd.ZData);
                bkgHandles(bIdx).mu.ZData  = abs(bkgHandles(bIdx).mu.ZData);
        bkgHandles(bIdx).med.ZData  = abs(bkgHandles(bIdx).med.ZData);

        if applyFisherZ
            zMean = tanh(mean(atanh(bkgStats(bIdx).vals)));
             bkgHandles(bIdx).mu.YData = [zMean zMean];
        end
    end
end
uistack([bkgHandles(:).sd],'top')
uistack([bkgHandles(:).mu],'top')
uistack([bkgHandles(:).med],'top')

ax = gca;
xlim([min(posMod)-0.4 max(posMod)+0.5])
ylim([ min([bkgStats(:).mu])-1*std([bkgStats(:).mu])   1])
% ylim([0.98,1])
ax.YLabel.String = 'norm. coeff.';
ax.XLabel.String = 'PB glomerulus-pairing scheme';
ax.XTick = posMod;
ax.XTickLabel = pairingStr(pairingVal(mod(posMod,8)+1));
ax.XTickLabelRotation = 45;

offsetAxesXTickONLY(ax)
setAVPaxes(ax)
ax.YGrid = 'on';
pbaspect([2.5,1,1])
% uistack(gca,'bottom');
% ax.Layer = 'bottom';

set(gcf,'color','w')
tightfig(gcf)
addExportFigToolbar(gcf)

