function varargout = superTuningCurve(objarray,useCellMask)
% Generates plot with normalized tuning curves (for angle of polarization)
% using per pixel values. we shift to center all pixel's tuning curves on
% their pref angle before pooling and plotting. using all available angle
%resp data (ie typically 12 orientations), not pooled (6 unique angles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getAVPplotParams
superUseMSP(objarray,1) % Use MSP instead of MIP
if nargin<2
useCellMask = 1;
end
useShadedErrorBar = 0;
addMeanLine = 1;
superPolThreshold(objarray,-1)
inclObj = [];
tuningCurve = [];
for oidx = 1:length(objarray)
    if ~isempty(objarray(oidx).MIP) && ~isempty(objarray(oidx).MIP.polExp) && ismember( objarray(oidx).MIP.polExp,[4,8,10] )
        if isempty(objarray(oidx).MIP.polTuningImg)
            getPolMaps(objarray(oidx).MIP);
        end
        
        % get vector of above-threshold pixels
        if ~useCellMask
%             threshPix = find(~isnan(objarray(oidx).MIP.polPix));
            threshPix = find(objarray(oidx).MIP.layerMask.mask == 1);

        elseif useCellMask
            threshPix = find(objarray(oidx).MIP.cellMask==1);
        end
        
        % Get avg resp of each pixel to the 12 angles presented
        threshResp = objarray(oidx).MIP.polAngResp(threshPix,:);
        % For normalizing we'll use the pix values from the
        % inactivity frame (spontaneous activity)
        threshBaselineVals = objarray(oidx).MIP.InactivityFrame(threshPix);
        
        % normalize responses to baseline F vals
        threshCurves = threshResp./threshBaselineVals -1;
        
        % Get preferred angle of each pixel (mean resultant angle - can take any value, not just the 12 orientations presented)
        threshPrefAng = objarray(oidx).MIP.polTuningImg(threshPix);
        
        % interpolate in Fourier domain then shift tuning angles to be
        % centered on the tuning peak for each pix
        m = interpft(threshCurves,360,2);
        M = circshift(m,360/size(threshCurves,2),2); % first angle in polangresp is 30deg - or equivalent - regardless of motor direction
        
        C = [];
        for midx = 1:length(threshPix)
            C(midx,:) = circshift(M(midx,:),-round(threshPrefAng(midx))-90);% shift by an extra 90deg for aesthetics
        end
        tuningCurve(end+1,:) = nanmean(C,1);
        inclObj(end+1) = oidx;
        
    end
end

% avg tuning curves which came from the same fly, if any exist
uniqueFlies = unique([objarray(inclObj).Fly]);
tuningCurves = [];
for uidx = 1:length(uniqueFlies)
    tuningCurves(uidx,:) = mean( tuningCurve([objarray(inclObj).Fly] == uniqueFlies(uidx) ,:) ,1);
end

if ~useShadedErrorBar
    %% Plot all normalized responses and their mean
%     printpath = fig8path;
%     savename = 'TuBu_R4m';
%     
%     prefix = 'normalized_resp';
    
    figure('color','w')
    hold on
    for sidx = 1:size(tuningCurves,1)
        
        linecolor = ROIcolor(1);
        brightcolor = linecolor+(1-linecolor)*0.55;
        brightcolor = linecolor;
        % set transparency:
        brightcolor(4) = 0.2;
        p(sidx) = plot(tuningCurves(sidx,:),'color',brightcolor,'LineWidth',0.5);
    end
    % pm = plot(mean(tuningCurve,1),'color',linecolor);
    % uistack(pm,'top')
    if addMeanLine
        P = plot(nanmean(tuningCurves,1),'color',linecolor,'LineWidth',1);
    end
    addExportFigToolbar
    % ylim([-1 1.5])
    xlim([0 360])
    xticks([0:90:360])
    xticklabels({'';'\Psi';'';['\Psi' char(177) '180\circ'];''})
    daspect([360,6,1])
    xlabel('AoP');
    ylabel(['\Delta' 'F/F']);
    offsetAxes(gca)
    scalebarF(gca)
    setAVPaxes(gca,[],2)
    tightfig(gcf)
    suffix = 'traces';
    % printAVP
    
    %%
else
    
    %% Plot mean normalized response with a shaded error bar
    % printpath = fig8path;
    % savename = 'TuBu_R4m';
    %
    % prefix = 'normalized_resp';
    %
    figure('color','w')
    hold on
    % for sidx = 1:2
    lineprops.col = {[ROIcolor(1)]};
    lineprops.width = 1;
    p = mseb([],mean(tuningCurves,1),std(tuningCurves,[],1)./sqrt(length(uniqueFlies)),lineprops,1);
    % end
    addExportFigToolbar
    % ylim([-0.2 1.2])
    xlim([0 360])
    xticks([0:90:360])
    xticklabels({'';'\Psi';'';['\Psi' char(177) '180\circ'];''})
    daspect([360,6,1])
    xlabel('AoP');
    ylabel(['\Delta' 'F/F']);
    offsetAxes(gca)
    scalebarF(gca)
    setAVPaxes(gca,[],2)
    tightfig(gcf)
    suffix = 'errorbar';
    % printAVP
    %%
end
%% print info
% disp(['[TuBua (blue)]: ' num2str(size(rs{1},1)) ' ROIs, N=' num2str(numFlies(1))])
% disp(['[R4m (red)]: ' num2str(size(rs{2},1)) ' ROIs, N=' num2str(numFlies(2))])
if nargout
    output = struct;
    if ~useShadedErrorBar
        output.Lines = p;
        if addMeanLine
           output.meanLine = P; 
        end
    else
        output.mseb = p;
    end
    varargout{1} = output;
end

end