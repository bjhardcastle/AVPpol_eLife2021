function superXY(obj,Yoff,Xoff,rotateOFF,darkMode,applyTuningCols,plotIndividFits,plotPooledFit,UseMSPifLayerMasksEmpty)
% For AVP super objects: (superclass version of plotXYcorr)
% Take pol selective pixels from the tuning map and plot x-position or
% y-position vs tuning angle. Find the circular-linear correlation and plot
% the line and R-value
if nargin<2 || isempty(Yoff) || ~any(Yoff==0 || Yoff==1)
    Yoff = 0; % set to 1 to exclude the y-correlation plot
end
if nargin<3 || isempty(Xoff) || ~any(Xoff==0 || Xoff==1)
    Xoff = 0; % set to 1 to exclude the x-correlation plot
end
if nargin<4 || isempty(rotateOFF) || ~any(rotateOFF==0 || rotateOFF==1)
    rotateOFF = 0;
end
if nargin<5 || isempty(darkMode) || ~any(darkMode==0 || darkMode==1)
end

if nargin < 6 || isempty(applyTuningCols) || ~any(applyTuningCols==0 || applyTuningCols==1)
    applyTuningCols = 0; % Color scatter points according to pol tuning ( redundant info)
end
if nargin < 7 || isempty(plotIndividFits)
    plotIndividFits = 1;
end
if nargin < 8 || isempty(plotPooledFit)
    plotPooledFit = ~plotIndividFits;
end
if nargin < 9 || isempty(UseMSPifLayerMasksEmpty)
    UseMSPifLayerMasksEmpty = 1;
    % 0, whole frame will be used if Layer layerMasks are empty.
    % 1, we revert to MSP instead.
    % 2, force MSP to be used
end

applyWeights =0; % Weighted correlation using pol selectivity for each pixel
if ~applyTuningCols
    applyPSICols = 1; % If zero, a solid color is applied to all points
end
applyPtsLimit = 1; % Plot only a subset of data pooled across recordings (too many points otherwise)
ptsLimit = 1000;

ptSize = 2; % size of dots in scatter plot

convertAngles4Pub = 1; % original polarizer angles are used for correlation
% calculations etc, but scatter plots modify the angles to fit Weir,Henze
% et al 2016

runPermutationTest = 0; % shuffles data and recalculates correlation multiple times

corrStatMethodPooled = 0; % 0' regular corr coef | 1* bootstrap corr coef (all data) | 2 same as 1 (smaller subset)
corrStatMethodIndiv = 5;  % 5' Fisher mean coef, SEM  | 4'* Fisher mean with bootstrap CI | 3' regular mean with bootstrap CI
                          %  'quick   *preferred 
                          if corrStatMethodPooled~=1 || corrStatMethodIndiv~=4
                              warning('superXY: running in fast mode, stats may be slightly different to published versions')
                          end
showPlot = 0;
if showPlot
        figure(919),clf,cla
end

limitToSignifRecs = 1; % recordings with only a couple of
% pixels have a meaningless regression line/coeff so we won't plot them or
% analyze their individual coefficients, which have p>0.05
% However their data are still pooled for an analysis across all recordings
%

MIPoff =1; %  if on a tuning map is plot alongside the scatter plots (was useful to check rotations)

darkMode = 0; % was used to plot black bkgrnd/white text, now defunct

getAVPplotParams

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup axes for plotting

figure
addExportFigToolbar(gcf)
if darkMode
    colstr = 'w';
    acolstr = 'k';
    greyFitLines = [.5 .5 .5];
    meanFitLine = [1 1 1];
else
    colstr = 'k';
    acolstr = 'w';
    greyFitLines = [0 0 0];
    meanFitLine = [.7 .7 .7];
    
end
set(gcf,'color',acolstr)


% Plot x- correlation
if ~Xoff
    if ~Yoff && ~MIPoff
        ax(1) = subplot(1,3,1);
    elseif Yoff || MIPoff
        ax(1) = subplot(1,2,1);
    elseif Yoff && MIPoff
        ax(1) = subplot(1,1,1);
    end
    hold(ax(1),'on')
end
% Plot y- correlation
if ~Yoff
    if ~Xoff && ~MIPoff
        ax(2) = subplot(1,3,2);
    elseif Xoff || MIPoff
        ax(2) = subplot(1,2,2);
    elseif Xoff && MIPoff
        ax(2) = subplot(1,1,1);
    end
    hold(ax(2),'on')
end
% Add example image of tuning map
if ~MIPoff
    if ~Yoff && ~Xoff
        ax(3) = subplot(1,3,3);
    elseif Yoff || Xoff
        ax(3) = subplot(1,2,2);
    elseif Yoff && Xoff
        ax(3) = subplot(1,1,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% go through all superclass objects and extract pixel positions and tunings

% cross-recording pooled vectors:
posX = [];
posY = [];
angD = [];
selW = [];
magW = [];
indAllFly = [];

rhotemp = nan;

% individual x- or y- only variables:
phi_0X = [];
amaxX = [];
rhoX = [];
indXFly = [];
inclX = [];
pvalX = nan;
posXindiv = {};
angXindiv = {};
selXindiv = {};

phi_0Y = [];
amaxY = [];
rhoY = [];
indYFly = [];
inclY = [];
pvalY = nan;
angYindiv = {};
posYindiv = {};
selYindiv = {};

for oidx = 1:length(obj)
    % First check that object included a pol tuning experiment
    if ~any(obj(oidx).Layers(1).TrialSeqNum==4) && ~any(obj(oidx).Layers(1).TrialSeqNum==8)
        disp('No pol exp (4,8) available')
        continue
    end
    
    % check that there are some pixels above the designated pol threshold
    if ~any(obj(oidx).MIP.polSelImg(:)>obj(oidx).MIP.polSelThreshold)
        continue
    end
    %useMIP
    
    if contains(obj(oidx).Name,' L ') % recording made on fly's left side
        leftFlag = 1;
    else
        leftFlag = 0;
    end
    
    % We might use the MIP obj for finding image rotation - if it's missing
    % because only a single layer was recorded, use that layer
    clearMIPflag = 0;
    if isempty(obj(oidx).MIP) && length(length(obj(oidx).Layers)) == 1
        obj(oidx).MIP = obj(oidx).Layers(1);
        clearMIPflag = 1;
    end
    
    % Fetch all points in object's MIP (no thresholding):
    % temporarily set thrshold to zero
    pst = obj(oidx).MIP.polSelThreshold;
    obj(oidx).MIP.polSelThreshold = 0 ;
    % Get x and y position of each pixel in the MIP mask
    % (image(x,y) are swapped here to correct for 90deg rotation of 2p
    % images: now ventral-dorsal is y(min:max)):
    [i,j] = ind2sub(size(obj(oidx).MIP.polPix),find(~isnan(obj(oidx).MIP.polPix)));
    % reset threshold
    obj(oidx).MIP.polSelThreshold = pst;
    
    % Get rotation angle for this experiment to align the major axis of the
    % neuropil with either x- (medial->lateral on rhs) or y- (ventral->dorsal)
    % (the 90deg rotation of all images from the 2p is taken care of
    % separately above - don't consider here)
    if ~rotateOFF
        
%         if strcmp(obj(oidx).Area,'Me')
%             % For recordings across a layer of the medulla the points fall
%             % approximately in a line when viewed dorsally: fit a line to
%             % all points in the MIP image
%             p = polyfit(i,j,1);
%             rotateAngle = -atand(p(1)) ;
%             
%         else
if strcmp(obj(oidx).Area,'AOTu')
            
            rotateAngle = 20;

            if leftFlag
                rotateAngle = -rotateAngle;
            end
        else
            rotateAngle = 0;
        end
    else
        rotateAngle = 0;
    end
    
    %     % Also modify the scale appropriately
    %     ImicronPix = obj(oidx).MIP.micronsPerPixel*(max(j) - min(j) )*cosd(rotateAngle) + obj(oidx).MIP.micronsPerPixel*(max(i) - min(i) )*sind(rotateAngle);
    %     JmicronPix = -obj(oidx).MIP.micronsPerPixel*(max(j) - min(j) )*sind(rotateAngle) + obj(oidx).MIP.micronsPerPixel*(max(i) - min(i) )*cosd(rotateAngle);
    I = []; J = []; D = []; sel = []; mag=[];
    
    loadLayerMasks(obj(oidx).Layers)
    
    if ~isempty(obj(oidx).Layers(1).layerMask) && ~isempty(obj(oidx).Layers(end).layerMask) && UseMSPifLayerMasksEmpty~=2
        
        % Run through all layers in superobj
        for Lidx = 1:length(obj(oidx).Layers)
            
            if isempty(obj(oidx).Layers(Lidx).layerMask)
                disp(['obj(' num2str(oidx) ').Layer ' num2str(Lidx) ' has no layerMask! Whole frame is included, which can distort results'])
            end
            
            % Rotate image and extract data
            [layerI,layerJ,layerD,layersel,layermag,imrotateAngle]=getPixRotate(obj(oidx).Layers(Lidx),rotateAngle,rotateOFF);
            
            D = [D;layerD];
            sel = [sel;layersel];
            mag = [mag;layermag];
            
            I = [I;layerI];
            J = [J;layerJ];
            
        end
        
    elseif ~isempty(obj(oidx).MIP) && UseMSPifLayerMasksEmpty>0
        
        % Use a projection (MIP or MSP, depeding on obj.MIP.UseMSP setting)
        loadLayerMasks(obj(oidx).MIP)
        
        if isempty(obj(oidx).MIP.layerMask)
            disp(['obj(' num2str(oidx) ').MIP has no layerMask! Whole frame is included, which can distort results'])
        elseif UseMSPifLayerMasksEmpty == 1
            disp(['obj(' num2str(oidx) ') Layers have no Masks - Using MSP instead'])
        end
        
        % Rotate image and extract data
        [layerI,layerJ,layerD,layersel,layermag,imrotateAngle]=getPixRotate(obj(oidx).MIP,rotateAngle,rotateOFF);
        
        D = [D;layerD];
        sel = [sel;layersel];
        mag = [mag;layermag];
        
        I = [I;layerI];
        J = [J;layerJ];
        
    else
        continue
    end
    
    
    if ~strcmp(obj(oidx).Area,'PB')
        
        % Scale to within range [0:1]
        X = (I - min(I))./(max(I) - min(I));
        Y = (J - min(J))./(max(J) - min(J));
        
    else

        % special case for E-PG recordings in the bridge:
        % map x,y pixel coordinates onto a curve (from left to right across
        % the bridge) and use distance along the curve as position x (and
        % ignore y)
        X = getCurvilinearDist(obj(oidx).MIP,I,J);
 
        % Y remains the same, but is no longer meaningful in relation to X.
        % Do not use
        Y = (J - min(J))./(max(J) - min(J));
        leftFlag = 0;
        
        % To check against regular projection of pixel position onto horizontal axis
        %{
        figure
        
        subplot(1,2,1)
        cm = flipud(hsv(100));
        cm = [cm; cm];
        Yup = round(( (I - min(I))./(max(I) - min(I))).*199);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
        scatter(I,J,100,C,'filled')
        axis image
        title('horizontal position')
        
        subplot(1,2,2)
        cm = flipud(hsv(100));
        cm = [cm; cm];
        Yup = round(X.*199);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
        scatter(I,J,100,C,'filled')
        axis image
         title('position along curve')
       %}
        

    end

    if strcmp(obj(oidx).Area,'AOTu') && leftFlag
        % Reverse X Y directions to pool with r.h.s. data
        X = 1-X;
        Y = 1-Y;
        leftFlag = 0;
    elseif strcmp(obj(oidx).Area,'Bu') && leftFlag
        X = 1-X;
        Y = 1-Y;
        leftFlag = 0;
    end
    
    mag(isnan(mag)) = 0;
    sel(isnan(sel)) = 0;
    
    
    if clearMIPflag
        obj(oidx).MIP = [];
    end
    
    incl = 0;
    % Plot individual x- fit line
    if ~Xoff
        if applyWeights
            selIndiv = sel;
        else
            selIndiv = ones(size(sel));
        end
        [rhotemp, pval, phi_0, amax] = runCorr_AntPost(obj(oidx),D,X,selIndiv,convertAngles4Pub);
        
        if pval < 0.05 || ~limitToSignifRecs
            
            if plotIndividFits
                % if initial fit appears significant, plot fit line
                [lineH] = plotPhi_AntPost(obj(oidx),amax,linspace(0,1,100),phi_0,ax(1),convertAngles4Pub);
                
                
                for n = 1:length(lineH)
                    lineH(n).Color = greyFitLines;
                    lineH(n).LineWidth = 0.25;
                end
            end
            
        end
        if ( limitToSignifRecs && pval < 0.05) || ~limitToSignifRecs
            % if option is enabled: add data to pool. Otherwise, skip.
            % If option disabled, all signif + n.s. data is retained
            
            posXindiv{end+1} = X;
            angXindiv{end+1} = D;
            selXindiv{end+1} = selIndiv;
             
            inclX = [inclX;oidx];
            
            phi_0X(end+1) = phi_0;
            amaxX(end+1) = amax;
            rhoX(end+1) = rhotemp;
            indXFly(end+1) = obj(oidx).Fly;
            incl = 1;
        end
    end
    
    
    % Plot individual y- fit line
    if ~Yoff
        if applyWeights
            selIndiv = sel;
        else
            selIndiv = ones(size(sel));
        end
        [rhotemp, pval, phi_0, amax] = runCorr_VentDors(obj(oidx),D,Y,selIndiv,convertAngles4Pub);
        
        if pval < 0.05 || ~limitToSignifRecs
            
            if plotIndividFits
                % if initial fit appears significant, plot fit line
                [lineH] = plotPhi_VentDors(obj(oidx),amax,linspace(0,1,100),phi_0,ax(2),convertAngles4Pub);
                
                for n = 1:length(lineH)
                    lineH(n).Color = greyFitLines;
                    lineH(n).LineWidth = 0.25;
                end
            end
        end
        
        if ( limitToSignifRecs && pval < 0.05) || ~limitToSignifRecs
            % if option is enabled: add data to pool. Otherwise, skip.
            % If option disabled, all signif + n.s. data is retained
            
            posYindiv{end+1} = Y;
            angYindiv{end+1} = D;
            selYindiv{end+1} = selIndiv;
           
            inclY = [inclY;oidx];
            
            phi_0Y(end+1) = phi_0;
            amaxY(end+1) = amax;
            rhoY(end+1) = rhotemp;
            indYFly(end+1) = obj(oidx).Fly;
            incl = 1;
        end
        
    end
    
    % Add to cross-object vectors:
    posX = [posX;X];
    posY = [posY;Y];
    angD = [angD;D];
    selW = [selW;sel];
    magW = [magW;mag];
    indAllFly = [indAllFly;obj(oidx).Fly];
    
    if incl
        inclIdx = oidx; % we might need common info from one of the objects or to plot an example tuning image
    end
end

if isempty(angD) 
    close gcf
    return
end
if length(indAllFly)==1
    inclIdx = oidx;
end

% Weightings (if used) will be pol selectivity
W = selW;
[W,sortIdx] = sort(W,'ascend'); % in order to plot darker scatter dots over lighter ones
posX = posX(sortIdx);
angD = angD(sortIdx);
posY = posY(sortIdx);
origSelW = W;
if ~applyWeights
    W = ones(size(W)); % set all weights to one
end

origAngD = angD;
if convertAngles4Pub
    angD = rec2weir(origAngD);
end

% Display rotated image:
if ~MIPoff
    if isempty(obj(inclIdx).MIP) && length(length(obj(inclIdx).Layers)) == 1
        obj(inclIdx).MIP = obj(inclIdx).Layers(1);
        clearMIPflag = 1;
    end
    pst = obj(inclIdx).MIP.polSelThreshold;
    obj(inclIdx).MIP.polSelThreshold = 0 ;
    
    
    rotK = imrotate(obj(inclIdx).MIP.polPix,imrotateAngle);
    obj(inclIdx).MIP.polSelThreshold = pst ;
    
    imagesc(ax(3),rotK);
    if applyTuningCols
        colormap([[1 1 1];flipud(hsv)])
    else
        colormap([[1 1 1];flipud(hsv)])
    end
    set(ax(3),'Color',colstr)
    
    axis image
    [Kj,Ki] = ind2sub(size(rotK),find(~isnan(rotK)&(rotK~=0)));
    
    ylim([min(Kj),max(Kj)])
    xlim([min(Ki),max(Ki)])
    
    axis off
    if clearMIPflag
        obj(inclIdx).MIP = [];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot x- correlation

if ~Xoff
    
    %%% Create scatter plot:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if applyTuningCols % Colormap for each point based on tuning
        cm = flipud(hsv(18001));
        Yup = round(origAngD*100);
        C = cm(Yup+1,:);
    elseif applyPSICols % Colormap for each point based on selectivty
        if darkMode
            cm = [0 0 0;(magma(255))];
        else
            cm = [1 1 1;flipud(magma(255))];
        end
        Yup = round(origSelW.*255);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
    else % No colormap
        cm = ones(18001,3).*ROIcolor(1);
        Yup = round(origAngD*100);
        C = cm(Yup+1,:);
    end
    
    % Decimate arrays before plotting to limit number of pts in plot
    if applyPtsLimit && length(posX) > ptsLimit
        f = floor(length(posX)/ptsLimit);
        posXdec = posX(1:f:end);
        angDdec = angD(1:f:end);
        Wdec = origSelW(1:f:end);
        if ~isempty(C)
            Cdec = C(1:f:end,:);
        end
        sHX = scatter(ax(1),posXdec,angDdec,ptSize,Cdec,'filled');
    else
        sHX = scatter(ax(1),posX,angD,ptSize,C,'filled');
    end
    if darkMode
        set(gca,'color',[1 1 1])
    end
    uistack(sHX,'bottom')
    
    
    %%% Find correlation and plot fit line:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  ~isempty(indXFly)
        Xflies = unique(indXFly);
        flymean_phi_0Y = [];
        flymean_rhoX = [];
        for n = 1:length(Xflies)
            [~, flymean_phi_0X(n)] = circ_axialmean(phi_0X(indXFly==Xflies(n)),1,2);
            flymean_rhoX(n) = tanh(mean(atanh(abs(rhoX(indXFly==Xflies(n))))));
        end
        rhoIndivMeanX = tanh(mean(atanh(flymean_rhoX)));
        rhoIndivSEMX = tanh(std(atanh(flymean_rhoX)/sqrt(length(flymean_rhoX))));
        
        if plotIndividFits
            [~, mean_phi_0X] = circ_axialmean(flymean_phi_0X,1,2);
            [lineHX2, phi_fitX] = plotPhi_AntPost(obj(inclIdx),mean(amaxX),posX,mean_phi_0X,ax(1),convertAngles4Pub);
            
            for n = 1:length(lineHX2)
                lineHX2(n).Color = meanFitLine;
                lineHX2(n).LineWidth = 2;
            end
            
        end
    end
    
    [rhoPooledX, pvalPooledX, phi_0, amax] = runCorr_AntPost(obj(inclIdx),origAngD,posX,W,convertAngles4Pub);
    if plotPooledFit
        [lineHX, phi_fitX] = plotPhi_AntPost(obj(inclIdx),amax,posX,phi_0,ax(1),convertAngles4Pub);
        for n = 1:length(lineHX)
            lineHX(n).Color = 'k';
            lineHX(n).LineWidth = 1;
            lineHX(n).LineStyle = ':';
        end
    end
    
    
    corrDirection = 'x';
    
    dataIndiv = {};
    dataIndiv(1,:) = posXindiv;
    dataIndiv(2,:) = angXindiv;
    dataIndiv(3,:) = selXindiv;
    
    flyIdx = indXFly;
    rhoIndiv = rhoX;
    
    switch corrStatMethodPooled
        case 0
            % use existing values
        case 1
            % bootstrap circlin corr distribution from all paired pos/tuning data
            useAll = 1;
            [rhoPooledX,CIpooledX] = boot_Corr( obj(inclIdx), corrDirection , useAll, dataIndiv, flyIdx,rhoIndiv,showPlot );
            
        case 2
            % bootstrap circlin corr distribution from subset of paired pos/tuning data - average rec length
            useAll = 0;
            [rhoPooledX,CIpooledX] = boot_Corr( obj(inclIdx), corrDirection , useAll, dataIndiv, flyIdx, rhoIndiv, showPlot );
    end
    
    switch corrStatMethodIndiv
        case 3
            % bootstrap mean corclincorr coef from the sample of
            % coefs found for individual recordings
            useFisherZ = 0;
            [rhoIndivMeanX,CIindivX] = boot_meanCorr( obj(inclIdx), corrDirection , useFisherZ, dataIndiv, flyIdx, rhoIndiv, showPlot );
            
        case 4
            % bootstrap mean corclincorr coef from the sample of
            % coef values found for individual recordings
            useFisherZ = 1;
            [rhoIndivMeanX,CIindivX] = boot_meanCorr( obj(inclIdx), corrDirection , useFisherZ, dataIndiv, flyIdx, rhoIndiv, showPlot );
            
        case 5
            % use existing values
    end
    
    
    fprintf('\n\n')
    disp(['[' obj(inclIdx).Line '_' obj(inclIdx).Area '] ant-post (X):'])
    fprintf('\n')
    disp('mean of individual fly values')
    indivNfliesX = length(unique(indXFly));
    fprintf('N = %d flies\n',indivNfliesX)
    fprintf('correlation coefficient = %2.2f\n',rhoIndivMeanX)
    if exist('CIindivX')
        fprintf('CI95= [%2.2f, %2.2f]\n',CIindivX(1),CIindivX(2))
    elseif exist('rhoIndivSEMX')
        fprintf('SEM = %2.2f\n',rhoIndivSEMX)
    end
    
    fprintf('\n')
    pooledNfliesX = length(unique(indAllFly));
    disp('all recordings pooled values')
    fprintf('n = %d recordings\n',pooledNfliesX)
    fprintf('correlation coefficient = %2.2f\n',rhoPooledX)
    if exist('CIpooledX')
        fprintf('CI95boot = [%2.2f, %2.2f]\n',CIpooledX(1),CIpooledX(2))
    end
    
    if runPermutationTest
        [~, pvalPooledX] = runCorr_AntPost(obj(inclIdx),origAngD,posX,W,showPlot,1); % second to last arg is showPlot (not pubAngs)
    else
        % use pooled val found above
    end
    fprintf('P-value = %2.8f\n',pvalPooledX)
    fprintf('\n\n')
    
    %  set(ax(1),'xlim',[0 1])
    %  set(ax(1),'xtick',[0 1])
    if convertAngles4Pub
        set(ax(1),'ylim',[-90 90])
        set(ax(1),'ytick',[-90:30:90])
        set(ax(1),'yticklabel',{['-90' char(176)];'';'';['0' char(176)];'';'';['90' char(176)]})
    else
        set(ax(1),'ylim',[0 180])
        set(ax(1),'ytick',[0:30:180])
        set(ax(1),'yticklabel',{['0' char(176)];'';'';['90' char(176)];'';'';['180' char(176)]})
    end
    set(ax(1),'Color',colstr)
    ylabel(ax(1),'preferred AoP')
    % offsetAxes(ax(1))
   
    % old verions printed stats in ax title
    %{
    if exist('CIpooledX') && exist('CIindivX')
        title(ax(1),sprintf('individual: r=%1.2f, CI[%2.2f,%2.2f], N=%d\npooled: r=%1.2f, CI[%2.2f,%2.2f], p=%1.5f, N=%d',...
            rhoIndivMeanX,CIindivX(1),CIindivX(2),indivNfliesX,rhoPooledX,CIpooledX(1),CIpooledX(2),pvalPooledX,pooledNfliesX))
    elseif exist('CIindivX')
        title(ax(1),sprintf('individual: r=%1.2f, CI[%2.2f,%2.2f], N=%d\npooled: r=%1.2f, p=%1.5f, N=%d',...
            rhoIndivMeanX,CIindivX(1),CIindivX(2),indivNfliesX,rhoPooledX,pvalPooledX,pooledNfliesX))
    elseif exist('rhoIndivSEMX')
          title(ax(1),sprintf('individual: r=%1.2f, SEM%2.2f, N=%d\npooled: r=%1.2f, p=%1.5f, N=%d',...
            rhoIndivMeanX,rhoIndivSEMX,indivNfliesX,rhoPooledX,pvalPooledX,pooledNfliesX))
    end
    %}
    
    set(ax(1),'xlim',[0 1])
    set(ax(1),'xtick',[0,0.25,0.5,0.75 1])
    xlabel(ax(1),'normalized position')
    if strcmp(obj(inclIdx).Area,'Me')
        set(ax(1),'xticklabel',{'anterior';'';'';'';'posterior'})
    else
        set(ax(1),'xticklabel',{'medial';'';'';'';'lateral'})
    end
    
    pbaspect(ax(1),[1,1,1])
    ax(1).XAxis.Color = colstr;
    ax(1).YAxis.Color = colstr;
    ax(1).Color = acolstr;
    % ax(1).FontWeight = 'bold';
    if strcmp(obj(inclIdx).Cell,'R7R8')
        text(1,80,' R8','FontSize',axisLabelFontSize)
        text(1,-10,' R7','FontSize',axisLabelFontSize)
    end
    
    if MIPoff
        setAVPaxes(ax(1),2)
    else
        ax(1).FontWeight = 'bold';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot y- correlation

if ~Yoff
    
    %%% Create scatter plot:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if applyTuningCols % Colormap for each point based on tuning
        cm = flipud(hsv(18001));
        Yup = round(origAngD*100);
        C = cm(Yup+1,:);
    elseif applyPSICols % Colormap for each point based on selectivty
        if darkMode
            cm = [0 0 0;(magma(255))];
        else
            cm = [1 1 1;flipud(magma(255))];
        end
        Yup = round(origSelW.*255);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
    else % No colormap
        cm = ones(18001,3).*ROIcolor(1);
        Yup = round(origAngD*100);
        C = cm(Yup+1,:);
    end
    % Decimate arrays before plotting: limit to 1000 pts
    if applyPtsLimit && length(posY) > ptsLimit
        f = floor(length(posY)/ptsLimit);
        posYdec = posY(1:f:end);
        angDdec = angD(1:f:end);
        Wdec = origSelW(1:f:end);
        if ~isempty(C)
            Cdec = C(1:f:end,:);
        end
        sHY = scatter(ax(2),posYdec,angDdec,ptSize,Cdec,'filled');
    else
        sHY = scatter(ax(2),posY,angD,ptSize,C,'filled');
    end
    uistack(sHY,'bottom')
    
    %%% Find correlation and plot fit line:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isempty(indYFly)
        
        Yflies = unique(indYFly);
        flymean_phi_0Y = [];
        flymean_rhoY = [];
        for n = 1:length(Yflies)
            [~, flymean_phi_0Y(n)] = circ_axialmean(phi_0Y(indYFly==Yflies(n)),1,2);
            flymean_rhoY(n) = tanh(mean(atanh((rhoY(indYFly==Yflies(n))))));
        end
        rhoIndivMeanY = tanh(mean(atanh(flymean_rhoY)));
        rhoIndivSEMY = tanh(std(atanh(flymean_rhoY)/sqrt(length(flymean_rhoY))));
        
        if plotIndividFits
            [~, mean_phi_0Y] = circ_axialmean(flymean_phi_0Y,1,2);
            [lineHY2, phi_fitY] = plotPhi_VentDors(obj(inclIdx),mean(amaxY),posY,mean_phi_0Y,ax(2),convertAngles4Pub);
            for n = 1:length(lineHY2)
                lineHY2(n).Color = meanFitLine;
                lineHY2(n).LineWidth = 2;
            end
        end
    
    
    [rhoPooledY, pvalPooledY, phi_0, amax] = runCorr_VentDors(obj(inclIdx),origAngD,posY,W,convertAngles4Pub);
    if plotPooledFit
        [lineHY, phi_fitY] = plotPhi_VentDors(obj(inclIdx),amax,posY,phi_0,ax(2),convertAngles4Pub);
        for n = 1:length(lineHY)
            lineHY(n).Color = 'k';
            lineHY(n).LineWidth = 1;
            lineHY(n).LineStyle = ':';
        end
    end
    
    corrDirection = 'y';
    
    dataIndiv = {};
    dataIndiv(1,:) = posYindiv;
    dataIndiv(2,:) = angYindiv;
    dataIndiv(3,:) = selYindiv;
    
    
    flyIdx = indYFly;
    rhoIndiv = rhoY;
    
    switch corrStatMethodPooled
        case 0 
            % use existing values
        case 1
            % bootstrap circlin corr distribution (from all paired pos/tuning data)
            useAll = 1;
            [rhoPooledY,CIpooledY] = boot_Corr( obj(inclIdx), corrDirection , useAll, dataIndiv, flyIdx,rhoIndiv,showPlot );
            
        case 2
            % bootstrap circlin corr distribution (from subset of paired pos/tuning data - average rec length)
            useAll = 0;
            [rhoPooledY,CIpooledY] = boot_Corr( obj(inclIdx), corrDirection , useAll, dataIndiv, flyIdx, rhoIndiv, showPlot );
    end
    switch corrStatMethodIndiv
        case 3
            % bootstrap mean corclincorr coef from the sample of
            % coefs found for individual recordings
            useFisherZ = 0;
            [rhoIndivMeanY,CIindivY] = boot_meanCorr( obj(inclIdx), corrDirection , useFisherZ, dataIndiv, flyIdx, rhoIndiv, showPlot );
            
        case 4
            % bootstrap mean corclincorr coef from the sample of
            % coef values found for individual recordings
            useFisherZ = 1;
            [rhoIndivMeanY,CIindivY] = boot_meanCorr( obj(inclIdx), corrDirection , useFisherZ, dataIndiv, flyIdx, rhoIndiv, showPlot );
            
        case 5
            % use existing values
    end
    
    if runPermutationTest
        [~, pvalPooledY] = runCorr_VentDors(obj(inclIdx),origAngD,posY,W,showPlot,1); % second to last arg is showPlot (not pubAngs)
    else
        % use pooled val found above
    end
    
    fprintf('\n\n')
    disp(['[' obj(inclIdx).Line '_' obj(inclIdx).Area '] vent-dors (Y):'])
    fprintf('\n')
    
    disp('mean of individual fly values')
    indivNfliesY = length(unique(indYFly));
    fprintf('N = %d flies\n',indivNfliesY)
    fprintf('correlation coefficient = %2.2f\n',rhoIndivMeanY)
    if exist('CIindivY')
        fprintf('CI95= [%2.2f, %2.2f]\n',CIindivY(1),CIindivY(2))
    elseif exist('rhoIndivSEMY')
        fprintf('SEM = %2.2f\n',rhoIndivSEMY)
    end
    
    fprintf('\n')
    pooledNfliesY = length(unique(indAllFly));
    disp('all recordings pooled values')
    fprintf('N = %d recordings\n',pooledNfliesY)
    fprintf('correlation coefficient = %2.2f\n',rhoPooledY)
    if exist('CIpooledY')
        fprintf('CIboot95= [%2.2f, %2.2f]\n',CIpooledY(1),CIpooledY(2))
    end
    

    fprintf('P-value = %2.8f\n',pvalPooledY)
    fprintf('\n\n')
    end
    %% y plot
    if convertAngles4Pub
        set(ax(2),'ylim',[-90 90])
        set(ax(2),'ytick',[-90:30:90])
        set(ax(2),'yticklabel',{['-90' char(176)];'';'';['0' char(176)];'';'';['90' char(176)]})
    else
        set(ax(2),'ylim',[0 180])
        set(ax(2),'ytick',[0:30:180])
        set(ax(2),'yticklabel',{['0' char(176)];'';'';['90' char(176)];'';'';['180' char(176)]})
    end
    ax(2).XAxis.Color = colstr;
    ax(2).YAxis.Color = colstr;
    ax(2).Color = acolstr;
    %     xlabel(ax(2),{['ventral '  char(8594) ' dorsal'];'position'})
    set(ax(2),'xlim',[0 1])
    set(ax(2),'xtick',[0,0.25,0.5,0.75, 1])
    set(ax(2),'xticklabel',{'ventral';'';'';'';'dorsal'})
    
    ax(2).XLabel.String =   'position';
    %     ax(2).XLabel.String =  ['ventral '  char(8594) ' dorsal'];
    
    ylabel(ax(2),'preferred AoP')
    
    % old verions printed stats in ax title
    %{
    if exist('CIpooledY') && exist('CIindivY')
        title(ax(2),sprintf('individual: r=%1.2f, CI[%2.2f,%2.2f], N=%d\npooled: r=%1.2f, CI[%2.2f,%2.2f], p=%1.5f, N=%d',...
            rhoIndivMeanY,CIindivY(1),CIindivY(2),indivNfliesY,rhoPooledY,CIpooledY(1),CIpooledY(2),pvalPooledY,pooledNfliesY))
    elseif exist('CIindivY')
        title(ax(2),sprintf('individual: r=%1.2f, CI[%2.2f,%2.2f], N=%d\npooled: r=%1.2f, p=%1.5f, N=%d',...
            rhoIndivMeanY,CIindivY(1),CIindivY(2),indivNfliesY,rhoPooledY,pvalPooledY,pooledNfliesY))
    elseif exist('rhoIndivSEMY')
          title(ax(2),sprintf('individual: r=%1.2f, SEM%2.2f, N=%d\npooled: r=%1.2f, p=%1.5f, N=%d',...
            rhoIndivMeanY,rhoIndivSEMY,indivNfliesY,rhoPooledY,pvalPooledY,pooledNfliesY))
    end
    %}
    pbaspect(ax(2),[1,1,1])
    
    if MIPoff
        
        setAVPaxes(ax(2),2)
    else
        ax(2).FontWeight = 'bold';
        
    end
end

if MIPoff
    tightfig(gcf);
end

end

function [Ii,Jj,D,sel,mag,imrotateAngle]=getPixRotate(obj,rotateAngle,rotateOFF)

% Find the (y,x) position of each polselective pixel in tuning map 'k'
[i,j] = ind2sub(size(obj.polPix),find(~isnan(obj.polPix)));

% Option to rotate all pixels first:
% Auto-rotation uses a linear fitted line through all points to find the
% orientation of the cells in the image, then rotates to align with x-axis.
% 90 degree rotation to correct for 2p orientation might be sufficient
% here.
if ~rotateOFF
    
    
    % We don't simply rotate the image becuase that method interpolates and
    % creates new data points
    Jj = [j*cosd(rotateAngle) + i*sind(rotateAngle)];
    Ii = [-j*sind(rotateAngle) + i*cosd(rotateAngle)];
    
    %
    %         % Also modify the scale appropriately
    %         ImicronPix = obj(oidx).Layers(Lidx).micronsPerPixel*(max(i) - min(i) + 1)*cosd(rotateAngle) + obj(oidx).Layers(Lidx).micronsPerPixel*(max(j) - min(j) + 1)*sind(rotateAngle);
    %         JmicronPix = -obj(oidx).Layers(Lidx).micronsPerPixel*(max(i) - min(i) + 1)*sind(rotateAngle) + obj(oidx).Layers(Lidx).micronsPerPixel*(max(j) - min(j) + 1)*cosd(rotateAngle);
    
    imrotateAngle = rotateAngle + 90;
else
    imrotateAngle = 0;
    Ii = i;
    Jj = j;
    %         ImicronPix = obj(oidx).Layers(Lidx).micronsPerPixel*(max(j) - min(j) + 1);
    %         JmicronPix = obj(oidx).Layers(Lidx).micronsPerPixel*(max(i) - min(i) + 1);
end

D = obj.polPix(find(~isnan(obj.polPix)));
sel = obj.polSelImg(find(~isnan(obj.polPix)));
mag = obj.fftMagImg(find(~isnan(obj.polPix)));
end
function [rho, pval, phi_0, amax] = runCorr_AntPost(obj,angularData,XpositionData,weights,convertAngles4Pub,permTest)
if contains(obj.Name,' L ') % recording made on fly's left side
    leftFlag = 1;
else
    leftFlag = 0;
end

% Modify for particular brain areas
if strcmp(obj.Cell,'R7R8')
%     if leftFlag
%         bounds = [0 2];
%     else
        bounds = [-4 0.2];
%     end
    if exist('permTest') && permTest
        showPlot = convertAngles4Pub; % alternative use for argument.. will convert angles by default
        [rho, pval ] = circlin_wcorr_permute(deg2rad(4*angularData), XpositionData,bounds, weights,showPlot);
        phi_0 = [];
        amax = [];
    else
        [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(4*angularData), XpositionData,bounds, weights);
    end
elseif strcmp(obj.Area,'PB')

        bounds = [-2 2];
        posX2 = 2*XpositionData;
        
        if exist('permTest') && permTest
            showPlot = convertAngles4Pub; % alternative use for argument.. will convert angles by default
            [rho, pval ] = circlin_wcorr_permute(deg2rad(2*angularData),  mod(posX2,1),bounds, weights,showPlot);
            phi_0 = [];
            amax = [];
        else
            [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*angularData), mod(posX2,1),bounds, weights);
        end

else
    if leftFlag
        bounds = [-2 0.2]; % changed from [0 1]
    else
        bounds = [-2 0.2];% changed from [-1 0]
    end
    if exist('permTest') && permTest
        showPlot = convertAngles4Pub; % alternative use for argument.. will convert angles by default
        [rho, pval] = circlin_wcorr_permute(deg2rad(2*angularData), XpositionData,bounds, weights,showPlot);
        phi_0 = [];
        amax = [];
    else
        [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*angularData), XpositionData,bounds, weights);
    end
end
if convertAngles4Pub
    rho = -rho;
end
end
function [lineHandle, phi_fit_returned] =  plotPhi_AntPost(obj,amax,XpositionData,phi_0,ax,convertAngles4Pub)
colstr = 'k';
acolstr = 'w';
if convertAngles4Pub
    lowerLim = -90;
    upperLim=90;
else
    lowerLim = 0;
    upperLim = 180;
end
if nargin < 5 || isempty(ax),ax=gca;end
% Modify for particular brain areas
n=0;
if strcmp(obj.Cell,'R7R8')
    
    phi_fit_returned =180*(amax/2)*XpositionData + rad2deg(phi_0/4);
    if convertAngles4Pub
        phi_fit_plot = -180*(amax/2)*XpositionData + rec2weir(rad2deg(phi_0/4));%;
    else
        phi_fit_plot  = phi_fit_returned;
    end
    hold(ax,'on')
    for phi_shift = [-360:90:360]
        
        this_phi_fit = phi_fit_plot + phi_shift;
        if mod(phi_shift,180) == 90
            h = plot(ax,XpositionData(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),this_phi_fit(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),'--','color',colstr);
        else
            h = plot(ax,XpositionData(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),this_phi_fit(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),colstr);
        end
        if ~isempty(h)
            n=n+1;
            lineHandle(n) = h;
        end
    end
elseif strcmp(obj.Area,'PB')
    XpositionData = mod(2*XpositionData,1);
    phi_fit_returned = 180*amax*XpositionData + rad2deg(phi_0/2);
    
    hold(ax,'on')
    if convertAngles4Pub
        phi_fit_plot = - 180*amax*XpositionData + rec2weir(rad2deg(phi_0/2));%;
    else
        phi_fit_plot  = phi_fit_returned;
    end
    for phi_shift = [-2160:180:2160]
        this_phi_fit = phi_fit_plot + phi_shift;
        h = plot(ax,XpositionData(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),this_phi_fit(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),colstr);
        if ~isempty(h)
            n=n+1;
            lineHandle(n) = h;
        end
    end

    
else
    
    phi_fit_returned = 180*amax*XpositionData + rad2deg(phi_0/2);
    
    hold(ax,'on')
    if convertAngles4Pub
        phi_fit_plot = - 180*amax*XpositionData + rec2weir(rad2deg(phi_0/2));%;
    else
        phi_fit_plot  = phi_fit_returned;
    end
    for phi_shift = [-360:180:360]
        this_phi_fit = phi_fit_plot + phi_shift;
        h = plot(ax,XpositionData(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),this_phi_fit(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),colstr);
        if ~isempty(h)
            n=n+1;
            lineHandle(n) = h;
        end
    end
end
end
function [rho, pval, phi_0, amax] = runCorr_VentDors(obj,angularData,YpositionData,weights,convertAngles4Pub,permTest)
% Modify for particular brain areas
if contains(obj.Name,' L ') % recording made on fly's left side
    leftFlag = 1;
else
    leftFlag = 0;
end
if strcmp(obj.Cell,'R7R8') % || strcmp(obj(oidx).Cell,'MeTu')
%     if leftFlag
%         bounds = [0 2];
%     else
        bounds = [-4 0.2];
%     end
    if exist('permTest') && permTest
        showPlot = convertAngles4Pub; % alternative use for argument.. will convert angles by default
        [rho,pval] = circlin_wcorr_permute(deg2rad(4*angularData), YpositionData,bounds, weights,showPlot);
        phi_0 = [];
        amax = [];
    else
        [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(4*angularData), YpositionData,bounds, weights);
    end
    
else
%     if leftFlag
%         bounds = [-0.5 1];
%     else
        bounds = [-2 2]; 
%     end
    if exist('permTest') && permTest
        showPlot = convertAngles4Pub; % alternative use for argument.. will convert angles by default
        [rho,pval] = circlin_wcorr_permute(deg2rad(2*angularData), YpositionData,bounds, weights,showPlot);
        phi_0 = [];
        amax = [];
    else
        [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*angularData), YpositionData,bounds, weights);
  
    end
end

if convertAngles4Pub
    rho = -rho;
end
end
function [lineHandle, phi_fit_returned] = plotPhi_VentDors(obj,amax,YpositionData,phi_0,ax,convertAngles4Pub)
colstr = 'k';
acolstr = 'w';
if convertAngles4Pub
    lowerLim = -90;
    upperLim=90;
else
    lowerLim = 0;
    upperLim = 180;
end
if nargin < 5 || isempty(ax),ax=gca;end
% Modify for particular brain areas
n=0;
if strcmp(obj.Cell,'R7R8')%% || strcmp(obj.Cell,'MeTu')
    
    
    hold(ax,'on')
    phi_fit_returned = 180*(amax/2)*YpositionData + rad2deg(phi_0/4);
    
    %    phi_fit_plot  = phi_fit_returned;
    if convertAngles4Pub
        phi_fit_plot = -180*(amax/2)*YpositionData + rec2weir(rad2deg(phi_0/4));
    else
        phi_fit_plot  = phi_fit_returned;
    end
    
    for phi_shift = [-720:90:720]
        
        this_phi_fit = phi_fit_plot + phi_shift;
        if mod(phi_shift,180) == 90
            h = plot(ax,YpositionData(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),this_phi_fit(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),'--','color',colstr);
        else
            h = plot(ax,YpositionData(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),this_phi_fit(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),colstr);
        end
        
        if ~isempty(h)
            n=n+1;
            lineHandle(n) = h;
        end
    end
    
else
    phi_fit_returned = 180*amax*YpositionData + rad2deg(phi_0/2);
    
    if convertAngles4Pub
        phi_fit_plot = -180*amax*YpositionData + rec2weir(rad2deg(phi_0/2));%;
    else
        phi_fit_plot  = phi_fit_returned;
    end
    hold(ax,'on')
    for phi_shift =  [-720:180:720]
        
        this_phi_fit = phi_fit_plot  + phi_shift;
        h = plot(ax,YpositionData(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),this_phi_fit(this_phi_fit>=lowerLim & this_phi_fit<=upperLim),colstr);
        if ~isempty(h)
            n=n+1;
            lineHandle(n) = h;
        end
    end
    
    
end
end
function [rho,perm_pval] = circlin_wcorr_permute(phi, x, bounds, weights, showPlot)
% with help from https://courses.washington.edu/matlab1/Bootstrap_examples.html#1
nreps = 1000000;

P = phi;
X = x;
W = weights;

B = [-max(abs(bounds)) max(abs(bounds))];
perm = zeros(nreps,1);
rng default
fprintf('permutation loop')
for i=1:nreps
    if mod(i,1000)==0
        fprintf('.')
    end
    %shuffle the lsat scores and recalculate the correlation
    perm(i) = circlin_wcorr(P,X(randperm(length(X))),B,W);
    perm(i) = -perm(i); % to express in published angles (see rec2Weir)
end
fprintf('\n')

[rho,pval,~,~] = circlin_wcorr(P,X,bounds,W);
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

function [anglesPublished] = rec2weir(anglesRecorded)
anglesPublished = wrapTo180(-anglesRecorded-270);
end

function [ci] = bootbca(stat,bstat,jstat,alpha,weights,varargin)
% corrected and accelerated percentile bootstrap CI
% T.J. DiCiccio and B. Efron (1996), "Bootstrap Confidence Intervals",
% statistical science, 11(3)
% From MATLAB built-in function bootci.m, modified to provide bstat and
% jstat after manually resampling and running stat and jackknife

% bstat = bootstrp(nboot,bootfun,varargin{:},'weights',weights,'Options',bootstrpOptions);

% same as bootcper, this is the bias correction
z_0 = fz0(bstat,stat);

% % apply jackknife
% try
%     jstat = jackknife(bootfun,varargin{:},'Options',bootstrpOptions);
% catch ME
%     m = message('stats:bootci:JackknifeFailed',func2str(bootfun));
%     throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
% end
N = size(jstat,1);
if isempty(weights)
    weights = repmat(1/N,N,1);
else
    weights = weights(:);
    weights = weights/sum(weights);
end

% acceleration finding, see DiCiccio and Efron (1996)
mjstat = sum(bsxfun(@times,jstat,weights),1); % mean along 1st dim.
score = bsxfun(@minus,mjstat,jstat); % score function at stat; ignore (N-1) factor because it cancels out in the skew
iszer = all(score==0,1);
skew = sum(bsxfun(@times,score.^3,weights),1) ./ ...
    (sum(bsxfun(@times,score.^2,weights),1).^1.5) /sqrt(N); % skewness of the score function
skew(iszer) = 0;
acc = skew/6;  % acceleration

% transform back with bias corrected and acceleration
z_alpha1 = norminv(alpha/2);
z_alpha2 = -z_alpha1;
pct1 = 100*normcdf(z_0 +(z_0+z_alpha1)./(1-acc.*(z_0+z_alpha1)));
pct1(z_0==Inf) = 100;
pct1(z_0==-Inf) = 0;
pct2 = 100*normcdf(z_0 +(z_0+z_alpha2)./(1-acc.*(z_0+z_alpha2)));
pct2(z_0==Inf) = 100;
pct2(z_0==-Inf) = 0;

% inverse of ECDF
m = numel(stat);
lower = zeros(1,m);
upper = zeros(1,m);
for i=1:m
    lower(i) = prctile(bstat(:,i),pct2(i),1);
    upper(i) = prctile(bstat(:,i),pct1(i),1);
end

% return
ci = sort([lower;upper],1);
end % bootbca()
% -------------------------
function z0=fz0(bstat,stat)
% Compute bias-correction constant z0
% MATLAB built-in function bootci.m
z0 = norminv(mean(bsxfun(@lt,bstat,stat),1) + mean(bsxfun(@eq,bstat,stat),1)/2);
end   % fz0()

function [rhoMean,CI] = boot_Corr( obj, corrDirection , useAll, dataIndiv, flyIdx, rhoIndiv, showPlot )
if ~strcmp(corrDirection,'x') && ~strcmp(corrDirection,'y')
    disp('corrDireciton input to bootstrap func must be ''x'' or ''y''')
    return
end
uniqueFlies = unique(flyIdx);
posindiv = dataIndiv(1,:);
angindiv = dataIndiv(2,:);
selindiv =dataIndiv(3,:);

% find avg num of pixels in individual recordings
avgRecLength = ceil(mean(cellfun(@min,[(cellfun(@length,posindiv(~cellfun(@isempty,posindiv)),'UniformOutput',0))])));

% make pooled vectors from all individual recordings
allPos = []; allAng = []; allSel = [];
for nr =  1:length(posindiv)
    allPos = [allPos; posindiv{nr}()];
    allAng = [allAng; angindiv{nr}()];
    allSel = [allSel; selindiv{nr}()];
end

samplength = 10000;
bootRho = zeros(1,samplength);
rng default  % For reproducibility
fprintf('bootstrap loop')
for n = 1:10000
    if length(uniqueFlies)==1
        randFlies = uniqueFlies;
    else
        randFlies = datasample(uniqueFlies,length(uniqueFlies)); % sample flies, with replacement
    end
    
    % pool indices of all recordings from re-sampled flies (different size each time)
    selectedRecs = [];
    for m =1:length(randFlies)
        selectedRecs = [selectedRecs,find(flyIdx==randFlies(m))];
    end
    
    randRecs = posindiv(selectedRecs); % get data for those recordings
    randAngD = angindiv(selectedRecs); % and associated tunings
    randWeights = selindiv(selectedRecs); % and associated selectivites
    
    % assemble vectors with pooled data from resampled recordings
    rPos = []; rAng = []; rSel = [];
    for nr =  1:length(randRecs)
        rPos = [rPos; randRecs{nr}()];
        rAng = [rAng; randAngD{nr}()];
        rSel = [rSel; randWeights{nr}()];
    end
    
    % calculate our parameter of interest, the circ-linear corr coef:
    if useAll
        % sample size will be same as orginal pooled dataset (no. of pixels
        % from all recordings)
        sampSizeAll = length(allPos);
        [~,allResampIdx] = datasample(rPos,sampSizeAll);
        
        if strcmp(corrDirection,'x')
            [rho, ~, ~, ~] = runCorr_AntPost(obj,rAng(allResampIdx),rPos(allResampIdx),rSel(allResampIdx),1,0);
        elseif strcmp(corrDirection,'y')
            [rho, ~, ~, ~] = runCorr_VentDors(obj,rAng(allResampIdx),rPos(allResampIdx),rSel(allResampIdx),1,0);
        end
        
        bootRho(n) = rho;
        
    else
        % Sample size will be same as the average recording (no. of pixels
        % in average recording)
        sampSize = avgRecLength;
        % use a subset of data (Resampled with replacement)
        [rPosSub,rrResampIdx] = datasample(rPos,sampSize);
        rAngSub = rAng(rrResampIdx);
        rSelSub = rSel(rrResampIdx);
        
        if strcmp(corrDirection,'x')
            [rho, ~, ~, ~] = runCorr_AntPost(obj,rAngSub,rPosSub,rSelSub,1,0);
        elseif strcmp(corrDirection,'y')
            [rho, ~, ~, ~] = runCorr_VentDors(obj,rAngSub,rPosSub,rSelSub,1,0);
        end
        
        bootRho(n) = rho;
        
    end
    
    if mod(n,1000)==0
        fprintf('.')
    end
    
end
fprintf('\n')

% To find skew in bootstrapped coef distribution use jackknife func (see bootci.m)
% Also find the corr coef for the original pooled dataset (to plot on
% bootstrapped distribution histogram)
if strcmp(corrDirection,'x')
    if useAll
        jstat = jackknife(@runCorr_AntPost,obj,rAng(allResampIdx),rPos(allResampIdx),rSel(allResampIdx),1,0,'Options',statset('UseParallel',false));
    else
        jstat = jackknife(@runCorr_AntPost,obj,rAngSub,rPosSub,rSelSub,1,0,'Options',statset('UseParallel',false));
    end
    [pooledrho, ~, ~, ~] = runCorr_AntPost(obj,allAng,allPos,allSel,1,0);
elseif strcmp(corrDirection,'y')
    if useAll
        jstat = jackknife(@runCorr_VentDors,obj,rAng(allResampIdx),rPos(allResampIdx),rSel(allResampIdx),1,0,'Options',statset('UseParallel',false));
    else
        jstat = jackknife(@runCorr_VentDors,obj,rAngSub,rPosSub,rSelSub,1,0,'Options',statset('UseParallel',false));
    end
    [pooledrho, ~, ~, ~] = runCorr_VentDors(obj,allAng,allPos,allSel,1,0);
end
alpha = 0.05;
[CI] = bootbca(pooledrho,bootRho',jstat,alpha,[],[]);

if showPlot
    figure(919)
    subplot(2,1,2)
    
    histbins = [-1:0.05:1];
    
    h = histogram(bootRho,histbins);
    
    h.FaceColor = [0.7 0.7 0.7];
    h.EdgeAlpha = 0;
    
    hold on
    plot(CI(1)*[1,1],ylim,'r-','LineWidth',2);
    plot(CI(2)*[1,1],ylim,'r-','LineWidth',2);
    % Add individual recording coef values
    plot(rhoIndiv'.*ones(length(rhoIndiv),2),ylim,'k:','LineWidth',0.5);
    
    
    l1 = plot(mean(bootRho)*[1,1],ylim,'c-','LineWidth',2);
    l2 =  plot(pooledrho*[1,1],ylim,'g-','LineWidth',2);
    legend('boxoff',[l1 l2],'mean of bs coeffs','orig coeff','Location','northwest')
    
    
    if useAll
        title('all data')
    else
        title('subset')
    end
    
    xlim([-1 1])
    
    pbaspect([2 1 1])
end

rhohat = mean(jstat(:));
n = length(rPos);
biasrho = (n-1) * (pooledrho-rhohat);

rhoMean = pooledrho;

end

function [rhoMean,CI] = boot_meanCorr( obj, corrDirection , useFisherZ, dataIndiv, flyIdx, rhoIndiv, showPlot )

flies = unique(flyIdx);
flymean_rho = [];
for n = 1:length(flies)
    flymean_rho(n) = tanh(mean(atanh(rhoIndiv(flyIdx==flies(n)))));
end
meanFlyMeanRho =  tanh(mean(atanh(flymean_rho)));

% Bootstrap ciclin corr coef distribution from individual corr coef vals
samplength = 10000;
bootRho = zeros(1,samplength);
rng default  % For reproducibility
for n = 1:samplength
    if length(flies) == 1
        randFlies = flies;
    else
    randFlies = randsample(flies,length(flies),1); % sample flies, with replacement
    end
    r = rhoIndiv(ismember(flyIdx,unique(randFlies)));
    randRho = randsample(r,length(rhoIndiv),1);
    if useFisherZ
        bootRho(n) =  tanh(mean(atanh(randRho)));
    else
        bootRho(n) =  mean(randRho);
    end
end

if useFisherZ
    jstatfunc = @(rRhosamp) tanh(mean(atanh(rhoIndiv)));
    jstat = jackknife(jstatfunc,rhoIndiv);
    meanRho = tanh(mean(atanh(rhoIndiv)));
    
else
    jstat = jackknife(@mean,rhoIndiv,1);
    meanRho = (mean((rhoIndiv)));
end


alpha = 0.05;
[CI] = bootbca(meanRho,bootRho',jstat,alpha,[],[]);
rhoMean = meanFlyMeanRho;

if showPlot
    figure(919)
    subplot(2,1,2)
    histbins = [-1:0.05:1];
    h = histogram(bootRho,histbins);
    h.FaceColor = [0.7 0.7 0.7];
    h.EdgeAlpha = 0;
    hold on
    plot(CI(1)*[1,1],ylim,'r-','LineWidth',2);
    plot(CI(2)*[1,1],ylim,'r-','LineWidth',2);
    plot(rhoIndiv'.*ones(length(rhoIndiv),2),ylim,'k:','LineWidth',0.5);
    l1 =  plot(mean(bootRho)*[1,1],ylim,'c-','LineWidth',2);
    l2 =  plot(meanRho*[1,1],ylim,'y-','LineWidth',2);
    l3 = plot(meanFlyMeanRho*[1,1],ylim,'g-','LineWidth',2);
    if useFisherZ
        legend('boxoff',[l1 l2 l3],'mean of bs coeffs','mean of pooled z-transformed coefs','[mean of flymean coefs]','Location','northwest')
    else
        legend('boxoff',[l1 l2 l3],'mean of bs coeffs','mean of pooled coefs','mean of flymean coefs','Location','northwest')
    end
    title(sprintf('CI%d: %5.2f %5.2f',100*(1-alpha),CI(1),CI(2)));
    
    xlim([-1 1])
    pbaspect([2 1 1])
end
end