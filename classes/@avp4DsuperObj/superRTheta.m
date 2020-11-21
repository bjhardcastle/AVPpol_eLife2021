function superRTheta(obj,Yoff,Xoff,rotateOFF, darkMode)
% For AVP super objects: (superclass version of plotRThetacorr)
% Take pol selective pixels from the tuning map and plot polar Theta
% coords vs tuning angle. Find the circular-circular correlation and plot the
% line and R-value

applyTuningCols = 1; % Color scatter points according to pol tuning ( redundant info)
applyWeights = 0; % Weighted correlation using pol selectivity + normalized fft magnitude for each pixel
applyPtsLimit = 1;
ptsLimit = 1000;

UseMSPifLayerMasksEmpty = 2;
% 0, whole frame will be used if Layer layerMasks are empty.
% 1, we revert to MSP instead.
% 2, force MSP to be used

if nargin<2 || isempty(Yoff) || ~any(Yoff==0 || Yoff==1)
    Yoff = 1; % set to 1 to exclude the y-correlation plot
end
if nargin<3 || isempty(Xoff) || ~any(Xoff==0 || Xoff==1)
    Xoff = 0; % set to 1 to exclude the x-correlation plot
end
if nargin<4 || isempty(rotateOFF) || ~any(rotateOFF==0 || rotateOFF==1)
    rotateOFF = 0;
end
if nargin<5 || isempty(darkMode) || ~any(darkMode==0 || darkMode==1)
    darkMode = 0;
end
if darkMode
    Xoff = 1;
end


posX = [];
posY = [];
angD = [];
selW = [];
magW = [];

for oidx = 1:length(obj)
    % First check that object included a pol tuning experiment
    if ~any(obj(oidx).Layers(1).TrialSeqNum==4) && ~any(obj(oidx).Layers(1).TrialSeqNum==2)
        continue
    end
    if contains(obj(oidx).Name,' L') % recording made on fly's left side
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
    % neuropil with either x- (medial->lateral) or y- (ventral->dorsal)
    % (the 90deg rotation of all images from the 2p is taken care of
    % separately above - don't consider here)
    if ~rotateOFF
        
        if strcmp(obj(oidx).Area,'Me')
            % For recordings across a layer of the medulla the points fall
            % approximately in a line when viewed dorsally: fit a line to
            % all points in the MIP image
            p = polyfit(i,j,1);
            rotateAngle = -atand(p(1)) ;
            
        elseif strcmp(obj(oidx).Area,'AOTu')
            
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
    loadLayerMasks(obj(oidx).Layers)
    I = []; J = []; D = []; sel = []; mag=[];
    
    if ~isempty(obj(oidx).Layers(1).layerMask) && ~isempty(obj(oidx).Layers(end).layerMask) && UseMSPifLayerMasksEmpty~=2
        
        for Lidx = 1:length(obj(oidx).Layers)
            
            
            if isempty(obj(oidx).Layers(Lidx).layerMask)
                disp(['obj(' num2str(oidx) ').Layer ' num2str(Lidx) ' has no layerMask! Whole frame is included, which can distort results'])
            end
            
            
            % Find the (y,x) position of each polselective pixel in tuning map 'k'
            [i,j] = ind2sub(size(obj(oidx).Layers(Lidx).polPix),find(~isnan(obj(oidx).Layers(Lidx).polPix)));
            
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
                %         imrotateAngle = 0;
                Ii = i;
                Jj = j;
                %         ImicronPix = obj(oidx).Layers(Lidx).micronsPerPixel*(max(j) - min(j) + 1);
                %         JmicronPix = obj(oidx).Layers(Lidx).micronsPerPixel*(max(i) - min(i) + 1);
            end
            
            
            D = [D;obj(oidx).Layers(Lidx).polPix(find(~isnan(obj(oidx).Layers(Lidx).polPix)))];
            sel = [sel;obj(oidx).Layers(Lidx).polSelImg(find(~isnan(obj(oidx).Layers(Lidx).polPix)))];
            mag = [mag;obj(oidx).Layers(Lidx).fftMagImg(find(~isnan(obj(oidx).Layers(Lidx).polPix)))];
            
            I = [I;Ii];
            J = [J;Jj];
            
        end
        
    elseif ~isempty(obj(oidx).MIP) && UseMSPifLayerMasksEmpty>0
        loadLayerMasks(obj(oidx).MIP)
        
        if isempty(obj(oidx).MIP.layerMask)
            disp(['obj(' num2str(oidx) ').MIP has no layerMask! Whole frame is included, which can distort results'])
        elseif UseMSPifLayerMasksEmpty == 1
            disp(['obj(' num2str(oidx) ') Layers have no Masks - Using MSP instead'])
        end
        
        [layerI,layerJ,layerD,layersel,layermag,imrotateAngle]=getPixRotate(obj(oidx).MIP,rotateAngle,rotateOFF);
        
        D = [D;layerD];
        sel = [sel;layersel];
        mag = [mag;layermag];
        
        I = [I;layerI];
        J = [J;layerJ];
        
    else
        continue
    end
    
    
    % Scale to within range [0:1]
    X = (I - min(I))./(max(I) - min(I));
    Y = (J - min(J))./(max(J) - min(J));
    
    if strcmp(obj(oidx).Area,'AOTu') && leftFlag
        % Reverse X Y directions to pool with r.h.s. data
        X = 1-X;
        Y = 1-Y;
        leftFlag = 0;
        
    elseif strcmp(obj(oidx).Area,'Bu') && leftFlag
        X = 1-X;
        leftFlag = 0;
    end
    
    mag(isnan(mag)) = 0;
    sel(isnan(sel)) = 0;
    
    
    % Add to cross-object vectors:
    posX = [posX;X];
    posY = [posY;Y];
    angD = [angD;D];
    selW = [selW;sel];
    magW = [magW;mag];
    
    if clearMIPflag
        obj(oidx).MIP = [];
    end
    
    
    
end


posT = deg2rad(atan2d(posY-0.5,posX-0.5));


magW = magW./max(magW(:));

W = magW.*selW;

[W,sortIdx] = sort(W,'ascend');
posT = posT(sortIdx);
angD = angD(sortIdx);





% Display rotated image:
figure('units','normalized','position',[0.2 0.2 0.5 0.3])
addExportFigToolbar(gcf)
if darkMode
    colstr = 'w';
    acolstr = 'k';
    
else
    colstr = 'k';
    acolstr = 'w';
end
set(gcf,'color',acolstr)

if ~Yoff && ~Xoff
    subplot(1,3,1)
else
    subplot(1,2,1)
end

if isempty(obj(oidx).MIP) && length(length(obj(oidx).Layers)) == 1
    obj(oidx).MIP = obj(oidx).Layers(1);
    clearMIPflag = 1;
end
pst = obj(oidx).MIP.polSelThreshold;
obj(oidx).MIP.polSelThreshold = 0 ;

rotK = imrotate(obj(oidx).MIP.polPix,imrotateAngle);
obj(oidx).MIP.polSelThreshold = pst ;

imagesc(rotK);
if applyTuningCols
    colormap([[1 1 1];flipud(hsv)])
else
    colormap([[1 1 1];flipud(hsv)])
end
set(gca,'Color',colstr)

axis image
[Kj,Ki] = ind2sub(size(rotK),find(~isnan(rotK)&(rotK~=0)));

ylim([min(Kj),max(Kj)])
xlim([min(Ki),max(Ki)])

axis off
if clearMIPflag
    obj(oidx).MIP = [];
end



% Plot x- correlation
if ~Xoff
    if ~Yoff
        subplot(1,3,2)
    else
        subplot(1,2,2)
    end
end

%%% Create scatter plot:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if applyTuningCols % Colormap for each point based on tuning
    cm = flipud(hsv(18001));
    Yup = round(angD*100);
    C = cm(Yup+1,:);
elseif applyWeights % Colormap for each point based on selectivty
    if darkMode
        cm = [0 0 0;(magma(255))];
    else
        cm = [1 1 1;flipud(magma(255))];
    end
    Yup = round(W.*255);
    Yup(isnan(Yup)) = 0;
    C = cm(Yup+1,:);
else % No colormap
    C = [];
end
% Decimate arrays before plotting: limit to 1000 pts
if applyPtsLimit && length(posT) > ptsLimit
    f = floor(length(posT)/ptsLimit);
    posTdec = posT(1:f:end);
    angDdec = angD(1:f:end);
    if ~isempty(C)
        Cdec = C(1:f:end,:);
    end
    scatter(posTdec,angDdec,10,Cdec,'filled')
else
    
    scatter(posT,angD,10,C,'filled')
    
end
if darkMode
    set(gca,'color',[1 1 1])
end

%
% lincorr = corr(phi,x);
%
% if sign(lincorr)==-1 && sign(rho)==1
% rho = abs(rho)*sign(lincorr);
% amax = abs(amax)*sign(lincorr);
% end

%{
% % % Find correlation and plot
if strcmp(obj(oidx).Cell,'R7R8') || strcmp(obj(oidx).Cell,'MeTu')
    if leftFlag
        bounds = [0 2];
    else
        bounds = [-2 0];
    end
    [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(4*angD), posX,bounds, W);
    phi_fit =180*(amax/2)*posX + rad2deg(phi_0/4);
    
    
    
    hold on
    
    for phi_shift = [-360:90:360]
        
        this_phi_fit = phi_fit + phi_shift;
        if mod(phi_shift,180) == 90
            L = plot(posX(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),'--','color',colstr);
        else
            L = plot(posX(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
        end
    end
    %     plot(X,phi_fit+360,'k')
    %     plot(X,phi_fit+270,'k--')
    %
    %     plot(X,phi_fit+180,'k')
    %     plot(X,phi_fit+90,'k--')
    %
    %     plot(X,phi_fit,'k')
    %     plot(X,phi_fit-90,'k--')
    %
    %     plot(X,phi_fit-180,'k')
    %     plot(X,phi_fit-270,'k--')
    %
    %     plot(X,phi_fit-360,'k')
    
else
    if leftFlag
        bounds = [0 1];
    else
        bounds = [-1 0];
    end
    [rho, pval, phi_0, amax] = circ_corrcc(deg2rad(2*angD), posT); %bounds, W);
    phi_fit = 180*amax*posT + rad2deg(phi_0/2);
    
    
    
    hold on
    for phi_shift = [-360:180:360]
        
        this_phi_fit = phi_fit + phi_shift;
        L = plot(posT(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
        
    end
    %     plot(X,phi_fit+360,'k')
    %     plot(X,phi_fit+180,'k')
    %     plot(X,phi_fit,'k')
    %     plot(X,phi_fit-180,'k')
    %     plot(X,phi_fit-360,'k')
end
%}
%  set(gca,'xlim',[0 1])
%  set(gca,'xtick',[0 1])
set(gca,'ylim',[0 180])
set(gca,'ytick',[0:30:180])
set(gca,'yticklabel',{'0';'';'';'90';'';'';'180'})
set(gca,'xlim',[-pi pi])
set(gca,'xtick',[-pi,0,pi])
set(gca,'xticklabel',{'-pi';'0';'pi'});
xlabel(['circular position (rad)'])

set(gca,'Color',colstr)
% if strcmp(obj(oidx).Area,'Me') && leftFlag
%     xlabel('p  \rightarrow  a','Interpreter','tex')
% elseif strcmp(obj(oidx).Area,'Me')
%     xlabel('a  \rightarrow  p','Interpreter','tex')
% elseif leftFlag
%     xlabel('l  \rightarrow  m','Interpreter','tex')
% else
%     xlabel('m  \rightarrow  l','Interpreter','tex')
% end


ylabel('orientation preference (deg)')
offsetAxes
%title(['R = ' num2str(rho,2) ', p = ' num2str(pval)])
pbaspect([1,1,1])
ax = gca;
ax.XAxis.Color = colstr;
ax.YAxis.Color = colstr;
ax.Color = acolstr;
ax.FontWeight = 'bold';


% Plot y- correlation
if ~Yoff
    if Xoff
        subplot(1,2,2)
    else
        subplot(1,3,3)
    end
    
    % Scale to within range [0:1]
    posY = posY(sortIdx);
    
    %%% Create scatter plot:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if applyTuningCols % Colormap for each point based on tuning
        cm = flipud(hsv(18001));
        Yup = round(angD*100);
        C = cm(Yup+1,:);
    elseif applyWeights % Colormap for each point based on selectivty
        if darkMode
            cm = [0 0 0;(magma(255))];
        else
            cm = [1 1 1;flipud(magma(255))];
        end
        Yup = round(W.*255);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
    else % No colormap
        C = [];
    end
    % Decimate arrays before plotting: limit to 1000 pts
    if applyPtsLimit && length(posY) > ptsLimit
        f = floor(length(posY)/ptsLimit);
        posYdec = posY(1:f:end);
        angDdec = angD(1:f:end);
        if ~isempty(C)
            Cdec = C(1:f:end,:);
        end
        scatter(posYdec,angDdec,10,Cdec,'filled')
    else
        
        scatter(posY,angD,10,C,'filled')
        
    end
    
    if strcmp(obj(oidx).Cell,'R7R8') % || strcmp(obj(oidx).Cell,'MeTu')
        if leftFlag
            bounds = [0 2];
        else
            bounds = [-2 0];
        end
        [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(4*angD), posY,bounds, W);
        phi_fit =180*(amax/2)*posY + rad2deg(phi_0/4);
        
        
        
        hold on
        
        for phi_shift = [-360:90:360]
            
            this_phi_fit = phi_fit + phi_shift;
            if mod(phi_shift,180) == 90
                L = plot(posY(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),'--','color',colstr);
            else
                L = plot(posY(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
            end
        end
        %     plot(X,phi_fit+360,'k')
        %     plot(X,phi_fit+270,'k--')
        %
        %     plot(X,phi_fit+180,'k')
        %     plot(X,phi_fit+90,'k--')
        %
        %     plot(X,phi_fit,'k')
        %     plot(X,phi_fit-90,'k--')
        %
        %     plot(X,phi_fit-180,'k')
        %     plot(X,phi_fit-270,'k--')
        %
        %     plot(X,phi_fit-360,'k')
        
    else
        if leftFlag
            bounds = [0 1];
        else
            bounds = [-1 0];
        end
        [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*angD), posY,bounds, W);
        phi_fit = 180*amax*posY + rad2deg(phi_0/2);
        
        
        
        hold on
        for phi_shift = [-360:180:360]
            
            this_phi_fit = phi_fit + phi_shift;
            L = plot(posY(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
            
        end
    end
    
    set(gca,'ylim',[0 180])
    set(gca,'ytick',[0:30:180])
    set(gca,'yticklabel',{'0';'';'';'90';'';'';'180'})
    
    ax = gca;
    ax.XAxis.Color = colstr;
    ax.YAxis.Color = colstr;
    ax.Color = acolstr;
    ax.FontWeight = 'bold';
    xlabel('v  \rightarrow  d','Interpreter','tex')
    ylabel('orientation preference (deg)')
    offsetAxes
    title(['R = ' num2str(rho,2) ', p = ' num2str(pval)])
    pbaspect([1,1,1])
    
    
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
    %         imrotateAngle = 0;
    Ii = i;
    Jj = j;
    %         ImicronPix = obj(oidx).Layers(Lidx).micronsPerPixel*(max(j) - min(j) + 1);
    %         JmicronPix = obj(oidx).Layers(Lidx).micronsPerPixel*(max(i) - min(i) + 1);
end

D = obj.polPix(find(~isnan(obj.polPix)));
sel = obj.polSelImg(find(~isnan(obj.polPix)));
mag = obj.fftMagImg(find(~isnan(obj.polPix)));
end
