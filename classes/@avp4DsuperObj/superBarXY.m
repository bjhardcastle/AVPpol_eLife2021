function superBarXY(obj,Yoff,Xoff,rotateOFF, darkMode)
% For AVP super objects: (superclass version of plotXYcorr)
% Take pol selective pixels from the tuning map and plot x-position or
% y-position vs tuning angle. Find the circular-linear correlation and plot
% the line and R-value

applyTuningCols = 0; % Color scatter points according to pol tuning ( redundant info)
applyWeights = 1; % Weighted correlation using pol selectivity + normalized fft magnitude for each pixel
applyPtsLimit = 1;
ptsLimit = 20000;

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
    % First check that object included a bar experiment
    if ~any(obj(oidx).Layers(1).TrialSeqNum==5) 
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
    [i,j] = ind2sub(size(obj(oidx).MIP.barPix),find(~isnan(obj(oidx).MIP.barPix)));
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
    I = []; J = []; D = []; sel = []; 
    for Lidx = 1:length(obj(oidx).Layers)
        
        
        if isempty(obj(oidx).Layers(Lidx).layerMask)
            disp(['obj(' num2str(oidx) ').Layer ' num2str(Lidx) ' has no layerMask! Whole frame is included, which can distort results'])
        end
        
        
        % Find the (y,x) position of each polselective pixel in tuning map 'k'
        [i,j] = ind2sub(size(obj(oidx).Layers(Lidx).barPix),find(~isnan(obj(oidx).Layers(Lidx).barPix)));
        
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
        
        
        D = [D;obj(oidx).Layers(Lidx).barPix(find(~isnan(obj(oidx).Layers(Lidx).barPix)))];
        sel = [sel;obj(oidx).Layers(Lidx).barSelImg(find(~isnan(obj(oidx).Layers(Lidx).barPix)))];
        
        I = [I;Ii];
        J = [J;Jj];
        
    end
    % Scale to within range [0:1]
    X = (I - min(I))./(max(I) - min(I));
    Y = (J - min(J))./(max(J) - min(J));
    
    if strcmp(obj(oidx).Area,'AOTu') && leftFlag
       % Reverse X Y directions to pool with r.h.s. data
       X = 1-X;
       Y = 1-Y;
       leftFlag = 0;
    end
    
    sel(isnan(sel)) = 0;
    
    % if leftFlag
    %     X = 1-X; % reverse ant-post direction
    % end
    
    % Add to cross-object vectors:
    posX = [posX;X];
    posY = [posY;Y];
    angD = [angD;D];
    selW = [selW;sel];
    
    if clearMIPflag
        obj(oidx).MIP = [];
    end
    
    
    
end





W = selW;

[W,sortIdx] = sort(W,'ascend');
posX = posX(sortIdx);
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


rotK = imrotate(obj(oidx).MIP.barPix,imrotateAngle);
obj(oidx).MIP.polSelThreshold = pst ;

imagesc(rotK);
if applyTuningCols
    colormap([[1 1 1];flipud(aris)])
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
    cm = flipud(aris(18001));
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
if applyPtsLimit && length(posX) > ptsLimit
    f = floor(length(posX)/ptsLimit);
    posXdec = posX(1:f:end);
    angDdec = angD(1:f:end);
    if ~isempty(C)
        Cdec = C(1:f:end,:);
    end
    scatter(posXdec,angDdec,10,Cdec,'filled')
else
    
    scatter(posX,angD,10,C,'filled')
    
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
    [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*angD), posX,bounds, W);
    phi_fit = 180*amax*posX + rad2deg(phi_0/2);
    
    
    
%     hold on
%     for phi_shift = [-360:180:360]
%         
%         this_phi_fit = phi_fit + phi_shift;
%         L = plot(posX(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
%         
%     end
    %     plot(X,phi_fit+360,'k')
    %     plot(X,phi_fit+180,'k')
    %     plot(X,phi_fit,'k')
    %     plot(X,phi_fit-180,'k')
    %     plot(X,phi_fit-360,'k')
end

%  set(gca,'xlim',[0 1])
%  set(gca,'xtick',[0 1])
set(gca,'ylim',[-90 90])
set(gca,'ytick',[-90:45:90])
set(gca,'yticklabel',{'-90';'';'0';'';'90'})

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
if strcmp(obj(oidx).Area,'Me')
    xlabel('a  \rightarrow  p','Interpreter','tex')
else
    xlabel('m  \rightarrow  l','Interpreter','tex')
end

ylabel('orientation preference (\circ)')
offsetAxes
title(['R = ' num2str(rho,2) ', p = ' num2str(pval)])
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
        cm = flipud(aris(18001));
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
%     for phi_shift = [-360:180:360]
%         
%         this_phi_fit = phi_fit + phi_shift;
%         L = plot(posY(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
%         
%     end
    end
    
    set(gca,'ylim',[-90 90])
set(gca,'ytick',[-90:45:90])
set(gca,'yticklabel',{'-90';'';'0';'';'90'})

ax = gca;
ax.XAxis.Color = colstr;
ax.YAxis.Color = colstr;
ax.Color = acolstr;
ax.FontWeight = 'bold';
        xlabel('v  \rightarrow  d','Interpreter','tex')
    ylabel('orientation preference (\circ)')
    offsetAxes
    title(['R = ' num2str(rho,2) ', p = ' num2str(pval)])
    pbaspect([1,1,1])
    
    
end
