function superMIPXY(obj,Yoff,Xoff,rotateOFF, darkMode)
% For AVP super objects: (superclass version of plotXYcorr)
% Take pol selective pixels from the tuning map and plot x-position or
% y-position vs tuning angle. Find the circular-linear correlation and plot
% the line and R-value

applyTuningCols = 0; % Color scatter points according to pol tuning ( redundant info)
applyWeights = 1; % Weighted correlation using pol selectivity + normalized fft magnitude for each pixel
applyPtsLimit = 1;
ptsLimit = 1000;
runPermutationTest = 0;
convertAngles4Pub = 1;
ptSize = 1;
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
    % First check that object included a pol tuning experiment
    if isempty(obj(oidx).Layers(1).TrialSeqNum)
        getParameters(obj(oidx).Layers(1))
    end
%     if ~any(obj(oidx).Layers(1).TrialSeqNum==4) && ~any(obj(oidx).Layers(1).TrialSeqNum==2) && ~any(obj(oidx).Layers(1).TrialSeqNum==8)
    if ~any(obj(oidx).Layers(1).TrialSeqNum==4) && ~any(obj(oidx).Layers(1).TrialSeqNum==8)

        disp('No pol exp (4,8) available')
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
%     for Lidx = 1:length(obj(oidx).Layers)
        
        
        if isempty(obj(oidx).MIP.layerMask)
            disp(['obj(' num2str(oidx) ').MIP has no layerMask! Whole frame is included, which can distort results'])
        end
        
        
        % Find the (y,x) position of each polselective pixel in tuning map 'k'
        [i,j] = ind2sub(size(obj(oidx).MIP.polPix),find(~isnan(obj(oidx).MIP.polPix)));
        
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
        
        
        D = [D;obj(oidx).MIP.polPix(find(~isnan(obj(oidx).MIP.polPix)))];
        sel = [sel;obj(oidx).MIP.polSelImg(find(~isnan(obj(oidx).MIP.polPix)))];
        mag = [mag;obj(oidx).MIP.fftMagImg(find(~isnan(obj(oidx).MIP.polPix)))];
        
        I = [I;Ii];
        J = [J;Jj];
        
%     end
    % Scale to within range [0:1]
    X = (I - min(I))./(max(I) - min(I));
    Y = (J - min(J))./(max(J) - min(J));
    
    if (strcmp(obj(oidx).Area,'AOTu') || strcmp(obj(oidx).Area,'AOTu') )&& leftFlag
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
    
    % if leftFlag
    %     X = 1-X; % reverse ant-post direction
    % end
    
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




magW = magW./max(magW(:));

W = magW.*selW;

[W,sortIdx] = sort(W,'ascend');
posX = posX(sortIdx);
angD = angD(sortIdx);

origAngD = angD;
if convertAngles4Pub
    angD = rec2weir(angD);
end




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
    Yup = round(origAngD*100);
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
    cm = ones(18001,3).*ROIcolor(1);
    Yup = round(origAngD*100);
    C = cm(Yup+1,:);
end


% Decimate arrays before plotting: limit to 1000 pts
if applyPtsLimit && length(posX) > ptsLimit
    f = floor(length(posX)/ptsLimit);
    posXdec = posX(1:f:end);
    angDdec = angD(1:f:end);
    if ~isempty(C)
        Cdec = C(1:f:end,:);
    end
    scatter(posXdec,angDdec,ptSize,Cdec,'filled')
else
    scatter(posX,angD,ptSize,C,'filled')
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
if strcmp(obj(oidx).Cell,'R7R8') %|| strcmp(obj(oidx).Cell,'MeTu')
    if leftFlag
        bounds = [0 2];
    else
        bounds = [-2 0];
    end
    [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(4*angD), posX,bounds, W);
    phi_fit =180*(amax/2)*posX + rad2deg(phi_0/4);
     if runPermutationTest
    [rho, pval, phi_0, amax,perm_pval] = circlin_wcorr_permute(deg2rad(4*angD), posX,bounds, W,[]);
     else
         perm_pval = 999;
     end
    
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
elseif strcmp(obj(oidx).Area,'PB')

        bounds = [-2 2];

    posX2 = 2*posX;
    [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*angD), mod(posX2,1),bounds, W);
    phi_fit = 180*(1*amax)*posX*2 + rad2deg(phi_0/2);
     if runPermutationTest
    [rho, pval, phi_0, amax,perm_pval] = circlin_wcorr_permute(deg2rad(8*angD), posX,bounds, W,[]);
     else
         perm_pval = 999;
     end
    
    hold on
    
    for phi_shift = [-1080:180:1080]
        
        this_phi_fit = phi_fit + phi_shift;
            L = plot(posX(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
        
    end
else
    if leftFlag
        bounds = [0 1];
    else
        bounds = [-1 0];
    end
    [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*angD), posX,bounds, W);
    phi_fit = 180*amax*posX + rad2deg(phi_0/2);
    
    if runPermutationTest
        [rho, pval, phi_0, amax,perm_pval] = circlin_wcorr_permute(deg2rad(2*angD), posX,bounds, W,[]);
        
    else
        perm_pval = 999;
    end
    
    
    hold on
    for phi_shift = [-360:180:360]
        
        this_phi_fit = phi_fit + phi_shift;
        L = plot(posX(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
        
    end
    %     plot(X,phi_fit+360,'k')
    %     plot(X,phi_fit+180,'k')
    %     plot(X,phi_fit,'k')
    %     plot(X,phi_fit-180,'k')
    %     plot(X,phi_fit-360,'k')
end

%  set(gca,'xlim',[0 1])
%  set(gca,'xtick',[0 1])
if convertAngles4Pub
    set(gca,'ylim',[-90 90])
    set(gca,'ytick',[-90:30:90])
    set(gca,'yticklabel',{'-90';'';'';'0';'';'';'90'})
    
else
    set(gca,'ylim',[0 180])
    set(gca,'ytick',[0:30:180])
    set(gca,'yticklabel',{'0';'';'';'90';'';'';'180'})

end

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
    xlabel('a  ->  p','Interpreter','tex')
else
    xlabel('m  ->  l','Interpreter','tex')
end

ylabel('orientation preference (deg)')
offsetAxes
    title(['R = ' num2str(rho,2) ', p = ' num2str(pval) ', permute p = ' num2str(perm_pval)])
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
        Yup = round(origAngD*100);
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
        cm = ones(18001,3).*ROIcolor(1);
    Yup = round(origAngD*100);
    C = cm(Yup+1,:);
    end
    
    % Decimate arrays before plotting: limit to 1000 pts
    if applyPtsLimit && length(posY) > ptsLimit
        f = floor(length(posY)/ptsLimit);
        posYdec = posY(1:f:end);
        angDdec = angD(1:f:end);
        if ~isempty(C)
            Cdec = C(1:f:end,:);
        end
        scatter(posYdec,angDdec,ptSize,Cdec,'filled')
    else
        
        scatter(posY,angD,ptSize,C,'filled')
        
    end
    
    if strcmp(obj(oidx).Cell,'R7R8') % || strcmp(obj(oidx).Cell,'MeTu')
    if leftFlag
        bounds = [0 2];
    else
        bounds = [-2 0];
    end
     if runPermutationTest
    [rho, pval, phi_0, amax,perm_pval] = circlin_wcorr_permute(deg2rad(4*angD), posY,bounds, W);
     else 
         perm_pval = 999;
     end
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
        bounds = [-1 1];
end
     if runPermutationTest
    [rho, pval, phi_0, amax,perm_pval] = circlin_wcorr_permute(deg2rad(2*angD), posY,bounds, W);
     else
         perm_pval = 999;
     end
    phi_fit = 180*amax*posY + rad2deg(phi_0/2);
    
    
    
    hold on
    for phi_shift = [-360:180:360]
        
        this_phi_fit = phi_fit + phi_shift;
        L = plot(posY(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),colstr);
        
    end
    end
    
if convertAngles4Pub
    set(gca,'ylim',[-90 90])
    set(gca,'ytick',[-90:30:90])
    set(gca,'yticklabel',{'-90';'';'';'0';'';'';'90'})
    
else
    set(gca,'ylim',[0 180])
    set(gca,'ytick',[0:30:180])
    set(gca,'yticklabel',{'0';'';'';'90';'';'';'180'})

end

ax = gca;
ax.XAxis.Color = colstr;
ax.YAxis.Color = colstr;
ax.Color = acolstr;
ax.FontWeight = 'bold';
        xlabel('v  ->  d','Interpreter','tex')
    ylabel('orientation preference (deg)')
    offsetAxes
    title(['R = ' num2str(rho,2) ', p = ' num2str(pval) ', permute p = ' num2str(perm_pval)])
    pbaspect([1,1,1])
    
    
end
end
function [rho, pval, phi_0, amax,perm_pval] = circlin_wcorr_permute(phi, x, bounds, weights,nreps,plotflag)
% with help from https://courses.washington.edu/matlab1/Bootstrap_examples.html#1
if nargin < 5 || isempty(nreps)
    nreps = 10000;
end
if nargin <6 || isempty(plotflag)
    plotflag =0;
end
perm = zeros(nreps,1);
for i=1:nreps
    if mod(i,10000)==0
        disp(['permutation count:' num2str(i)])
    end
    %shuffle the lsat scores and recalculate the correlation
    perm(i) = circlin_wcorr(phi,x(randperm(length(x))),bounds,weights);
end
[rho, pval, phi_0, amax] = circlin_wcorr(phi,x,bounds,weights);
perm_pval = sum(abs(perm)>abs(rho))/nreps;
if plotflag
figure(99)
hist(perm,50);
ylim = get(gca,'YLim');
hold on
plot(rho*[1,1],ylim,'r-','LineWidth',2);
xlabel('Correlation');
title(sprintf('Proportion of permutations exceeding observed correlation: %5.5f',perm_pval));
end
end
function [anglesPublished] = rec2weir(anglesRecorded)
    anglesPublished = wrapTo180(-anglesRecorded-270);
end