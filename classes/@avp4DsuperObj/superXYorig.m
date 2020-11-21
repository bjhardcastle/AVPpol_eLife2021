function superXYorig(obj,Yoff,Xoff,rotateOFF)
%
% For AVP super objects: (superclass version of plotXYcorr)
% Take pol selective pixels from the tuning map and plot x-position or
% y-position vs tuning angle. Find the circular-linear correlation and plot
% the line and R-value
%  Ypff    set to 1 to exclude the y-correlation plot
if nargin<4 || isempty(rotateOFF)
    rotateOFF = 0;
elseif ~any(rotateOFF==0 || rotateOFF==1)
    rotateOFF = 0;
end
clearMIPflag = 0;
if isempty(obj.MIP) && length(length(obj.Layers)) == 1
   obj.MIP = obj.Layers(1);
   clearMIPflag = 1;
end

pst = obj.MIP.polSelThreshold;
obj.MIP.polSelThreshold = 0 ; 
[j,i] = ind2sub(size(obj.MIP.polPix),find(~isnan(obj.MIP.polPix)));
obj.MIP.polSelThreshold = pst ; 
if contains(obj.Name,' L ')
    leftFlag = 1;
else
    leftFlag = 0;
end
if ~rotateOFF
   
    if strcmp(obj.Area,'Me')
        p = polyfit(j,i,1);
        rotateAngle = -atand(p(1)) ;
        
    elseif strcmp(obj.Area,'AOTu')
        
        rotateAngle = 20;
        
        if leftFlag
            rotateAngle = rotateAngle*-1;
        end
        
    end
else
    rotateAngle = 0;
end




% Also modify the scale appropriately
ImicronPix = obj.MIP.micronsPerPixel*(max(i) - min(i) + 1)*cosd(rotateAngle) + obj.MIP.micronsPerPixel*(max(j) - min(j) + 1)*sind(rotateAngle);
JmicronPix = -obj.MIP.micronsPerPixel*(max(i) - min(i) + 1)*sind(rotateAngle) + obj.MIP.micronsPerPixel*(max(j) - min(j) + 1)*cosd(rotateAngle);

Ii = [];Y = []; Ysel = []; Ymag=[];  Jj = [];
for Lidx = 1:length(obj.Layers)
    if nargin<2 || isempty(Yoff)
        Yoff = 0;
    elseif ~any(Yoff==0 || Yoff==1)
        Yoff = 0;
    end
    
    if nargin<3 || isempty(Xoff)
        Xoff = 0;
    elseif ~any(Xoff==0 || Xoff==1)
        Xoff = 0;
    end
    
    if isempty(obj.Layers(Lidx).layerMask)
       disp(['Layer ' num2str(Lidx) ' has no layerMask! Whole frame is included, which can distort results']) 
    end
    
    
    % Find the (y,x) position of each polselective pixel in tuning map 'k'
    [j,i] = ind2sub(size(obj.Layers(Lidx).polPix),find(~isnan(obj.Layers(Lidx).polPix)));
    
    % Option to rotate all pixels first:
    % Auto-rotation uses a linear fitted line through all points to find the
    % orientation of the cells in the image, then rotates to align with x-axis.
    % 90 degree rotation to correct for 2p orientation might be sufficient
    % here.
    if ~rotateOFF
        
        
        % We don't simply rotate the image becuase that method interpolates and
        % creates new data points
        J = [i*cosd(rotateAngle) + j*sind(rotateAngle)];
        I = [-i*sind(rotateAngle) + j*cosd(rotateAngle)];
        
        %
        %         % Also modify the scale appropriately
        %         ImicronPix = obj.Layers(Lidx).micronsPerPixel*(max(i) - min(i) + 1)*cosd(rotateAngle) + obj.Layers(Lidx).micronsPerPixel*(max(j) - min(j) + 1)*sind(rotateAngle);
        %         JmicronPix = -obj.Layers(Lidx).micronsPerPixel*(max(i) - min(i) + 1)*sind(rotateAngle) + obj.Layers(Lidx).micronsPerPixel*(max(j) - min(j) + 1)*cosd(rotateAngle);
        
        imrotateAngle = rotateAngle + 90;
    else
%         imrotateAngle = 0;
        J = j;
        I = i;
        ImicronPix = obj.Layers(Lidx).micronsPerPixel*(max(i) - min(i) + 1);
        JmicronPix = obj.Layers(Lidx).micronsPerPixel*(max(j) - min(j) + 1);
    end
    
    
    Y = [Y;obj.Layers(Lidx).polPix(find(~isnan(obj.Layers(Lidx).polPix)))];
    Ysel = [Ysel;obj.Layers(Lidx).polSelImg(find(~isnan(obj.Layers(Lidx).polPix)))];
    Ymag = [Ymag;obj.Layers(Lidx).fftMagImg(find(~isnan(obj.Layers(Lidx).polPix)))];

    Ii = [Ii;I];
    Jj = [Jj;J];
    
end


% Scale to within range [0:1]
X = (Ii - min(Ii))./(max(Ii) - min(Ii));
% if leftFlag
%     X = 1-X; % reverse ant-post direction
% end
Ymag(isnan(Ymag)) = 0;
Ysel(isnan(Ysel)) = 0;

Ymag = Ymag./max(Ymag(:));

W = Ymag.*Ysel;

[W,sortIdx] = sort(W,'ascend');
X = X(sortIdx);
Y = Y(sortIdx);

% Display rotated image:
figure,
if ~Yoff && ~Xoff
    subplot(1,3,1)
else
    subplot(1,2,1)
end

pst = obj.MIP.polSelThreshold;
obj.MIP.polSelThreshold = 0 ; 

rotK = imrotate(obj.MIP.polPix,imrotateAngle);
obj.MIP.polSelThreshold = pst ; 

imagesc(rotK);
colormap([[1 1 1];flipud(hsv)])
axis image
[Kj,Ki] = ind2sub(size(rotK),find(~isnan(rotK)&(rotK~=0)));

ylim([min(Kj),max(Kj)])
xlim([min(Ki),max(Ki)])

axis off

% Plot x- correlation
if ~Xoff
    if ~Yoff
        subplot(1,3,2)
    else
        subplot(1,2,2)
    end
end


% Get colors for each point based on tuning
cm = flipud(hsv(18001));
Yup = round(Y*100);
C = cm(Yup+1,:);
nocol=1;
if nocol
%  cm = [1 1 1;flipud(magma(255))];
% Yup = round(W.*255);
% Yup(isnan(Yup)) = 0;
% C = cm(Yup+1,:);      
% scatter(X,Y,10,'filled')
   scatter(X,Y,10,C,'filled')

else
    
   scatter(X,Y,10,C,'filled')

end

% 
% lincorr = corr(phi,x);
% 
% if sign(lincorr)==-1 && sign(rho)==1
% rho = abs(rho)*sign(lincorr);
% amax = abs(amax)*sign(lincorr);
% end


% % % Find correlation and plot
if strcmp(obj.Cell,'R7R8')
    if leftFlag
        bounds = [0 2];
    else
        bounds = [-2 0];
    end
    [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(4*Y), X,bounds, W);
    phi_fit =180*(amax/2)*X + rad2deg(phi_0/4);
    
    
    
    hold on
    
    for phi_shift = [-360:90:360]
        
        this_phi_fit = phi_fit + phi_shift;
        if mod(phi_shift,180) == 90
            L = plot(X(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),'k--');
        else
            L = plot(X(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),'k');
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
    [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*Y), X,bounds, W);
    phi_fit = 180*amax*X + rad2deg(phi_0/2);
    
    
    
    hold on
      for phi_shift = [-360:180:360]
        
        this_phi_fit = phi_fit + phi_shift;
        L = plot(X(this_phi_fit>=0 & this_phi_fit<=180),this_phi_fit(this_phi_fit>=0 & this_phi_fit<=180),'k');

    end
%     plot(X,phi_fit+360,'k')
%     plot(X,phi_fit+180,'k')
%     plot(X,phi_fit,'k')
%     plot(X,phi_fit-180,'k')
%     plot(X,phi_fit-360,'k')
end

%  set(gca,'xlim',[0 1])
%  set(gca,'xtick',[0 1])
set(gca,'ylim',[0 180])
set(gca,'ytick',[0:60:180])
xlabel(['X position (\times' num2str(ImicronPix,4) '\mum)'])
ylabel('orientation preference (\circ)')
set(gcf,'color','w');
offsetAxes
title(['R = ' num2str(rho,2) ', p = ' num2str(pval)])
pbaspect([1,1,1])



% Plot y- correlation
if ~Yoff
    if Xoff
        subplot(1,2,2)
    else
        subplot(1,3,3)
    end
    
    % Scale to within range [0:1]
    X = (Jj - min(Jj))./(max(Jj) - min(Jj));
    X = X(sortIdx);

    cm = flipud(hsv(18001));
    Yup = round(Y*100);
    C = cm(Yup+1,:);
    
    nocol=1;
    if nocol
%         cm = [1 1 1;flipud(magma(255))];
%         Yup = round(W.*255);
%         C = cm(Yup+1,:);
        % scatter(X,Y,10,'filled')
        scatter(X,Y,10,C,'filled')
        
    else
        
        scatter(X,Y,10,C,'filled')
        
    end

        
        bounds = [-1 1];
    [rho, pval, phi_0, amax] = circlin_wcorr(deg2rad(2*Y), X,bounds, W);
    
    hold on
    
    phi_fit = 180*amax*X + rad2deg(0.5*phi_0);  % linear model of circular variable
    
    plot(X,phi_fit+360,'k')
    plot(X,phi_fit+180,'k')
    plot(X,phi_fit,'k')
    plot(X,phi_fit-180,'k')
    plot(X,phi_fit-360,'k')
    
    set(gca,'ylim',[0 180])
    set(gca,'ytick',[0:60:180])
    xlabel(['Y position (\times' num2str(JmicronPix,4) '\mum)'])
    ylabel('orientation preference (\circ)')
    set(gcf,'color','w');
    offsetAxes
    title(['R = ' num2str(rho,2) ', p = ' num2str(pval)])
    pbaspect([1,1,1])

    
    if clearMIPflag
        obj.MIP = [];
    end
end
