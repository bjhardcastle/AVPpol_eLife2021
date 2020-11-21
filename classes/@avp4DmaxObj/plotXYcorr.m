function plotXYcorr(obj,Yoff,Xoff,rotateOFF)
%
% For AVP MIP objects:
% Take pol selective pixels from the tuning map and plot x-position or
% y-position vs tuning angle. Find the circular-linear correlation and plot
% the line and R-value
%  Ypff    set to 1 to exclude the y-correlation plot

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

if nargin<4 || isempty(rotateOFF)
    rotateOFF = 0;
elseif ~any(rotateOFF==0 || rotateOFF==1)
    rotateOFF = 0;
end

if isempty(obj.polPix)
    getPolROIs(obj)
end

% Find the (y,x) position of each polselective pixel in tuning map 'k'
[j,i] = ind2sub(size(obj.polPix),find(~isnan(obj.polPix)));

% Option to rotate all pixels first:
% Auto-rotation uses a linear fitted line through all points to find the
% orientation of the cells in the image, then rotates to align with x-axis.
% 90 degree rotation to correct for 2p orientation might be sufficient
% here.
if ~rotateOFF
    p = polyfit(j,i,1);
    
    rotateAngle = -atand(p(1));
    
    % We don't simply rotate the image becuase that method interpolates and
    % creates new data points
    J = i*cosd(rotateAngle) + j*sind(rotateAngle);
    I = -i*sind(rotateAngle) + j*cosd(rotateAngle);
    
    % Also modify the scale appropriately
    ImicronPix = obj.micronsPerPixel*(max(i) - min(i) + 1)*cosd(rotateAngle) + obj.micronsPerPixel*(max(j) - min(j) + 1)*sind(rotateAngle);
    JmicronPix = -obj.micronsPerPixel*(max(i) - min(i) + 1)*sind(rotateAngle) + obj.micronsPerPixel*(max(j) - min(j) + 1)*cosd(rotateAngle);

    imrotateAngle = 90 - rotateAngle;
else
    imrotateAngle = 0;
    J = j;
    I = i;
    ImicronPix = obj.micronsPerPixel*(max(i) - min(i) + 1);
    JmicronPix = obj.micronsPerPixel*(max(j) - min(j) + 1);
end


% Display rotated image:
figure,
if ~Yoff && ~Xoff
    subplot(1,3,1)
else
    subplot(1,2,1)
end
rotK = imrotate(obj.polPix,imrotateAngle);
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

% Scale to within range [0:1]
X = (I - min(I))./(max(I) - min(I));
Y = obj.polPix((~isnan(obj.polPix)));
% Get colors for each point based on tuning
cm = flipud(hsv(18001));
Yup = round(Y*100);
C = cm(Yup+1,:);

scatter(X,Y,10,C,'filled')

% Find correlation and plot
[rho, pval, phi_0, amax] = circlin_corr(deg2rad(2*Y), X);
phi_fit = 180*amax*X + rad2deg(0.5*phi_0);

hold on

plot(X,phi_fit+360,'k')
plot(X,phi_fit+180,'k')
plot(X,phi_fit,'k')
plot(X,phi_fit-180,'k')
plot(X,phi_fit-360,'k')

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
end


% Plot y- correlation
if ~Yoff
    if Xoff
        subplot(1,2,2)
    else
    subplot(1,3,3)
    end
    
    X = (J - min(J))./(max(J) -min(J));
    Y = obj.polPix((~isnan(obj.polPix)));
    cm = flipud(hsv(18001));
    Yup = round(Y*100);
    C = cm(Yup+1,:);
    
    
    scatter(X,Y,10,C,'filled')
    
    [rho, pval, phi_0, amax] = circlin_corr(deg2rad(2*Y), X);
    
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
end
