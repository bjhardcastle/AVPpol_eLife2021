function plotRThetacorr(obj,Yoff,Xoff,rotateOFF)
%
% For AVP MIP objects:
% Take pol selective pixels from the tuning map and plot x-position or
% y-position vs tuning angle. Find the circular-linear correlation and plot
% the line and R-value
%  Ypff    set to 1 to exclude the y-correlation plot

if nargin<2 || isempty(Yoff)
    Yoff = 1;
elseif ~any(Yoff==0 || Yoff==1)
    Yoff = 1;
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

if isempty(obj.ROI)
    getPolROIs(obj)
end

% Find the (y,x) position of each polselective pixel in tuning map 'k'
[j,i] = ind2sub(size(obj.k),find(~isnan(obj.k)));

% Find the centroid of each ROI mask 
for roiIdx = 1:length(obj.ROI)
    if ~isnan(obj.ROI(roiIdx).angle)
       rp = regionprops(obj.ROI(roiIdx).mask,'area','centroid');
       [~,sortIdx] = sort([rp.Area],'descend');
       rp = rp(sortIdx);
       
       angleList(roiIdx) = obj.ROI(roiIdx).angle;
       centroidList(roiIdx,:) = rp.Centroid;
       
       % Find mean response within k map.*ROI mask (in radians)
       respList(roiIdx) = wrapTo2Pi(circ_mean(deg2rad(2.*obj.k(obj.ROI(roiIdx).mask))))./2;

    end
end

% Find center of all centroids in (x,y)
ROIcenter(1) = mean(centroidList(:,1));
ROIcenter(2) = mean(centroidList(:,2));

% For each ROI centroid, convert to polar coordinates and keep its angle,
% theta, in radians
for angIdx = 1:length(angleList)
    thetaList(angIdx) = deg2rad(atan2d(centroidList(angIdx,2)-ROIcenter(2),centroidList(angIdx,1)-ROIcenter(1)));
end
[R,pval] = circ_corrcc(thetaList,respList);



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

    imrotateAngle = 90 + rotateAngle;
    
    ROIcenterRot(2) = ROIcenter(1)*cosd(rotateAngle) + ROIcenter(2)*sind(rotateAngle);
    ROIcenterRot(1) = -ROIcenter(1)*sind(rotateAngle) + ROIcenter(2)*cosd(rotateAngle);
    
else
    imrotateAngle = 0;
    J = j;
    I = i;
    ImicronPix = obj.micronsPerPixel*(max(i) - min(i) + 1);
    JmicronPix = obj.micronsPerPixel*(max(j) - min(j) + 1);
    ROIcenterRot(1) = ROIcenter(2);
    ROIcenterRot(2) = ROIcenter(1);
end


% Display rotated image:
figure,
if ~Yoff && ~Xoff
    subplot(1,3,1)
else
    subplot(1,2,1)
end
rotK = imrotate(obj.k,imrotateAngle);
imagesc(rotK);
colormap([[1 1 1];flipud(hsv)])
axis image
hold on
if rotateOFF
plot(ROIcenterRot(2),ROIcenterRot(1),'k+')
end
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
% X = (I - min(I))./(max(I) - min(I));
% Y = obj.k((~isnan(obj.k)));

% Get colors for each point based on tuning
cm = flipud(hsv(18001));
Yup = round(rad2deg(respList)*100);
C = cm(Yup+1,:);
scatter(thetaList,rad2deg(respList),100,C,'filled')

% Find correlation and plot
[rho, pval, phi_0, amax] = circlin_corr(respList',thetaList');
p = polyfit(thetaList,rad2deg(respList),1);
X = linspace(-pi,pi,100);
phi_fit = p(1)*X +p(2);
[rho,pval] = circ_corrcc(respList,thetaList);
% phi_fit = 180*amax*X/(pi) + rad2deg(phi_0);
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
set(gca,'xlim',[-pi pi])
set(gca,'xtick',[-pi,0,pi])
set(gca,'xticklabel',{'-\pi';'0';'\pi'});
xlabel(['circular position (rad)'])
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
    Y = obj.k((~isnan(obj.k)));
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
