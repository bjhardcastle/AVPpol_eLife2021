function superPCA2(obj)
if isempty(obj.Layers(2).k)
    Layersummary(obj)
end
s = floor(obj.MIP.micronsPerZStep./obj.MIP.micronsPerPixel);
% s=1;
clear R4
for oidx = 1:length(obj.Layers)
    R4(:,:,oidx) = obj.Layers(oidx).polPix;%(1:s:end,1:s:end);
end

[I,J,K] = ind2sub(size(R4),find(~isnan(R4)));



%%
cols = obj.Layers(1).polColMap;
   
P = nan(length(vertcat(I(:))) , 4);
angs = obj.Layers(1).expPolAngs(4);

    % load parameters
    % X position:
    P(:,1) = J;
    % Y position:
    P(:,2) = I;
    % Z position:
    P(:,3) = K;
    % angle idx:
    P(:,4) = deg2rad(R4(~isnan(R4(:))));
    
%
[rX,pX,phi0X,aX] = circlin_corr(P(:,4),P(:,1));
[rY,pY,phi0Y,aY] = circlin_corr(P(:,4),P(:,2));
[rZ,pZ,phi0Z,aZ] = circlin_corr(P(:,4),P(:,3));

% shuffle same data and find correlation coeficient
for n = 1:1000
    sIdx = randperm(size(P,1));
Ps = P;
Ps(:,1:3) = (P(sIdx,1:3));
[rXt,pX] = circ_corrcl(Ps(:,4),Ps(:,1));
[rYt,pY] = circ_corrcl(Ps(:,4),Ps(:,2));
[rZt,pZ] = circ_corrcl(Ps(:,4),Ps(:,3));
rXs(n) = rXt;
rYs(n) = rYt;
rZs(n) = rZt;

end

figure,
ax1 = subplot(1,3,1);
histogram(rXs,20)
    pbaspect([1 1 1])
    xlim([0 1])
line([rX rX],[ylim],'color','r')

ax2 =subplot(1,3,2);
histogram(rYs,20)
    pbaspect([1 1 1])
        xlim([0 1])

line([rY rY],[ylim],'color','r')

ax3 = subplot(1,3,3);
histogram(rZs,20)
    pbaspect([1 1 1])
        xlim([0 1])

line([rZ rZ],[ylim],'color','r')

% if any([pX pY]<0.1) && any(abs([rX rY])>0)
    if isempty(obj.MIP.fullAngMaskImg)
        getPolROIs(obj.MIP)
    end
    figure,image(obj.MIP.fullAngMaskImg)
    axis image
    XYratio = abs(rX/rY);
    midpoints = floor(size(obj.MIP.fullAngMaskImg)./2);
    %     if abs(XYratio) > 1 % X is larger
    Xscale = rX;
    %         Yscale = sign(rY)/XYratio;
    %     elseif abs(XYratio) < 1 % Y is larger
    %         Xscale = sign(rX)*XYratio;
    Yscale = rY;
    %     end
    halflengths = floor(midpoints(1:2)./4);
    
    line([ midpoints(2) midpoints(2)+Xscale*halflengths(2)], ...
        [ midpoints(1) midpoints(1)+Yscale*halflengths(1)], ...
        'LineWidth',2,'color','k')
    line([ midpoints(2) midpoints(2)-Xscale*halflengths(2)], ...
        [ midpoints(1) midpoints(1)-Yscale*halflengths(1)], ...
        'LineWidth',2,'color','k')
    
    if sign(rX) > 0
        markerX = max([midpoints(2)-Xscale*halflengths(2) midpoints(2)+Xscale*halflengths(2)]);
    else
        markerX = min([midpoints(2)-Xscale*halflengths(2) midpoints(2)+Xscale*halflengths(2)]);
    end
    if sign(rY) > 0
        markerY = max([midpoints(1)-Yscale*halflengths(1) midpoints(1)+Yscale*halflengths(1)]);
    else
        markerY = min([midpoints(1)-Yscale*halflengths(1) midpoints(1)+Yscale*halflengths(1)]);
    end
    hold on,
    plot(markerX,markerY,'kp','MarkerSize',12,'MarkerFaceColor','k');
% end


varWeights = [1,1,1,1];
[coeff,score,~,~,explained] =  pca(P,'VariableWeights',varWeights);

% Pc = score;
% figure
% hold on
% for n = 1:6
%     if n ==1
%         idx1 = 1;
%     else
%         idx1 = length(vertcat(I{1:n-1})) + 1;
%     end
%     idx2 = idx1 + length(I{n}) -1;
%
%     plot(Pc(idx1:idx2,1),Pc(idx1:idx2,2),'o','color',cols(n,:),'markerfacecolor',cols(n,:))
% end
%

figure,hbi = biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'X','Y','Z','ang'},'ObsLabels',num2str((1:size(score,1))') );
 
cm = flipud(hsv(18001));
    Yup = round(P(:,4)*100);
    C = cm(Yup+1,:);
    
for ii = 1:length(hbi)
    userdata = hbi(ii).UserData;
    if ~isempty(userdata) && strcmp( hbi(ii).Tag, 'obsmarker')
        angIdx = int32( rad2deg(P(userdata,4)));
        hbi(ii).Color=cm(round(angIdx*100)+1,:);
        hbi(ii).MarkerSize = 12;
    else
        set(hbi(ii), 'Color', 'k');
        
    end
end
title(['weights: [x y z ang] [' num2str(varWeights) ']'])
set(gcf,'color','w');

%
% dirCoeff = coeff(1:3,1);
% dirCoeff = dirCoeff./max(dirCoeff(:));
% figure,title('xyz coeffs')
% stem([1 2 3],dirCoeff)
% % hold on
% % stem([1 2 3],coeff(1:3,2))
% set(gca,'xlim',[0.5 3.5],'ylim',[-1 1],'XTick',[1 2 3],'XTickLabel',{'x';'y';'z'})
% pbaspect([1 4 1])

