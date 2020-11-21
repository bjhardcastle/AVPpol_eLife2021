function plotStacks(objarray)

if isempty([objarray.TrialPatNum])
    runAcrossPlanes(objarray,'getParameters')
end

if isempty([objarray.ROI])
    assignStackROIs(objarray)
end

angles = [30 60 90];
angleStk = zeros(size(objarray(1).AverageFrame,1),size(objarray(1).AverageFrame,2),objarray(1).numZPlanes,length(angles));

for oidx = 1:length(objarray)
    objarray(oidx).Unattended = 1;
    getFrames(objarray(oidx))
    
    for aidx = 1:length(angles)
            % Get activity frame for each angle and put into 'angle stack'
            angleStk(:,:,objarray(oidx).ZPlane,aidx) = getPolAngleActivityFrame(objarray(oidx),angles(aidx));
            % SOMETIMES FAILS ON THIS LINE: JUST RUN AGAIN!
    end
    
    objarray(oidx).Frames = [];
    objarray(oidx).Unattended = 0;
end

for n = 1:2
    if n == 1
        nocrop = 0;
    elseif n ==2 
        nocrop = 1;
    end 
       
[boundingBox,planeLim] = getCropBox(objarray);

f = figure('Position',get(0, 'Screensize'), 'color',[1 1 1]);

for pidx = 1:(length(angles) + 1)
    
    if pidx <= length(angles)
                % draw line to indicate pol angles
               pax(pidx) = subplot(4,length(angles) + 1,pidx);
        lineRadius = 4;
        xyCenter = [0,0];       
        
        aop = angles(pidx);
        lineX = lineRadius*cosd(aop);
        lineY = lineRadius*sind(aop);
        xyEnd1 = xyCenter - [lineX, lineY];
        xyEnd2 = xyCenter + [lineX, lineY];
        L1 = line(pax(pidx) ,[xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);       
        L1.Color = 'r';
        L1.LineWidth = 3;
        
        aop = angles(pidx)+90;
        lineX = lineRadius*cosd(aop);
        lineY = lineRadius*sind(aop);
        xyEnd1 = xyCenter - [lineX, lineY];
        xyEnd2 = xyCenter + [lineX, lineY];
        L1 = line(pax(pidx) ,[xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
        L1.Color = 'b';
        L1.LineWidth = 3;
        
        hold(pax(pidx),'on')
        pax(pidx) .XLim = [-8 8];
                pax(pidx) .YLim = [-8 8];
                axis(pax(pidx) ,'equal')
                title(pax(pidx),['\fontsize{12}{\color{red}' num2str(angles(pidx)) '^{\circ}   \color{blue}' num2str(angles(pidx)+90) '^{\circ}}'],'interpreter','tex')

        axis(pax(pidx) ,'off')
        
        % Plot activity data
        ax(pidx) = subplot(4,length(angles)+1,[(length(angles)+1 +pidx),(2*length(angles)+2 +pidx)]);
        a = squeeze(angleStk(:,:,:,pidx));
        A = [];
        for sidx = 1:size(a,3)
             if ~nocrop
                 A(:,:,sidx) = imcrop(a(:,:,sidx),boundingBox);
             else
                 A(:,:,sidx) = a(:,:,sidx);
             end
        end
        
%             A(:,:,sidx) = a(:,:,sidx);            
            Anorm = flip(A./max(abs(A(:))),3); % Reverse z-direction
            Anorm = permute(Anorm,[2 1 3]);
%             Anorm = A;

        % Make (x,y,z) vectors for plotting data against
        Zpix = objarray(1).micronsPerZStep;
        XYpix = objarray(1).micronsPerPixel;
        xSlice = linspace(0,(size(Anorm,1)-1)*XYpix,size(Anorm,1));
        ySlice = linspace(0,(size(Anorm,2)-1)*XYpix,size(Anorm,2));
        zSlice = linspace(0,(size(Anorm,3)-1)*Zpix,size(Anorm,3));
        [xSliceVec,ySliceVec,zSliceVec] = meshgrid(ySlice,xSlice,zSlice);
%         
%         ax(pidx).DataAspectRatio = [XYpix2Z XYpix2Z 1];
%                 h = slice(Anorm,[],[],[1:objarray(1).numZPlanes]);

        h = slice(xSliceVec,ySliceVec,zSliceVec,Anorm,[],[],[zSlice]);
        
%         maskStack1  = double((a./max(a(:)))>0.1);
%         maskStack2  = double((a./min(a(:)))>0.1);
%         emptyMask = [];
%         for oidx = 1:length(objarray)
%             objarray(oidx).ROI4D(pidx).mask = maskStack1(:,:,oidx);
%             objarray(oidx).ROI4D(pidx + length(angles)).mask = maskStack2(:,:,oidx);
%             if ~sum(sum( maskStack1(:,:,oidx) ))  && ~sum(sum( maskStack2(:,:,oidx) ))
%                 emptyMask(oidx) = 0;
%             else
%                 emptyMask(oidx) = 1;
%             end
%         end
        
        %     figure
        %     h = slice( ,[],[],[1:10])
        %     alpha('color')
        %        set(h, 'EdgeColor','none')
        %     set(h, 'EdgeColor','none', 'FaceColor','interp','FaceAlpha','interp')
        %
        set(h, 'EdgeColor','none')
        set(h, 'EdgeColor','none', 'FaceColor','interp','FaceAlpha','interp')
        
        polarmap
        
        alpha('color')
%         alphamap('vdown')
        % Make sharper alphamap than 'vdown'
        rUp = exp([linspace(1,10,32)]);
        newAmap = round([ 1-(rUp./max(rUp)) fliplr(1-(rUp./max(rUp))) ],2);
        alphamap(ax(pidx),newAmap);
        
        % alphamap('decrease',.2)
        ax(pidx).ALim = [-1 1];
        % title(num2str(angles(pidx)))
%         title(['\fontsize{12}{\color{red}' num2str(angles(pidx)) '^{\circ}   \color{blue}' num2str(angles(pidx)+90) '^{\circ}}'],'interpreter','tex')
        XYpix2Z = objarray(1).micronsPerZStep/objarray(1).micronsPerPixel;
        ax(pidx).DataAspectRatio = [1 1 1];
        % ax(pidx).ZAxis.TickLabel = flipud(ax(pidx).ZAxis.TickLabel);
        ax(pidx).CLim = [-0.5 0.5];
        %     colormap([polarmap(64);[0.3 0.3 0.3]])
        
        ax(pidx).ZLim = [zSlice(planeLim(1)) zSlice(planeLim(end))];
        ax(pidx).YLim = [xSlice(1) xSlice(end)];
        ax(pidx).XLim = [ySlice(1) ySlice(end)];
        daspect(ax(pidx),[1 1 1])

% colorbar('horiz')
        % last and first are the 'wrong' way round because a was flipped.. 
    else
        
        ax(pidx) = subplot(4,length(angles)+1,[(length(angles)+1 +pidx),(2*length(angles)+2 +pidx)]);
        
        % Plot average image for anatomy of region, independent of activity
        % Plot average image for anatomy of region, independent of activity
        A2 = [];
        for sidx = 1:length(objarray)
            if ~nocrop
                A2(:,:,sidx) = imcrop(objarray(sidx).AverageFrame,boundingBox);
            else
                A2(:,:,sidx) = objarray(sidx).AverageFrame;
            end
        end
        A2norm = double(flip(( A2./max(abs(A2(:))) ) > 0.1,3));
        A2thresh = mean(A2(:)) + 0.5*std(A2(:));
        A2norm = double(flip( A2 > A2thresh,3));

            A2norm = permute(A2norm,[2 1 3]);
%             Anorm = A;

%        h2 = slice(xSliceVec,ySliceVec,zSliceVec,0.2.*A2norm,[],[],[1:objarray(1).numZPlanes]);
%        h2 = slice(0.2.*A2norm,[],[],[1:objarray(1).numZPlanes]);
 h2 = slice(xSliceVec,ySliceVec,zSliceVec,0.2.*A2norm,[],[],[zSlice]);
 
        set(h2, 'EdgeColor','none')
        set(h2, 'EdgeColor','none', 'FaceColor','interp','FaceAlpha','interp')
        %   alpha('color')
        pause(1)
        colormap(ax(pidx), flipud(gray))
        alpha(ax(pidx) ,'color')
        alphamap(ax(pidx), 'rampup')
        % alphamap('decrease',.2)
        ax(pidx).ALim = [0 0.3];
        ax(pidx).CLim = [0 1];
        % title(num2str(angles(pidx)))
%         XYpix2Z = objarray(1).micronsPerZStep/objarray(1).micronsPerPixel;
        ax(pidx).DataAspectRatio = [1 1 1];
        
        ax(pidx).ZLim = [zSlice(planeLim(1)) zSlice(planeLim(end))];
        ax(pidx).YLim = [xSlice(1) xSlice(end)];
        ax(pidx).XLim = [ySlice(1) ySlice(end)];
        
        % last and first are the 'wrong' way round because a was flipped.. 
    end
end
% % linkaxes(ax)
% hlink = linkprop([ax(1) ax(2) ax(3) ax(4)],{'XLim','YLim','ZLim','CameraPosition','CameraUpVector'});
% addprop(hlink,'PlotBoxAspectRatio')
% rotate3d on



% Add activity traces below 3d images
cols = [1 0 0;0 0 1];
for pidx = 1:length(angles)
    % Find the F0 region before Trial1 of the pol rotate experiment
    % Find all trials with SeqNum == 4
    firstTrial = find([objarray(1).TrialSeqNum == 4] ,1,'first');
    firstTrialFrame = objarray(1).TrialStartFrame(firstTrial);
    % Now subtract find range { firstTrial - 6s : firstTrial - 2s }
    FPS = objarray(1).AIrate/objarray(1).IFI;
    F0range = [floor(firstTrialFrame-(6*FPS)) floor(firstTrialFrame-(3*FPS)) ];
    
    fields = [];
    fields.TrialPatNum = 30;
    fields.TrialSeqNum = 4;
    cycStarts = objarray(1).TrialStartFrame(findTrials(objarray(1),fields));
    fields.TrialPatNum = 360;
    cycEnds = objarray(1).TrialEndFrame(findTrials(objarray(1),fields));
    
    ACax(pidx) = subplot(4,length(angles)+1,3*length(angles)+3+pidx);
    
    % method 1: find dF/F0, then plot both traces, or their mean
    % method 2: find F, then find the mean trace for each roi, and
    % normalize to the max of both means
    
    plotMethod = 2;
    switch plotMethod
        case 1
            for n = 1:2
                roiIdx = pidx + (n-1)*length(angles);
                F = scanROI4D(objarray,roiIdx);
                F(F0range(1):F0range(2))
                F0 = mean(F(F0range(1):F0range(2)));
                Ftrace = (F./F0)-1;
                cycTrace = [];
                cycTrace(:,1) = Ftrace(cycStarts(1):cycEnds(1));
                cycTrace(:,2) = Ftrace(cycStarts(2):cycStarts(2)+length(cycTrace(:,1))-1);
                
                
                plot(cycTrace,'color',cols(n,:),'LineWidth',1)
                hold on
            end
            
        case 2
            
            for n = 1:2
                meanTrace = [];
                
                roiIdx = pidx;
                F = scanROI4D(objarray,roiIdx);
                F(F0range(1):F0range(2))
                F0 = mean(F(F0range(1):F0range(2)));
                Ftrace = F-F0;
                cycTrace = [];
                cycTrace(:,1) = Ftrace(cycStarts(1):cycEnds(1));
                cycTrace(:,2) = Ftrace(cycStarts(2):cycStarts(2)+length(cycTrace(:,1))-1);
                meanTrace(:,1) = mean(cycTrace,2);
                
                plot(cycTrace(:,1)./(max(meanTrace(:,1)) - min(meanTrace(:,1))),'color',cols(1,:),'LineWidth',1)
                plot(cycTrace(:,2)./(max(meanTrace(:,1)) - min(meanTrace(:,1))),'color',cols(1,:),'LineWidth',1)

                hold on
                roiIdx = pidx + length(angles);
                F = scanROI4D(objarray,roiIdx);
                F(F0range(1):F0range(2))
                F0 = mean(F(F0range(1):F0range(2)));
                Ftrace = F-F0;
                cycTrace = [];
                cycTrace(:,1) = Ftrace(cycStarts(1):cycEnds(1));
                cycTrace(:,2) = Ftrace(cycStarts(2):cycStarts(2)+length(cycTrace(:,1))-1);
                meanTrace(:,2) = median(cycTrace,2);

                
%                 plot(cycTrace(:,1)./(max(meanTrace(:,1)) - min(meanTrace(:,1))),'color',cols(1,:),'LineWidth',1)
%                 plot(cycTrace(:,2)./(max(meanTrace(:,2)) - min(meanTrace(:,2))),'color',cols(2,:),'LineWidth',1)



                plot(cycTrace(:,1)./(max(meanTrace(:,2)) - min(meanTrace(:,2))),'color',cols(2,:),'LineWidth',1)
                plot(cycTrace(:,2)./(max(meanTrace(:,2)) - min(meanTrace(:,2))),'color',cols(2,:),'LineWidth',1)

                hold on
            end
            
    end

axis square
axis off
ACax(pidx).YAxis.Visible = 'off';
ACax(pidx).XAxis.Visible = 'on';
offsetAxesXTickONLY

angVector = nan(1,size(cycTrace,1));
testAngles = [30:30:360];
angleSpacing = linspace(0.5,length(testAngles)+0.5, length(angVector));
for an = 1:length(testAngles)
    [~,mIdx] = min(abs(angleSpacing - an));
    angVector(mIdx) = testAngles(an);
end
prefAngle1 = angles(pidx);
prefAngle2 = mod(angles(pidx)+180,360);
if prefAngle2 == 0
    prefAngle2 = 360;
end
aprefAngle1 = mod(angles(pidx)+90,360);
aprefAngle2 = mod(angles(pidx)+270,360);
if aprefAngle2 == 0
    aprefAngle2 = 360;
end



currXMarkers = ACax(pidx).XAxis.TickValues;

m1 = find(angVector == prefAngle1);
m2 = find(angVector == aprefAngle1);
m3 = find(angVector == prefAngle2);
m4 = find(angVector == aprefAngle2);

ACax(pidx).XAxis.TickValues = [m1 m2 m3 m4];
ACax(pidx).XAxis.TickLabels = {num2str(prefAngle1);num2str(aprefAngle1);num2str(prefAngle2);num2str(aprefAngle2)};
ACax(pidx).XAxis.FontWeight ='bold';
ACax(pidx).XAxis.FontSize = 12;

line(ACax(pidx),[1 length(angVector)], [0 0],'color',[0.5 0.5 0.5],'LineStyle',':')
ACax(pidx).XAxis.TickLength(1) = 2*ACax(pidx).XAxis.TickLength(1);

end
linkaxes([ACax])


for xidx = 1:4    
%     ax(xidx).FontSize = 8;
% 
%     if contains(objarray(1).File, ' L ')
%         ax(xidx).XAxis.Label.String = 'l \leftrightarrow m';        
%     else        
%         ax(xidx).XAxis.Label.String = 'm \leftrightarrow l';
%     end
%     ax(xidx).ZAxis.Label.String = 'v \leftrightarrow d';
%     ax(xidx).YAxis.Label.String = 'a \leftrightarrow p';
%         ax(xidx).FontWeight = 'bold';
%     ax(xidx).YAxis.Visible = 'on';
%         ax(xidx).XAxis.Visible = 'on';
% h3(xidx) = rotate3d(ax(xidx)); 
% set(h3(xidx), 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)') 
% % set(h3(xidx), 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')') 
%    set(h3(xidx),'ActionPostCallback',@axislabel_rotation);
% 
% set(f, 'ResizeFcn', @align_axislabel) 
% align_axislabel([], ax(xidx)) 
% % axislabel_rotation([],ax(xidx))
% % axislabel_translation_slider;
% view(ax(xidx),2)
% pause(2)
% view(ax(xidx),3)

        if contains(objarray(1).File, ' L ')
            axesAnatomyLabelsL(ax(xidx))
        else 
            axesAnatomyLabelsR(ax(xidx))
        end
% 
% Add scalebar3
% xScale = objarray(1).micronsPerPixel;
% yScale = objarray(1).micronsPerPixel;
% zScale = objarray(1).micronsPerZStep;
% scalebar3(ax(xidx),5)

end


% linkaxes(ax)
hlink = linkprop([ax(1) ax(2) ax(3) ax(4)],{'XLim','YLim','ZLim','CameraPosition','CameraUpVector'});
addprop(hlink,'PlotBoxAspectRatio')
% rotate3d on

%% PRINT

% filename = 'selectivityXYZ';
%     if n == 1
%         filename = [filename '_crop'];
%     end 
% set(gcf, 'InvertHardCopy', 'off');
% print(gcf, fullfile(objarray(1).Folder,filename) , '-dpng','-r450');
% 
% view(ax(1),2)
% view(ax(2),2)
% view(ax(3),2)
% view(ax(4),2)
% 
% filename = 'selectivityXY';
%     if n == 1
%         filename = [filename '_crop'];
%     end 
% print(gcf, fullfile(objarray(1).Folder,filename) , '-dpng','-r450');
%  close(gcf)



end
end

function scalebar3(ax,micronLength)
% Work out appropriate number of microns by taking a proportion of the
% plot limits
xrange = ax.XLim(2) - ax.XLim(1) + 1;
yrange = ax.YLim(2) - ax.YLim(1) + 1;
zrange = ax.ZLim(2) - ax.ZLim(1) + 1;

[maxVal,maxIdx] = max([xrange,yrange]);
pixLength = 0.2*maxVal;
% Convert pixels into microns, then round to closest 5 
% if maxIdx == 1
%     micronLength = round(pixLength*xs /5) *5;
% else
%     micronLength = round(pixLength*ys /5) *5;
% end
% mux = micronLength/xs;
% muy = micronLength/ys;
% muz = micronLength/zs;

mux = micronLength;
muy = micronLength;
muz = micronLength;

 x1 = ax.XLim(1);
 y1 = ax.YLim(2);
 z1 = ax.ZLim(1);
 x2 = x1 + mux;
 y2 = y1 - muy;
 z2 = z1 + muz;
  x = [x1,x2, x1,x1 , x1,x1  ];
  y = [y1,y1  , y1,y2,  y1,y1 ];
  z = [z1,z1  , z1,z1 , z1,z2];    
  
 
  yl = ax.YLim;
  
  
  hl = line(ax,x,y,z);
  hl.LineWidth = 3;
  hl.Color = [0.4 0.4 0.4];
  hl.MarkerFaceColor = [0.3 0.3 0.3];
  ht = text(x1-0.5*mux,y1-0.5*muy,z1,[num2str(micronLength), ' \mum']);
  ht.FontSize = ax.FontSize;
  ht.FontWeight = 'bold';
  ht.Color = [0 0 0];
  ht.VerticalAlignment = 'bottom';
  ht.HorizontalAlignment = 'left';
  ax.YLim = yl;
  ax.Clipping = 'off';
  

% if ~event.hasListener(ax, 'MarkedClean')
    
     addlistener (ax, 'MarkedClean', @(obj,event)deleteVertex(ax));
% end
end


function deleteVertex ( ax )


% extract the x axis vertext data
% X, Y and Z row of the start and end of the individual axle.
ax.XRuler.Axle.VertexData = single(zeros(3,2));
ax.YRuler.Axle.VertexData= single(zeros(3,2));
ax.ZRuler.Axle.VertexData= single(zeros(3,2));

% ax.XRuler.Axle.VertexData(1,2) = 0;
% ax.XRuler.Axle.VertexData(1,1) = 0 ;
% ax.XRuler.Axle.VertexData(3,1) = 0;
% ax.XRuler.Axle.VertexData(3,2) = 0;
% 
% 
% if strcmp(xoff, 'xoff')
%     ax.XRuler.Axle.VertexData(1,1) = 0;
%     ax.XRuler.Axle.VertexData(1,2) = 0;
% end

% % repeat for Y (set 2nd row)
% ax.YRuler.Axle.VertexData(2,1) = 0;
% ax.YRuler.Axle.VertexData(2,2) = 0;
% ax.ZRuler.Axle.VertexData(2,1) = 0;
% ax.ZRuler.Axle.VertexData(2,2) = 0;

ax.XTickLabel = [];
ax.YTickLabel = [];
ax.ZTickLabel = [];

end

%%
%{

f = figure('Position',get(0, 'Screensize'), 'color',[1 1 1]);

for pidx = 1:length(angles)
    clear h h2
    % Plot activity data
    ax(pidx) = subplot(1,length(angles),pidx);
    a = squeeze(angleStk(:,:,:,pidx));
    A = [];
    for sidx = 1:size(a,3)
        A(:,:,sidx) = imcrop(a(:,:,sidx),boundingBox);
        Anorm = flip(A./max(abs(A(:))),3);
    end
    h = slice(Anorm,[],[],[1:objarray(1).numZPlanes]);
    set(h, 'EdgeColor','none')
    set(h, 'EdgeColor','none', 'FaceColor','interp','FaceAlpha','interp')
        
    alpha('color')
    alphamap('vdown')
    % alphamap('decrease',.2)
    ax(pidx).ALim = [-1 1];
    % title(num2str(angles(pidx)))
    title(['\fontsize{12}{\color{blue}' num2str(angles(pidx)) '^{\circ}   \color{red}' num2str(angles(pidx)+90) '^{\circ}}'],'interpreter','tex')
    XYpix2Z = objarray(1).micronsPerZStep/objarray(1).micronsPerPixel;
    ax(pidx).DataAspectRatio = [XYpix2Z XYpix2Z 1];
    ax(pidx).YAxis.Visible = 'on';
    ax(pidx).XAxis.Visible = 'on';
    % ax(pidx).ZAxis.TickLabel = flipud(ax(pidx).ZAxis.TickLabel);

    hold on

    % Plot average image for anatomy of region, independent of activity
    A2 = [];
    for sidx = 1:length(objarray)
        A2(:,:,sidx) = imcrop(objarray(sidx).AverageFrame,boundingBox);
        A2norm = double(flip(( A2./max(abs(A2(:))) ) > 0.1,3));
    end
    h2 = slice(A2norm,[],[],[1:objarray(1).numZPlanes]);
    set(h2, 'EdgeColor','none')
    set(h2, 'EdgeColor','none', 'FaceColor','interp','FaceAlpha','interp')
    %   alpha('color')
    alphamap('vdown')
    
    
    colormap([polarmap(64); flipud(gray(64))])
   
    for hidx = 1:length(h2)
        
        
        % Below from:
    % https://www.mathworks.com/matlabcentral/answers/101346-how-do-i-use-multiple-colormaps-in-a-single-figure
    % Move the pcolor to Z = -10.
    % The 0*Z is in the statement below to insure that the size
    % of the ZData does not change.
%     set(h2(hidx),'CData',1 )
  
 
           % Scale the CData (Color Data) of each plot so that the
    % plots have contiguous, nonoverlapping values.  The range
    % of each CData should be equal. Here the CDatas are mapped
    % to integer values so that they are easier to manage;
    % however, this is not necessary.
    % Initially, both CDatas are equal to Z.
    m = 64;  % 64-elements is each colormap
    Avec = Anorm(:,:,hidx);
    A2vec = A2norm(:,:,hidx);
    cmin = min([ Avec(:);A2vec(:) ]);
    cmax = max([ Avec(:);A2vec(:) ]);
    % CData for surface
    C1 = min(m,round((m-1)*(Anorm(:,:,hidx)-cmin)/(cmax-cmin))+1);
    % CData for pcolor
    C2 = 64+C1;
    % Update the CDatas for each object.
    set(h(hidx),'CData',C1);
    set(h2(hidx),'CData',C2);
    
    % Change the CLim property of axes so that it spans the
    % CDatas of both objects.
    caxis([min(C1(:)) max(C2(:))])
    
    end
    



end
linkaxes(ax)
hlink = linkprop([ax(1) ax(2) ax(3)],{'CameraPosition','CameraUpVector'});
addprop(hlink,'PlotBoxAspectRatio')
rotate3d on


%}

