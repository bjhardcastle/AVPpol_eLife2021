function makeStackVid(objarray)

if isempty([objarray.TrialPatNum])
    runAcrossPlanes(objarray,'getParameters')
end

if isempty([objarray.ROI])
    assignStackROIs(objarray)
end

loops = 4;


% Open figure for displaying/capturing frames
% VidFig = figure('Color',[0 0 0]);
% set(VidFig,'units','normalized' ...
%     , 'outerposition',[0 0 1 1]  );

angles = [30 60 90 120 150 180 ];
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

for cidx = 1:2  % cropped then not cropped
    for midx = 1:2  % both +ve + -ve, then positive enchanced data only
        for nidx = 1:2 % 2 different views
            
            % Get path for saving video
            saveFile = '4stack';
        
            
            
            if midx == 1
                saveFile = [saveFile '_sep'];
            elseif midx == 2
                saveFile = [saveFile '_comb'];
            end
            if nidx == 1
                saveFile = [saveFile '_xy'];
            elseif nidx == 2
                saveFile = [saveFile '_xyz'];
            end
            if cidx == 1
                nocrop = 0;
                [boundingBox,planeLim] = getCropBox(objarray);
                 saveFile = [saveFile '_crop'];
            elseif cidx == 2
                nocrop = 1;
                boundingBox = [];
                planeLim= [];
            end
            
            savePath = objarray(1).Folder;
            % Make video file and open
            v = VideoWriter(fullfile(savePath,saveFile),'Motion JPEG AVI');
            % v.FrameRate = round(obj.AIrate/obj.IFI);
            v.FrameRate = 4;
            v.Quality = 100;
            open(v)
            
            VidFig = figure('Position',get(0, 'Screensize'), 'color',[1 1 1]);
            
            aidx = 1;
            ax(2) = subplot(1,3,1);
            ax(aidx) = subplot(1,3,2:3);

            for tidx = 1:loops
                for pidx = 1:(length(angles))
                    
                    if pidx <= length(angles)
                        
                        % Plot activity data
                        a = squeeze(angleStk(:,:,:,pidx));
                        A = [];
                        for sidx = 1:size(a,3)
                            if ~nocrop
                                A(:,:,sidx) = imcrop(a(:,:,sidx),boundingBox);
                            else
                                A(:,:,sidx) = a(:,:,sidx);
                            end
                        end
                        
                        Anorm = flip(A./max(abs(A(:))),3); % Reverse z-direction
                        Anorm = permute(Anorm,[2 1 3]);
                        
                        % Only show positive, enhanced response
                        if midx == 1
                            Anorm = Anorm.*(Anorm>0);
                        end
                        %{
        % Original
        h = slice(Anorm,[],[],[1:objarray(1).numZPlanes]);

                        %}
                        
                        % Make (x,y,z) vectors for plotting data against
                        Zpix = objarray(1).micronsPerZStep;
                        XYpix = objarray(1).micronsPerPixel;
                        xSlice = linspace(0,(size(Anorm,1)-1)*XYpix,size(Anorm,1));
                        ySlice = linspace(0,(size(Anorm,2)-1)*XYpix,size(Anorm,2));
                        zSlice = linspace(0,(size(Anorm,3)-1)*Zpix,size(Anorm,3));
                        [xSliceVec,ySliceVec,zSliceVec] = meshgrid(ySlice,xSlice,zSlice);
                        
                        h = slice(ax(1),xSliceVec,ySliceVec,zSliceVec,Anorm,[],[],[zSlice]);
                        
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
                        % colormap('parula')
                        
                        alpha('color')
                        %         alphamap('vdown')
                        % Make sharper alphamap than 'vdown'
                        rUp = exp([linspace(1,10,32)]);
                        newAmap = round([ 1-(rUp./max(rUp)) fliplr(1-(rUp./max(rUp))) ],2);
                        alphamap(ax(aidx),newAmap);
                        
                        % alphamap('decrease',.2)
                        ax(aidx).ALim = [-0.5 0.5];
                        % title(num2str(angles(pidx)))
                        %         title(['\fontsize{12}{\color{red}' num2str(angles(pidx)) '^{\circ}   \color{blue}' num2str(angles(pidx)+90) '^{\circ}}'],'interpreter','tex')
                        XYpix2Z = objarray(1).micronsPerZStep/objarray(1).micronsPerPixel;
                        ax(aidx).YAxis.Visible = 'on';
                        ax(aidx).XAxis.Visible = 'on';
                        ax(aidx).ZAxis.Visible = 'on';
                        % ax(pidx).ZAxis.TickLabel = flipud(ax(pidx).ZAxis.TickLabel);
                        ax(aidx).CLim = [-0.5 0.5];
                        %     colormap([polarmap(64);[0.3 0.3 0.3]])
                        
                        %{
        %Original
        ax(aidx).ZLim = planeLim;
        ax(aidx).XLim = [0 size(Anorm,1)];
        ax(aidx).YLim = [0 size(Anorm,2)];
                ax(aidx).DataAspectRatio = [XYpix2Z XYpix2Z 1];

                        %}
                        % New
                        if ~nocrop
                        ax(aidx).ZLim = [zSlice(planeLim(1)) zSlice(planeLim(end))];
                        end
                        ax(aidx).YLim = [xSlice(1) xSlice(end)];
                        ax(aidx).XLim = [ySlice(1) ySlice(end)];
                        ax(aidx).DataAspectRatio = [1 1 1];
                        
                        if contains(objarray(1).File, ' L ')
                            axesAnatomyLabelsL(ax(aidx))
                        else
                            axesAnatomyLabelsR(ax(aidx))
                        end
                        
                        if nidx == 1
                            view(2)
                        end
                        % colorbar('horiz')
                        % last and first are the 'wrong' way round because a was flipped..
                        %{
    else
        
        ax(pidx) = subplot(2,length(angles)+1,pidx);
        
        % Plot average image for anatomy of region, independent of activity
        % Plot average image for anatomy of region, independent of activity
        A2 = [];
        for sidx = 1:length(objarray)
            if ~nocrop
                A2(:,:,sidx) = imcrop(objarray(sidx).AverageFrame,boundingBox);
            else
                A2(:,:,sidx) =objarray(sidx).AverageFrame;
            end
        end
        A2norm = double(flip(( A2./max(abs(A2(:))) ) > 0.1,3));
        A2thresh = mean(A2(:)) + 0.5*std(A2(:));
        A2norm = double(flip( A2 > A2thresh,3));

            A2norm = permute(A2norm,[2 1 3]);
%             Anorm = A;

        h2 = slice(0.2.*A2norm,[],[],[1:objarray(1).numZPlanes]);
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
        XYpix2Z = objarray(1).micronsPerZStep/objarray(1).micronsPerPixel;
        ax(pidx).DataAspectRatio = [XYpix2Z XYpix2Z 1];
        ax(pidx).YAxis.Visible = 'off';
        ax(pidx).XAxis.Visible = 'off';
        
        ax(pidx).ZLim = planeLim;
        ax(pidx).XLim = [0 size(A2norm,1)];
        ax(pidx).YLim = [0 size(A2norm,2)];
        
        view(2)
        % last and first are the 'wrong' way round because a was flipped..
    end
                        %}
                        
                        % % linkaxes(ax)
                        % hlink = linkprop([ax(1) ax(2) ax(3) ax(4)],{'XLim','YLim','ZLim','CameraPosition','CameraUpVector'});
                        % addprop(hlink,'PlotBoxAspectRatio')
                        % rotate3d on
                        
                        
                        %{
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
    
    ACax(pidx) = subplot(2,length(angles)+1,length(angles)+1+pidx);
    
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
                
                roiIdx = pidx + length(angles);
                F = scanROI4D(objarray,roiIdx);
                F(F0range(1):F0range(2))
                F0 = mean(F(F0range(1):F0range(2)));
                Ftrace = F-F0;
                cycTrace = [];
                cycTrace(:,1) = Ftrace(cycStarts(1):cycEnds(1));
                cycTrace(:,2) = Ftrace(cycStarts(2):cycStarts(2)+length(cycTrace(:,1))-1);
                meanTrace(:,2) = mean(cycTrace,2);

                plot(meanTrace(:,1)./max(meanTrace(:,1)),'color',cols(1,:),'LineWidth',1)
                plot(meanTrace(:,2)./max(meanTrace(:,2)),'color',cols(2,:),'LineWidth',1)

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



% linkaxes(ax)
hlink = linkprop([ax(1) ax(2) ax(3) ax(4)],{'XLim','YLim','ZLim','CameraPosition','CameraUpVector'});
addprop(hlink,'PlotBoxAspectRatio')
rotate3d on
                        %}
                        
                        % Also draw line to indicate pol angle
                        cla(ax(2));
                        hold(ax(2),'off')
                        
                        lineRadius = 4;
                        xyCenter = [5,0];
                        
                        aop = angles(pidx);

                        lineX = lineRadius*cosd(aop);
                        lineY = lineRadius*sind(aop);
                        xyEnd1 = xyCenter - [lineX, lineY];
                        xyEnd2 = xyCenter + [lineX, lineY];
                        
                        %             L = line([xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                        L = line(ax(2),[xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                        
                        L.Color = 'r';
                        L.LineWidth = 4;
                        
                        if midx == 2                           
                            
                            aop = angles(pidx)+90;                          
                          
                            lineX = lineRadius*cosd(aop);
                            lineY = lineRadius*sind(aop);
                            xyEnd1 = xyCenter - [lineX, lineY];
                            xyEnd2 = xyCenter + [lineX, lineY];
                            
                            %             L = line([xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                            L = line(ax(2),[xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                            
                            L.Color = 'b';
                            L.LineWidth = 4;
                        end
                        
                        ax(2).XLim = [-10 10];
                        ax(2).YLim = [-10 10];
                        axis(ax(2),'equal')
                        axis(ax(2),'off')
                        
                    end
                    
                    % Grab plot frame and write to video object
                    frame = getframe(VidFig);
                    writeVideo(v,frame);
                    
                end
            end
            
            close(VidFig)
            close(v)
        end
                end
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
    
