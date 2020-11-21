function polSel = findPolSelImg(obj)

if isempty(obj.Frames)
    getFrames(obj)
end
if isempty(obj.TrialSeqNum)
    obj.getDaqData;
    getParameters(obj)
    obj.Daq = [];
end
  f = median(obj.Frames,3);
  
%Get array of responses for every pixel:

indLength = size(obj.Frames,1)*size(obj.Frames,2);
rows = reshape(repmat(1:size(obj.Frames,1),1,size(obj.Frames,2)),indLength,1);
cols = reshape(repmat(1:size(obj.Frames,2),size(obj.Frames,1),1),indLength,1);
inds = sub2ind(size(obj.Frames(:,:,1)), rows, cols);

resp = reshape(obj.Frames, indLength, size(obj.Frames,3));
resp(find(f<median(f(:))),:) = 0;

% Construct array [numPix x numPolAngles/2], into which we'll enter the
% average response for each pol angle, for every pixel:
if any(obj.TrialSeqNum==4)
    expNum = 4;
elseif any(obj.TrialSeqNum==2)
    expNum = 2;
else
    return
    error('No pol tuning exp (2 or 4) available')
end

polAngs = unique(obj.TrialPatNum(obj.TrialSeqNum == expNum));
pixPolResp = nan(indLength,0.5*length(polAngs));

for aIdx = 1:0.5*length(polAngs)
    angleFrameArray = getAngleLogicArray(obj,polAngs(aIdx));
    pixPolResp(:,aIdx) = mean(resp.*angleFrameArray,2);
end

% Find polarization selectivity value for each pixel:
[respMax,maxIdx] = max(pixPolResp,[],2);
maxIdx = maxIdx - 1;
respMin = pixPolResp( sub2ind(size(pixPolResp), [1:indLength]', 1+mod(maxIdx + 3,6) ) );

%%
% Threshold pixels by max intensity ***(improve this: use deltaF or std)***
% sortMax = sort(respMax);
% intThresh = sortMax(floor(length(sortMax)/10));
% respMin(respMax<intThresh) = 0;
% respMax(respMax<intThresh) = 0;

%%

polSelVec = ( respMax - respMin )./respMax;

polSelImg = reshape(polSelVec,size(obj.Frames(:,:,1)));

polColMap = hsv(0.5*length(polAngs));
polColMap = [0.*ones(1,3); polColMap];

polAngImg = reshape(maxIdx+1,size(obj.Frames(:,:,1)));
polAngImg = polAngImg.*(polSelImg>0.5);
% polAngImg(find(median(resp,2)<intThresh)) = 0;
% figure,image((polAngImg+1).*(polSelImg>0.1))
% axis('image')
% colormap(polColMap)

% figure,imshow(polSelImg)



    
  

polStart = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first'));
polStop = obj.TrialEndFrame(find(obj.TrialSeqNum == expNum,1,'last'));

autoROI = struct('mask',[],'color',[],'position',[],'response',[]);



saveplot =0;
switch saveplot
    case 0
        
        mFig = figure('Color',[0.5 0.5 0.5],'SizeChangedFcn',@resetAxes);
        subplot(2,2,3)
        avgAx = gca;
        imagesc(f)
        axis off
        
        
        subplot(2,2,1)
        maskAx = gca;
        axis off
        maskAx.Color = [0 0 0];
        
       
        

         dotAx = subplot(2,2,4);
        hold(dotAx,'on');
        dotAx.Color = [0.9 0.9 0.9];

         subplot(2,2,[2])
        mAx = gca;
        axis(mAx,'square');
        hold on
        
        % Add listener to keep images at max size
        % addlistener(mAx, 'Resize', @(obj,event)resetAxes(maskAx,avgAx));
        
        ridx = 0;
        u = unique(polAngImg);
        u=u(u>0);
        r = zeros( size(f,1),size(f,2),length(u) );
        masks = zeros( size(f,1),size(f,2));
        masksbw = zeros( size(f,1),size(f,2));
        maskscale = masksbw;
        
%         polAngImg = smoothPolAngImg(polAngImg);

        
        f=double(f);
        ctr = 0;
        for l = 1:length(u)
            rk_scale =  mean(mean(f(polAngImg==u(l))))./max(f(:));
            if rk_scale*max(f(:)) > 0%median(f(:))
                ridx = ridx+1;
                rk = polAngImg==u(l);
                col = ROIcolor(ridx);
                maskcol = cat(3,col(1)*rk,col(2)*rk,col(3)*rk);
                masks = masks + maskcol;
                masksbw = u(l).*rk + masksbw;
                maskscale = rk_scale.*rk + maskscale;
                r(:,:,ridx) = rk;
                numHalfCycles = 2*length(find(obj.TrialPatNum(obj.TrialSeqNum == expNum) == 360));
                maskResponse = scanROI(obj,rk,[polStart, polStop]);
                CycleLims = floor(linspace(1,length(maskResponse),numHalfCycles+1));
              maskTrialResp = nan(numHalfCycles,max(diff(CycleLims)));
                for cycIdx = 1:numHalfCycles
                    thisResp = maskResponse(CycleLims(cycIdx):CycleLims(cycIdx+1)-1);
                   % Subtract off mean:
                   maskTrialResp(cycIdx,1:CycleLims(cycIdx+1)-CycleLims(cycIdx)) = thisResp - nanmean(thisResp); 
             

                end
               angleVector = linspace(30,180,size(maskTrialResp,2));
%                 plot(mAx,angleVector,nanmean(maskTrialResp),'Color',col,'LineWidth',3)
%                 hold on
%                 plot(mAx,angleVector,(maskTrialResp)','Color',col,'LineWidth',1)
                lineProps.col = {col};
                mAx;
                mseb(angleVector,nanmean(maskTrialResp),nanstd(maskTrialResp),lineProps,1);            
                
                ctr = ctr + 1;
                polSel(ctr) = mean(polSelImg(find(rk)));

            end
        end
        ylim = mAx.YLim(2);
        mAx.XLim = [30 180];
%         if ~isempty(obj.TrialStartFrame)
%             for t = 1:length(obj.TrialStartFrame)
% %                 plotStimPatch(obj,mAx,[obj.TrialStartFrame(t) obj.TrialEndFrame(t)])
%                 xpos = obj.TrialStartFrame(t) + 0.5*(obj.TrialEndFrame(t) - obj.TrialStartFrame(t));
%                 if ~isempty(obj.TrialPatNum)
%                     text(mAx,xpos,ylim,[num2str(obj.TrialPatNum(t))],'HorizontalAlignment','center');
%                 end
%             end
%         end
        
        for kS = 1:ctr
            plot(dotAx,kS,polSel(kS),'o','color',ROIcolor(kS),'MarkerFaceColor',ROIcolor(kS))
            hold on
        end
        dotAx.YLim = [0 1];


        roiIdx = 0;
        addflag = 1;
        addedmasks = zeros( size(f,1),size(f,2)); % store ROIs as they're added
        
        
 
        
        try
            while addflag == 1
                
                roiIdx = roiIdx + 1;
                
                continueflag = 1;
                while continueflag == 1
                    
                    %%% pick contiguous mask area as an ROI
                    
                    imshow(masks + cat(3,addedmasks*0.3,addedmasks*0.3,addedmasks*0.3) ,'Parent',maskAx)
                    
                    title(maskAx,'Select neighboring areas to join to make a single ROI, then press return');
                    
                    % Add any existing ROIs
                    
%  maskbwAng = masksbw.*30;  
%  maskbwAngVector = reshape(maskbwAng, indLength, 1);
% maskbwAngVectorData = maskbwAngVector(maskbwAngVector>0);
%  clear maskParams
%  maskParams(:,1) = sind(2*maskbwAngVectorData);
%  maskParams(:,2) = cosd(2*maskbwAngVectorData);
% rowData = rows( maskbwAngVector>0);
% rowScaled = 2.*(rowData - min(rowData(:)))./ (max(rowData(:)) - min(rowData(:)) );
%  maskParams(:,3) = 1.*rowScaled;
% colData = cols( maskbwAngVector>0);
% colScaled = 2.*(colData - min(colData(:)))./ (max(colData(:)) - min(colData(:)) );
%  maskParams(:,4) = 1.*colScaled;
%                   
%  
% eva = evalclusters(maskParams,'kmeans','CalinskiHarabasz','KList',[1:10]);
% 
% % Run kmeans 
%  kblank = zeros(size(maskbwAngVector,1),1);
% kidx = kmeans(maskParams,eva.OptimalK);
% kblank(maskbwAngVector>0) = kidx;
% 
%  k = reshape(kblank,size(f));
%         figure,image(k)
% colormap(polColMap)   
% axis('image')


                    W = [];
                    
                    getptsflag = 1;
                    while getptsflag == 1
                        [x,y] = getpts(maskAx);
                        
                        
                        
                        % only accept if within bounds of image
                        if any(any( [round(y) round(x)] > [size(f)] ))
                            h = warndlg('Select pixels within top-left image only.','Out of bounds') ;
                            uiwait(h)
                            continue
                        end
                        
                        
                        % only accept if within a non-zero area of masks ( above threshold, not background )
                        % and only accept if not-previously selected for an ROI
                        clear inds
                        inds(:,1) = sub2ind( size(masks) , round(y), round(x) , 1*ones(length(x),1) );
                        inds(:,2) = sub2ind( size(masks) , round(y), round(x) , 2*ones(length(x),1) );
                        inds(:,3) = sub2ind( size(masks) , round(y), round(x) , 3*ones(length(x),1) );
                        masksinds = masks(inds);
                        maskInt = sum(masksinds,2);
                        addInt = addedmasks(inds(:,1));
                        if any( ~maskInt  ) || any( addInt )
                            
                            h = warndlg('Select colored areas only','Out of bounds') ;
                            uiwait(h)
                            continue
                        end
                        
                        
                        
                        if ~isempty(x) && ~isempty(y)
                            getptsflag = 0;
                        end
                    end
                    
                    
                    
                    for pidx = 1:length(x)
                        W(:,:,pidx) = grayconnected(masksbw,round(y(pidx)),round(x(pidx)));
                    end
                    
                    mask = logical(sum(W,3));
                    
                    BW2 = imfill(mask,'holes');
                    %%% Make ROI a bit smaller
                    %     BW = bwmorph(bwconvhull(BW2), 'erode', 2);
                    B = bwboundaries(BW2);
                    
                    cla(maskAx)
                    
                    currentMask = addedmasks + 2*BW2;
                    imshowpair(currentMask,obj.AverageFrame,'Parent',maskAx)
                    title(maskAx,'New ROI')
                    
                    if ~isempty(B)
                        roix = [];
                        roiy = [];
                        for b = 1:length(B)
                            roix = [roix B{b}(:,2)'];
                            roiy = [roiy B{b}(:,1)'];
                        end
                        skipPts = 5;
                        roiPos = [ roix([1:skipPts:end])' roiy([1:skipPts:end])'];
                    end
                    
                    % Confirm ROI
                    qstring = 'Happy with this ROI?';
                    wintitle = 'Accept new ROI?';
                    str2 = 'Choose again';
                    str3 = 'Accept';
                    str1 = 'Reject and finish';
                    default = str3;
                    button = questdlg(qstring,wintitle,str1,str2,str3,default);
                    if strcmp(button,str3)
                        continueflag = 0;
                    elseif strcmp(button,str1)
                        continueflag = 0;
                        addflag = 0;
                    end
                    
                    
                end
                
                if addflag
                    % remove just-selected ROI from masks
                    masks = masks.*~mask;
                    masksbw = masksbw.*~mask;
                    
                    % and add it to the 'addedmasks'
                    addedmasks = addedmasks + mask;
                    
                    [bw] = impoly(maskAx,roiPos);
                    
                    autoROI(roiIdx).mask = bw.createMask;
                    autoROI(roiIdx).color = [1 1 1];
                    autoROI(roiIdx).position = roiPos;
                    autoROI(roiIdx).response = [];
                    
                    % Confirm ROI
                    qstring = 'Continue adding?';
                    
                    wintitle = 'Add more ROIs?';
                    str1 = 'Finish';
                    str2 = 'Add more';
                    default = str2;
                    button = questdlg(qstring,wintitle,str1,str2,default);
                    if strcmp(button,str1)
                        addflag = 0;
                    end
                    
                end
                
            end
            
        catch
            disp('No ROIs added to object.')
        end
        
        
        
        
    case 1
        
        try
            
            mFig = figure('Color',[0.3 0.3 0.3],'SizeChangedFcn',@resetAxes,'Position',get(0, 'Screensize'));
            subplot(2,2,3)
            avgAx = gca;
            imagesc(f)
            axis off
            
            
            subplot(2,2,1)
            maskAx = gca;
            axis off
            maskAx.Color = [0 0 0];
            
            subplot(2,2,[2])
            mAx = gca;
            mAx.Color = [0.6 0.6 0.6];
            mAx.XAxis.Color = [1 1 1];
            mAx.YAxis.Color = [1 1 1];
            
            hold on
             subplot(2,2,4)
        dotAx = gca;    
        hold on
        dotAx.Color = [0.8 0.8 0.8];

            % Add listener to keep images at max size
            % addlistener(mAx, 'Resize', @(obj,event)resetAxes(maskAx,avgAx));
            
            ridx = 0;
            u = unique(polAngImg);
            r = zeros( size(f,1),size(f,2),length(u) );
            masks = zeros( size(f,1),size(f,2));
            masksbw = zeros( size(f,1),size(f,2));
            maskscale = masksbw;
            f=double(f);
            for l = 1:length(u)
                rk_scale =  mean(mean(f(polAngImg==u(l))))./max(f(:));
                if rk_scale*max(f(:)) > mean(f(:)) +std(f(:))
                    ridx = ridx+1;
                    rk = polAngImg==u(l);
                    col = ROIcolor(ridx);
                    maskcol = cat(3,col(1)*rk,col(2)*rk,col(3)*rk);
                    masks = masks + maskcol;
                    masksbw = u(l).*rk + masksbw;
                    maskscale = rk_scale.*rk + maskscale;
                    r(:,:,ridx) = rk;
                     numHalfCycles = 2*length(find(obj.TrialPatNum(obj.TrialSeqNum == expNum) == 360));
                maskResponse = scanROI(obj,rk,[polStart, polStop]);
                CycleLims = floor(linspace(1,length(maskResponse),numHalfCycles+1));
              maskTrialResp = nan(numHalfCycles,max(diff(CycleLims)));
                for cycIdx = 1:numHalfCycles
                    thisResp = maskResponse(CycleLims(cycIdx):CycleLims(cycIdx+1)-1);
                   % Subtract off mean:
                   maskTrialResp(cycIdx,1:CycleLims(cycIdx+1)-CycleLims(cycIdx)) = thisResp - nanmean(thisResp); 
             

                end
               angleVector = linspace(30,180,size(maskTrialResp,2));
%                 plot(mAx,angleVector,nanmean(maskTrialResp),'Color',col,'LineWidth',3)
%                 hold on
%                 plot(mAx,angleVector,(maskTrialResp)','Color',col,'LineWidth',1)
                lineProps.col = {col};
                mAx;
                mseb(angleVector,nanmean(maskTrialResp),nanstd(maskTrialResp),lineProps,1);
                    
                    autoROI(ridx).mask = rk;
                    autoROI(ridx).color = [1 1 1];
                    autoROI(ridx).position = [1 1 1 1];
                    autoROI(ridx).response = [];             
                end
                
            end
            
            imshow(masks,'Parent',maskAx)
            bht = bentitle(obj.File,'none');
            bht.Color = 'w';
            
            ylim = mAx.YLim(2);
            if ~isempty(obj.TrialStartFrame)
                for t = 1:length(obj.TrialStartFrame)
                    plotStimPatch(obj,mAx,[obj.TrialStartFrame(t) obj.TrialEndFrame(t)])
                    xpos = obj.TrialStartFrame(t) + 0.5*(obj.TrialEndFrame(t) - obj.TrialStartFrame(t));
                    if ~isempty(obj.TrialPatNum)
                        text(mAx,xpos,ylim,[num2str(obj.TrialPatNum(t))],'HorizontalAlignment','center');
                    end
                end
            end
            
%     filename = [obj.File];
%     export_fig( fullfile(savedir,filename), ...
%         '-png', '-m4', gcf);
%         filename2 = [obj.File ' autoROIs'];
%      export_fig( fullfile(savedir2,filename2), ...
%         '-png', '-m4', gcf);
%     
%     pause(5)
%     close gcf
%             
        catch
            disp('No ROIs added to object.')
            
        end
        
        
        
end


if ~isempty( [autoROI.mask] )
    if nargout
        varargout{1} = autoROI;
    else
        
        % assign to object
        objROIidx = length(obj.ROI);
        if objROIidx == 0
            obj.ROI = struct('mask',[],'color',[],'position',[],'response',[]);
        end
        for roiIdx = 1:length(autoROI)
            objROIidx = objROIidx + 1;
            obj.ROI(objROIidx) = autoROI(roiIdx);
            obj.ROI(objROIidx).color = ROIcolor(objROIidx);
            disp(['ROI assigned to object as ROI index ' ROIcolor(objROIidx,1) ])
        end
        
    end
end

% delete(mFig)
end

function angleFrameArray = getAngleLogicArray(obj,polAng)

% Extract the start:stop frame numbers for the pol tuning experiment:

if any(obj.TrialSeqNum==4)
    expNum = 4;
elseif any(obj.TrialSeqNum==2)
    expNum = 2;
else
    return
    error('No pol tuning exp (2 or 4) available')
end

polStart = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first'));
polStop = obj.TrialEndFrame(find(obj.TrialSeqNum == expNum,1,'last'));


% Find the frames recorded with polarizer at polAng:
angStarts = obj.TrialStartFrame(obj.TrialPatNum == polAng);
angStops = obj.TrialEndFrame(obj.TrialPatNum == polAng);

% Add the equivalent +180deg angle
angStarts = [angStarts, obj.TrialStartFrame(obj.TrialPatNum == polAng + 180)];
angStops = [angStops, obj.TrialEndFrame(obj.TrialPatNum == polAng + 180)];

% Construct logical array for frameAngles == polAng 
angleFrameArray = zeros(1,size(obj.Frames,3));
for n = 1:length(angStarts)
    if (angStarts(n) >= polStart) && (angStops(n) <= polStop)
        angleFrameArray(angStarts(n):angStops(n)) = 1;
    end
end


end

function resetAxes(src,~)

for n = 1:length( src.Children )
    if  strcmp( src.Children(n).Visible, 'off' )
        g = src.Children(n);
        axis(g,'equal');
        g = [];
    end
end
end

function connMask = smoothPolAngImg(polAngImg)

polAngImgCopy = polAngImg;
connMask = zeros(size(polAngImg));

minArea = 16;

uniqueAng = unique(polAngImg(polAngImg>0));
for n = 1:length(uniqueAng)
           I = polAngImgCopy == uniqueAng(n);
           connMask = connMask +  uniqueAng(n).*bwareaopen(I,minArea);
end

smallAreas = polAngImgCopy.*(polAngImgCopy & ~connMask);
plusminus1Ang(:,1) = circshift(uniqueAng,1);
plusminus1Ang(:,2) = circshift(uniqueAng,-1);
for n = 1:length(uniqueAng)
           I = smallAreas == uniqueAng(n);
           rp = regionprops(I,'PixelIdxList');
                          
           % Examine each small area individually
           for rIdx = 1:length(rp)
               % Get the surrounding perimeter of the area 
               areaPix = zeros(size(polAngImg));
               areaPix(rp(rIdx).PixelIdxList) = 1;
               perimPix = find(imdilate(areaPix, true(3)) - areaPix);
               perimPixVals = polAngImgCopy(perimPix);
               
               % If all surrounding pixels are 0, it's a small area of
               % noise, too small to be a synapse: set to zero.
               if mode(perimPixVals) == 0                 
                   connMask(rp(rIdx).PixelIdxList) = 0;
                   
                   % If non-zero surrounding pixels are within one
                   % angle-step, they are small noisy areas likely to
                   % belong to the surrounding area or vice versa:                  
               elseif any( mode(perimPixVals(perimPixVals>0)) == plusminus1Ang(n,:))
                   % Find area which is larger 
                   innerSize = length(rp(rIdx).PixelIdxList);                   
                   [arow,acol] = ind2sub(size(polAngImg),perimPix( perimPixVals==mode(perimPixVals(perimPixVals>0)) ));
                   outerPix = find(grayconnected(polAngImgCopy,acol(1),arow(1)));
                   outerSize = length(outerPix);
                   
                   if outerSize > innerSize
                       connMask(rp(rIdx).PixelIdxList) = mode(perimPixVals(perimPixVals>0));
                   else
                       connMask(outerPix) = uniqueAng(n);
                   end
                   % If non-zero pixels are 
               else
                   % If non-zero surrounding pixels are NOT within one
                   % angle-step, they are small noisy areas which we will
                   % not include in a mask:
                   connMask(rp(rIdx).PixelIdxList) = 0;
               end
           end
end
end