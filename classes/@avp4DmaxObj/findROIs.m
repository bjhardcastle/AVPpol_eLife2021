function varargout = findROIs(obj, numROIs, expNum, detrend, saveplot)
if nargin < 4
    detrend = 1;
end
if nargin < 2
    numROIs = 5;
end
if numROIs > 16
    numROIs = 16;
    disp('Using maximum numROIs: 16')
end
if isempty(obj.Frames)
    getFrames(obj)
end


if nargin <3 || ~any(obj.TrialSeqNum==expNum)
    expNum = [];
end
if isempty( expNum)
    f = median(obj.Frames,3);
else
    f = median(obj.Frames(:,:,obj.expStart(expNum):obj.expStop(expNum)),3);
end
if ~isempty(obj.polExp)
    s = obj.polSelImg;
end
if nargin < 5 || isempty(saveplot)
    
    saveplot = 0;
else
    savedir2 = obj.Folder;
    
end

[row2,col2] = find(f>-9999); % ALl pixels
if ~isempty(obj.layerMask)
    [row1,col1] = find(obj.layerMask.mask); % Above threshold pixels
    
%     if expNum== obj.polExp
%             [row1,col1] = find(~isnan(obj.polPix)); % Above threshold pixels
% end
else
    [row1,col1] = find(f>median(f(:))); % Above threshold pixels
end


inds2 = sub2ind ( size(f), row2, col2 );
inds1 = sub2ind ( size(f), row1, col1 );
inds = setdiff(inds2,inds1);
inds = inds1;
% Detrend Frames:
if detrend
    detrendedFrames = detrendFrames(obj);
else
    detrendedFrames = obj.Frames;
end
if ~isempty(expNum)
    scanlims = [obj.expStart(expNum) obj.expStop(expNum)];
        detrendedFrames = detrendedFrames(:,:,scanlims(1):scanlims(2));
else
scanlims = [];
end
resp = reshape(detrendedFrames, length(f(:)),size(detrendedFrames,3));
clear detrendedFrames
kblank = zeros(size(resp,1),1);
resp = resp(inds1,:);

% Run kmeans on raw time-series
rng default 

myfunc = @(X,K)(kmeans(X, K, 'MaxIter',1000, 'replicate',5));
eva = evalclusters(resp,myfunc,'CalinskiHarabasz',...
    'klist',[3:9])
rng default 
kidx = kmeans(resp,eva.OptimalK+1, 'MaxIter',1000,'Replicates',5);
kblank(inds1) = kidx;
% Reshape vector into image array for displaying masks
k = reshape(kblank,size(f));

autoROI = struct('mask',[],'color',[],'position',[],'response',[]);

switch saveplot
    case 0
        
        mFig = figure('Color',[0.5 0.5 0.5],'SizeChangedFcn',@resetAxes);
        subplot(2,2,3)
        avgAx = gca;
           imagesc( obj.polPix )
   colormap([1 1 1; flipud(hsv(180))])
set(avgAx,'cLim',[0 180])
axis(avgAx,'equal');
        axis off
        
        
        subplot(2,2,1)
        maskAx = gca;
        axis off
        maskAx.Color = [0 0 0];
        
        subplot(2,2,[2,4])
        mAx = gca;
        hold on
        
        
        % Add listener to keep images at max size
        % addlistener(mAx, 'Resize', @(obj,event)resetAxes(maskAx,avgAx));
        
        ridx = 0;
        u = unique(k(k>0));
        r = zeros( size(f,1),size(f,2),length(u) );
        masks = zeros( size(f,1),size(f,2));
        masksbw = zeros( size(f,1),size(f,2));
        maskscale = masksbw;
        f=double(f);
        for l = 1:length(u)
            rk_scale =  mean(mean(f(k==u(l))))./max(f(:));
            if  expNum== obj.polExp
                if sum(~isnan(obj.polPix(k==u(l)))) > sum(isnan(obj.polPix(k==u(l))))
                    %mean(s(obj.bkgMask)) + std(s(obj.bkgMask))
                    selFlag = 1;
                else
                    selFlag =0;
                end
            elseif  rk_scale*max(f(:)) > mean(f(:)) +std(f(:))
                selFlag = 1;
            else
                selFlag = 0;
            end
            if selFlag
                ridx = ridx+1;
                rk = k==u(l);
                col = ROIcolor(ridx);
                maskcol = cat(3,col(1)*rk,col(2)*rk,col(3)*rk);
                masks = masks + maskcol;
                masksbw = u(l).*rk + masksbw;
                maskscale = rk_scale.*rk + maskscale;
                r(:,:,ridx) = rk;
                plot(mAx,scanROI(obj,rk,scanlims),'Color',col)
                hold on
            end
            
        end
        ylim = mAx.YLim(2);
        if ~isempty(obj.TrialStartFrame)
            for t = 1:length(obj.TrialStartFrame)
%                 plotStimPatch(obj,mAx,[obj.TrialStartFrame(t) obj.TrialEndFrame(t)])
                xpos = obj.TrialStartFrame(t) + 0.5*(obj.TrialEndFrame(t) - obj.TrialStartFrame(t));
                if ~isempty(obj.TrialPatNum)
                    text(mAx,xpos,ylim,[num2str(obj.TrialPatNum(t))],'HorizontalAlignment','center');
                end
            end
        end
        
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
               imagesc( obj.polPix )
   colormap([1 1 1; flipud(hsv(180))])
set(gca,'cLim',[0 180])
axis(avgAx,'equal');
            axis off
            
            
            subplot(2,2,1)
            maskAx = gca;
            axis off
            maskAx.Color = [0 0 0];
            
            subplot(2,2,[2,4])
            mAx = gca;
            mAx.Color = [0.6 0.6 0.6];
            mAx.XAxis.Color = [1 1 1];
            mAx.YAxis.Color = [1 1 1];
            
            hold on
            
            
            % Add listener to keep images at max size
            % addlistener(mAx, 'Resize', @(obj,event)resetAxes(maskAx,avgAx));
            
            ridx = 0;
            u = unique(k(k>0));
            r = zeros( size(f,1),size(f,2),length(u) );
            masks = zeros( size(f,1),size(f,2));
            masksbw = zeros( size(f,1),size(f,2));
            maskscale = masksbw;
            f=double(f);
            figure,imagesc(k)
            for l = 1:length(u)
                rk_scale =  mean(mean(f(k==u(l))))./max(f(:));
                if  expNum== obj.polExp
                    if mean(s(k==u(l)))>mean(s(obj.bkgMask)) + std(s(obj.bkgMask))
                        selFlag = 1;
                    else
                        selFlag =0;
                    end
                elseif  rk_scale*max(f(:)) > mean(f(:)) +std(f(:))
                    selFlag = 1;
                else
                    selFlag = 0;
                end
                if selFlag
                    ridx = ridx+1;
                    rk = k==u(l);
                    col = ROIcolor(ridx);
                    maskcol = cat(3,col(1)*rk,col(2)*rk,col(3)*rk);
                    masks = masks + maskcol;
                    masksbw = u(l).*rk + masksbw;
                    maskscale = rk_scale.*rk + maskscale;
                    r(:,:,ridx) = rk;
                    plot(mAx,scanROI(obj,rk,scanlims),'Color',col,'LineWidth',2)
                    hold on
                    
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

                filename2 = [obj.File ' autoROIs'];
 export_fig( fullfile(savedir2,filename2), ...
        '-png', '-m4', gcf);
    
    pause(3)
    close gcf
                     filename = [obj.File '_patches'];
export_fig( fullfile(savedir2,filename), ...
        '-png', '-m4', gcf);
close gcf
    pause(3)
       
        catch
            disp('No ROIs added to object.')
            
        end
        
        
        
end


if ~isempty( [autoROI.mask] )
    if nargout
        varargout{1} = autoROI;
    else
        
        % assign to object
            obj.ROI = struct('mask',[],'color',[],'position',[],'response',[]);
objROIidx = 0;
for roiIdx = 1:length(autoROI)
    objROIidx = objROIidx + 1;
    obj.ROI(objROIidx) = autoROI(roiIdx);
    obj.ROI(objROIidx).color = ROIcolor(objROIidx);
    disp(['ROI assigned to object as ROI index ' ROIcolor(objROIidx,1) ])
end

    end
end

delete(mFig)


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


