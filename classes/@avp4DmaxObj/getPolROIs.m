function varargout = getPolROIs(obj, numROIs, minPixArea)
if nargin< 3 || isempty(minPixArea)
    minPixArea = 20;
end
detrendFlag = 0;

if any(obj.TrialSeqNum==4)
    expNum = 4;
elseif any(obj.TrialSeqNum==2)
    expNum = 2;
else
    disp('No pol tuning exp (2 or 4) available')
    return
end

if nargin< 2|| isempty(numROIs)
    numROIs = 32;     % maximum number of areas to be combined, per angle presented
elseif numROIs > 32
    numROIs = 32;
    disp('Using maximum numROIs: 32')
end

% For this experiment (whether pol attached or not) find start/stop frames:
polStart = obj.expStart(expNum);
polStop = obj.expStop(expNum);
polAngs = obj.expPolAngs(expNum);

discreteTuning = obj.polPixDiscrete;
continuousTuning = obj.polPix;

% check pol pix exist / mask exist and larger than 10pix area
if length(~isnan(discreteTuning(:)))< 10
    disp('Less than 10 pol pixels - skipped')
    return
end

loadLayerMasks(obj)
polColMap = flipud(hsv(180));

if ~isempty(obj.layerMask)
    maskFrame = obj.layerMask.mask;
else
    maskFrame =  obj.fftMagImg > 0.2 & obj.polSelImg > 0.2;
end

% For each angle:
u = unique(discreteTuning(~isnan(discreteTuning))); % angles in the tuning map
ridx = 0;
polROI = [];
for l = 1:length(u+1) % extra iteration to find some non-responsive control ROI
    
    
    if l < length(u+1)
        
        controlROIidx = 0;
        I = logical((discreteTuning == u(l)).*maskFrame);
        
    elseif l == length(u+1)
        controlROIidx = 1;
        I = logical((obj.polSelImg<0.1).*maskFrame);
    end
    % detect blobs
    
    rp = regionprops(I,'PixelIdxList','Area','Centroid');
    % threshold size
    rp([rp.Area]< minPixArea) = [];
    if ~isempty(rp)
        
        % examine variance of each blob's pixels in continuous tuning map
        for rpidx = 1:length(rp)
            rp(rpidx).Std = std(continuousTuning(rp(rpidx).PixelIdxList));
            %              rp(rpidx).Std = mean(obj.avgPolImg(rp(rpidx).PixelIdxList));
        end
        
        
        % sort and keep lowest variance for each angle
        [~,sortIdx]= sort([rp.Std]','ascend');
        % sort and keep largest area for each angle
                [~,sortIdx]= sort([rp.Area]','descend');

        rp = rp(sortIdx);
        
        rpidx = 0;
        while rpidx < numROIs && rpidx<length(rp)
            rpidx = rpidx + 1;
            newMask = zeros(size(I));
            newMask(rp(rpidx).PixelIdxList) = 1;
            
            % Fill holes in mask (just for tracing outline - won't be saved)
            bw = imfill(newMask,'holes');
            % Trace outline of filled mask
            B = bwboundaries(bw);
            if ~isempty(B)
                roix = [];
                roiy = [];
                for b = 1:length(B)
                    roix = [roix B{b}(:,2)'];
                    roiy = [roiy B{b}(:,1)'];
                end
                % Downsample the outline by skipping every nth point
                skipPts = 3;
                roiPos = [ roix([1:skipPts:end])' roiy([1:skipPts:end])'];
            end
            
            
            ridx = ridx+1;
            
            % Get the pixels included at this angle
            polROI(ridx).mask = logical(bw);
            polROI(ridx).position = roiPos;
            polROI(ridx).response = scanROI(obj,newMask)';
            polROI(ridx).centroidX = rp(rpidx).Centroid(1);
            polROI(ridx).centroidY = rp(rpidx).Centroid(2);
            if ~controlROIidx
                polROI(ridx).angle = u(l);
                polROI(ridx).color = polColMap(u(l),:);
            else
                polROI(ridx).angle = nan;
                polROI(ridx).color = ROIcolor(8);
            end
        end
    end
    
end


% assign to object as polROIs

if ~isempty(polROI)
%     obj.polROI = [];
    obj.polROI = polROI;
else
    disp('No polROIs made')
end


%{
% extract responses
%  make dependent value for (contains polROIs ) = 0 or 1


f = median(obj.Frames(:,:,polStart:polStop),3);
g = mean(obj.Frames(:,:,polStart:polStop),3);

% Detrend Frames:
if detrendFlag
    detrendedFrames = detrendFrames(obj);
else
    detrendedFrames = obj.Frames;
end

% Get indices for frame size
indLength = size(obj.Frames,1)*size(obj.Frames,2);
rows = reshape(repmat(1:size(obj.Frames,1),1,size(obj.Frames,2)),indLength,1);
cols = reshape(repmat(1:size(obj.Frames,2),size(obj.Frames,1),1),indLength,1);
inds = sub2ind(size(obj.Frames(:,:,1)), rows, cols);

% Extract every pixel's time-series (whole experiment)
resp = reshape(detrendedFrames, indLength, size(obj.Frames,3));
detrendedFrames = [];

% Extract every pixel's time-series (pol tuning experiment only)
% Detrend Frames:
if detrendFlag
    respPolExp = detrend(resp(inds,polStart:polStop)')';
else
    respPolExp = resp(inds,polStart:polStop);
end


%% FFT stuff
% This is mainly to obtain an array of preferred tunings for every pixel,
% without binning responses by angle (yet)
% Magnitude data can also be used later to find pixels that respond most to
% the stimulus

% Find FFT of every pixel's time-series:
L = size(respPolExp,2);
n = 2^nextpow2(L); % sampling of FFT
Y = fft(respPolExp,n,2); % FFT array [pixel index , frequency index]
Fs = obj.AIrate/obj.IFI;
F = Fs*(0:(n/2))/n;
P = abs(Y/n);

% plot frequency spectrum from FFT:
%{
figure
plot(F,P(:,1:n/2+1))
%}

% Now work out the frequency of responses caused by polarizer:
% (one cycle of polarizer rotation = two cycles of activity in
% pol-sensitive cells)
numHalfCycles = 2*length(find(obj.TrialPatNum(obj.TrialSeqNum == expNum) == 360));
halfCycFreq = 1/((1/numHalfCycles)*(polStop-polStart+1)/Fs);

% Find corresponding FFT frequency:
[~,fIdx]=min(abs(halfCycFreq - F));

% Array of fft magnitudes for each pixel:
% fftMagImg = reshape((abs(Y(:,fIdx))),size(f));
fftMagImg = reshape(sqrt(abs(Y(:,fIdx))),size(f));

% Find the pref. pol angle for each pixel:
%   +ve phase of FFT means shifted to 'earlier' angle, which is
%   counterintuitive when angle usually increases with time, so mulitply by -1
%   The pol stimulus might not start at zero, so add the first angle
%   (before division by 2, so must be doubled)
%   Then convert -90:90 to 0:180 with wrapTo360
%   Phase angle from FFT is in range -pi:pi, which maps to -90:90 of the
%   pol stimulus, so we must divide phase angle by 2
angPref = wrapTo360(polAngs(1)*2 + rad2deg(-angle(Y(:,fIdx))))./2;

fftAngImg = reshape(angPref,size(f));

%{
% Run kmeans
kidx = kmeans(angPref,numROIs);
kblank = nan(size(resp,1),1);
kblank = kidx;

% Reshape vector into image array for displaying masks
k = reshape(kblank,size(f));
%}

%{
% old code made during testing:
%  maskbwAng = kAng;
%  maskbwAngVector = reshape(maskbwAng, indLength, 1);
% maskbwAngVectorData = maskbwAngVector.*ktempMask(:);
%  clear maskParams
%  maskParams(:,1) = sind(2.*maskbwAngVectorData);
% Project angle data to two components:
%  maskParams(:,2) = cosd(2.*maskbwAngVectorData);
% Get scaled (x,y) pixel positions, to enter into k-means clustering (so
% clustering is on magnitude of other parameters AND spatial location):
% rowData = rows( :);
% rowScaled = 2.*(rowData - min(rowData(:)))./ (max(rowData(:)) - min(rowData(:)) );
%  maskParams(:,3) =0.5.*rowScaled;
% colData = cols( :);
% colScaled = 2.*(colData - min(colData(:)))./ (max(colData(:)) - min(colData(:)) );
%  maskParams(:,4) = 0.5.*colScaled;

% maskParams(:,5) = polSelVec;

% Find optimum number of clusters:
% eva = evalclusters(resp.*ktempMask(:),'kmeans','CalinskiHarabasz','KList',[1:10]);

% Run kmeans
%  kblank = nan(size(maskbwAngVector,1),1);
% kidx = kmeans(respExp4(rankedMagIdx,:),eva.OptimalK);
% kblank(rankedMagIdx) = kidx;
%}


%% Binned response stuff
% Here we get the average response of every pixel for each angle presented.
% From this we can find the pol selectivity: (max - min)/max
% which is a good parameter for excl./including pixels


% Construct array [numPix x numPolAngles/2], into which we'll enter the
% average response for each pol angle, for every pixel:
pixPolResp = nan(indLength,0.5*length(polAngs));

%{
% For each angle:
for aIdx = 1:0.5*length(polAngs)
    % Get logic array that indicates every frame that was recorded at the
    % angle (including angle+180)
    angleFrameArray = getAngleLogicArray(obj,polAngs(aIdx));
    
    % Extract mean response for this angle:
    pixPolResp(:,aIdx) = mean(respPolExp(:,angleFrameArray(polStart:polStop)),2);
end
% Find maximum response value for each pixel:
% [respMax,maxIdx] = nanmax(pixPolResp,[],2);

%}

% For each angle:
for aIdx = 1:length(polAngs)
    % Get logic array that indicates every frame that was recorded at the
    % angle (including angle+180)
    angleFrameArray = getAngleLogicArray2(obj,polAngs(aIdx));
    
    % Extract mean response for this angle:
    pixPolResp(:,aIdx) = mean(respPolExp(:,angleFrameArray(polStart:polStop)),2);
end

[r,axtheta] = circ_rangle(circ_axial(deg2rad(repmat(polAngs(1:length(polAngs)),size(pixPolResp,1),1)),2), pixPolResp, deg2rad(mode(diff(polAngs))), 2);
% [r,theta] = circ_rangle(deg2rad(repmat(polAngs(1:0.5*length(polAngs)),size(pixPolResp,1),1)), pixPolResp, deg2rad(mode(diff(polAngs))), 2);
theta = axtheta/2;
% convert angle into linear orientation index [ 1: 0.5*length(polAngs) ]

% If polAngs = [30 60 90], make bin edges [15 45 75 105]
angEdges = polAngs(1:0.5*length(polAngs)+1)-0.5*(mode(diff(polAngs)))  ;

% deal with 0/180:
Ang =0.5.*(wrapTo360((2.*(rad2deg(theta)))));
shiftAng =0.5.*(wrapTo360((2.*(rad2deg(theta)-15))))  + 15;

maxIdx = discretize(shiftAng,angEdges);  % fft mag threshold applied




% This section allows for discretization of the continuous tuning map with
% any bin size:
binWidth = 30;
idxFreq = mode(diff(polAngs))/binWidth; % Only examine these indices out of 'discretize'
if idxFreq~=1
    allBins = polAngs(1)-0.5*binWidth:binWidth:360;
    lastBin = find(allBins ==(polAngs(length(polAngs)*0.5))+allBins(1)-binWidth);
    angEdges = allBins(1:lastBin);
    Ang =0.5.*(wrapTo360((2.*(rad2deg(theta)))));
    shiftAng =0.5.*(wrapTo360((2.*(rad2deg(theta)- angEdges(1)))))  + angEdges(1);
    maxIdxNarrow = discretize(shiftAng,angEdges);  % fft mag threshold applied
    
    maxIdxNarrow(mod(maxIdxNarrow,idxFreq)~=1) = nan;
    maxIdxNarrow = floor(maxIdxNarrow/idxFreq)+1;
else
    maxIdxNarrow = maxIdx;
end








% For each angle:
pixPolResp = [];
for aIdx = 1:0.5*length(polAngs)
    % Get logic array that indicates every frame that was recorded at the
    % angle (including angle+180)
    angleFrameArray = getAngleLogicArray(obj,polAngs(aIdx));
    
    % Extract mean response for this angle:
    pixPolResp(:,aIdx) = mean(respPolExp(:,angleFrameArray(polStart:polStop)),2);
end
respMax = pixPolResp( sub2ind(size(pixPolResp), [1:indLength]', maxIdx ) );

% From the max angle, find the orthogonal angle, which will be used as the
% minimum response:
maxIdx = maxIdx - 1; % To allow 'mod' function below to work. Add 1 again after using.
respMin = pixPolResp( sub2ind(size(pixPolResp), [1:indLength]', 1+mod(maxIdx + 0.25*length(polAngs),0.5*length(polAngs)) ) );
maxIdx = maxIdx + 1;

% Get the pol selectivity vector:
polSelVec = ( respMax - respMin )./respMax;
% and image:
polSelImg = reshape(polSelVec,size(f));
polSelImg = polSelImg.*(f>20); % threshold just above 0 to get rid of empty border from registration

% We can also find the preferred angle from the binned average responses:
discreteTuning = reshape(maxIdx,size(f)).*mode(diff(polAngs));
discreteTuning = discreteTuning.*(polSelImg>obj.polSelThreshold);






% discreteTuning = reshape(maxIdxNarrow,size(f)).*mode(diff(polAngs));
% discreteTuning = discreteTuning.*(polSelImg>obj.polSelThreshold);








polColMap = obj.polColMap;
% polColMap = flipud(hsv(0.5*length(polAngs)));
% polColMap(1+0.5*length(polAngs),:) = ROIcolor(8); % Gray for control ROI in last place
%{
% Polar plot of polSel vs angPref
figure
polarplot(-(angle(Y(:,fIdx)))./2,polSelVec(:),'.','MarkerSize',1)
set(gca,'ThetaLim',[-90 90],'ThetaDir','clockwise','ThetaZeroLocation','right')
%}
% connMask = smoothdiscreteTuning(discreteTuning);


discreteTuning(discreteTuning == 0) = nan;
k = discreteTuning./mode(diff(polAngs));


figure,
subplot(1,2,1)
respOrientImg = nan(size(f));
respOrientImg(polSelVec>obj.polSelThreshold) = Ang(polSelVec>obj.polSelThreshold);
respOrientImg(isnan(discreteTuning)) = nan;
imagesc(reshape(respOrientImg,size(f)))
colormap([[1 1 1];flipud(hsv)])
set(gca,'CLim',[0 180])




%% Joint array construction
% Make a single array, using various thresholds to discard non-responsive
% pixels or pixels we're not interested in, so that we can detect areas
% with common properties (preferred orientation)

polSelThreshold = obj.polSelThreshold;   % polSel = (maxresp - minresp)/maxresp   only consider pixels with polSel> Threshold
polRespMagThreshold = obj.polMagThreshold;  %  top X percentage: 10 would be top 10%

% Using the magnitude of the FFT, find the indices of the top xx% for
% thresholding:
[~,sortIdx] = sort(P(:,fIdx),'descend');
rankedMagIdx = sortIdx(1:ceil(polRespMagThreshold*length(sortIdx)/100)); % top 10 per cent, when ranked by FFT magnitude at pol freq
rankedMagMask = nan(size(f));
rankedMagMask(rankedMagIdx) = 1;

% Discretize that FFT-derived preferred angle data, by binning according to
% angles presented:
% If polAngs = [30 60 90], make bin edges [15 45 75 105]
angEdges = polAngs(1:0.5*length(polAngs)+1)-0.5*(mode(diff(polAngs)))  ;
% deal with 0/180:
AngF =0.5.*(wrapTo360((2.*(fftAngImg))));

shiftAng =0.5.*(wrapTo360((2.*(fftAngImg -15))))  + 15;

kd = discretize(shiftAng.*rankedMagMask,angEdges);  % fft mag threshold applied

k = k.*rankedMagMask;
%{
% Threshold binned angle data, based on pol selectivity
k = nan(size(angPref,1),1);
k(polSelVec>polSelThreshold) = kd(polSelVec>polSelThreshold); % binned pol sel threshold applied
k = reshape(k,size(f));
% k is an array of nans, size of obj.Frames(:,:,1),
% pixels which respond above threshold, with max resp to polAngs(1) have value 1
% pixels which respond above threshold, with max resp to polAngs(2) have value 2
% etc.

% assign k too: convert to angles first
obj.k= k;
%}
obj.k= respOrientImg;
obj.k = obj.k.*(polSelImg>polSelThreshold).*rankedMagMask;
obj.k(obj.k==0) =  nan;
obj.k(isnan(discreteTuning)) = nan;

subplot(1,2,2)
fftOrientImg = nan(size(f));
fftOrientImg(polSelVec>polSelThreshold) = k(polSelVec>obj.polSelThreshold);
imagesc(obj.k)
colormap([[1 1 1];flipud(hsv)])
set(gca,'CLim',[0 180])

% subplot(1,2,2)
% fftOrientImg = nan(size(f));
% fftOrientImg(polSelVec>obj.polSelThreshold) = Ang(polSelVec>obj.polSelThreshold);
% imagesc(reshape(fftOrientImg,size(f)))
% colormap([[1 1 1];flipud(hsv)])
% set(gca,'CLim',[0 180])
%% Find 'control areas', below thresholds
u = unique(k(~isnan(k)));
rS  = (reshape(f(:).*(polSelVec<polSelThreshold),size(f)));
rSM = rS.*isnan(rankedMagMask);
rSMm = rSM>mean(obj.Frames(:))+std(obj.Frames(:));
%{
% figure,imagesc(rSMm)

% if any(rSMm(:))
%     controlROI = struct('mask',[],'color',[],'position',[],'response',[]);
%
%                     ck = rSMm;
%                 cp = regionprops(rSMm,'PixelIdxList','Area');
%                 for cpIdx = 1:length(cp)
%                     % Get the pixels within each area and sum their avg
%                     % intensity values. Then use it to find a single ROI area
%                     % which is large and bright
%                     cp(cpIdx).Sum = sum(f([cp(cpIdx).PixelIdxList]));
%                 end
%                 [~,cpSrtIdx] = sort([cp.Sum],'descend');
%
%                 % Sort by area size and for the top X areas, get the idx of
%                 % all pixels contained within:
%                 cp = cp(cpSrtIdx);
%
%                 % Insert pixels into empty array to make mask
%                 controlMask = zeros(size(f));
%                 controlMask(cp(1).PixelIdxList) = 1;
%
%       controlROI.mask = controlMask;
%       controlROI.color = ROIcolor(8);
%       controlROI.position = [1 1 1 1];
%       controlROI.response = scanROI(obj,controlMask);
% end
%}
% All pixels found will be added to an (and if requested, a single
% subset area will be identified below, as an example )
% Add control areas as a new 'angle' value to k (will be highest numerical
% value)
controlROIidx = u(end)+1;
k(rSMm) = controlROIidx;

%% Assign useful arrays to object:

obj.discreteTuning = discreteTuning;
obj.polSelImg = polSelImg;

obj.polTuningImg = reshape(Ang,size(f));
obj.fftAngImg = fftAngImg;
obj.fftMagImg = fftMagImg;
obj.avgPolImg = f;
obj.varImg = std(obj.Frames,[],3).*mean(obj.Frames,3);
obj.polActivityImg = g;
%% Save areas with all pixels first
u = unique(k(~isnan(k)));
ridx = 0;
for l = 1:length(u)
    
    ridx = ridx+1;
    
    % Get the pixels included at this angle
    polROI(ridx).mask = k==u(l);
    polROI(ridx).color = polColMap(u(l),:);
    polROI(ridx).position = [1 1 1 1];
    polROI(ridx).response = scanROI(obj,[k==u(l)]);
    if u(l) ~= controlROIidx
        polROI(ridx).angle = polAngs(u(l));
    else
        polROI(ridx).angle = nan;
    end
end

obj.fullROI = polROI;

%% Assign a number of example ROIs automatically

autoROI = struct('angle',[],'mask',[],'color',[],'position',[],'response',[]);


u = unique(k(~isnan(k)));
r = zeros( size(f,1),size(f,2),length(u) );
% Invert f for displaying
fInv =  (1-(f./max(f(:))));
% Create ~mask, setting pixels that will become colored, polSel
% regions to zero (since we will add values for color)
fInvMask = fInv.*(isnan(k));
fullMasks = cat(3,fInvMask,fInvMask,fInvMask);
examplemasks = cat(3,fInvMask,fInvMask,fInvMask);
%     masksbw = zeros( size(f,1),size(f,2));
%     maskscale = masksbw;
f=double(f);

ridx = 0;
for l = 1:length(u)
    
    ridx = ridx+1;
    
    col = polColMap(u(l),:);
    
    rk = k==u(l);
    maskcol = squeeze(cat(3,col(1)*rk,col(2)*rk,col(3)*rk));
    fullMasks = fullMasks + maskcol;
    
    % Get the pixels included at this angle
    rk = k==u(l);
    rp = regionprops(rk,'PixelIdxList','Area');
    [~,rpSrtIdx] = sort([rp.Area],'descend');
    % Sort by area size and for the top X areas, get the idx of
    % all pixels contained within:
    rp = rp(rpSrtIdx);
    % Group the top numROIs areas by size:
    %                 %{
    if length(rp) < numROIs
        pixForColor = vertcat(rp(1:end).PixelIdxList);
    else
        pixForColor = vertcat(rp(1:numROIs).PixelIdxList);
    end
    %                 %}
    % -- changed to take the top area, only:
    %                 pixForColor = vertcat(rp(1).PixelIdxList);
    
    % Get the idx of all pixels currently set to zero in
    % 'masks', awaiting color (because not all will be colored
    % now)
    pixAtZero = find(k==u(l));
    % Find the pixels which won't be colored and should be
    % reset to original average image (inverted)
    pixToReset = setdiff(pixAtZero,pixForColor);
    
    % Update polAng array to discard all areas but the top X
    k(pixForColor) = u(l);
    k(pixToReset) = nan;
    
    % Update masks array to discard all areas but the top X
    tempMask = zeros( size(f,1),size(f,2));
    tempMask(pixToReset) = 1;
    toAdd = tempMask.*fInv;
    toAdd3 = cat(3,toAdd,toAdd,toAdd);
    examplemasks = examplemasks + toAdd3;
    
    rk = k==u(l);
    maskcol = squeeze(cat(3,col(1)*rk,col(2)*rk,col(3)*rk));
    examplemasks = examplemasks + maskcol;
    
    maskResponse = scanROI(obj,rk);
    
    
    autoROI(ridx).mask = rk;
    autoROI(ridx).color = col;
    autoROI(ridx).position = [1 1 1 1];
    autoROI(ridx).response = maskResponse;
    if u(l) ~= controlROIidx
        autoROI(ridx).angle = polAngs(u(l));
    else
        autoROI(ridx).angle = nan;
    end
    
    
    % % %         numHalfCycles = 2*length(find(obj.TrialPatNum(obj.TrialSeqNum == expNum) == 360));
    % % %         maskPolResponse = maskResponse(polStart: polStop);
    % % %         CycleLims = floor(linspace(1,length(maskPolResponse),numHalfCycles+1));
    % % %         maskTrialResp = nan(numHalfCycles,max(diff(CycleLims)));
    % % %         for cycIdx = 1:numHalfCycles
    % % %             thisResp = maskPolResponse(CycleLims(cycIdx):CycleLims(cycIdx+1)-1);
    % % %             % Subtract off mean:
    % % %             maskTrialResp(cycIdx,1:CycleLims(cycIdx+1)-CycleLims(cycIdx)) = thisResp - nanmean(thisResp);
    % % %         end
    % % %         angleVector = linspace(polAngs(1),polAngs(0.5*length(polAngs)+1),size(maskTrialResp,2));
    % % %         lineProps.col = {col};
    % % %         mAx;
    % % %         mseb(angleVector,nanmean(maskTrialResp),nanstd(maskTrialResp),lineProps,1);
    % % %
    % % %         % add line to polar axes
    % % %         theta = deg2rad(polAngs(u(l)));
    % % %         rho = mean(polSelImg(k==u(l)));
    % % %         polarplot(pAx,[theta-pi theta],[rho rho],'-','color',col,'LineWidth',2)
    % % %
    
end
if ~isempty( [autoROI.mask] )
    if nargout
        varargout{1} = autoROI;
    else
        
        % assign to object
        % objROIidx = length(obj.ROI); % appends existing ROIs
        objROIidx = 0; % Overwrites existing ROIs
        
        if objROIidx == 0
            obj.ROI = struct('angle',[],'mask',[],'color',[],'position',[],'response',[]);
        end
        for roiIdx = 1:length(autoROI)
            objROIidx = objROIidx + 1;
            obj.ROI(objROIidx) = autoROI(roiIdx);
            %disp(['ROI assigned to object as ROI index ' num2str(objROIidx) ])
        end
        
    end
    
end

obj.fullAngMaskImg = fullMasks;
obj.exampleAngMaskImg = examplemasks;

%{
if manualROIs
    mFig = figure('Color','w','SizeChangedFcn',@resetAxes);
    
    % Frame with masks found for common pol angle/color
    maskAx = subplot(2,3,2);
    
    % pol selectivity image
    polSelAx = subplot(2,3,6);
    imagesc(polSelImg.*(f>10),'Parent',polSelAx) % threshold at 10 is just to get rid of edges
    colormap(polSelAx,flipud(gray))
    colormap(polSelAx,flipud(hot))
    
    title(['selectivity (Avg bins)[thresh>' num2str(polSelThreshold) ']'])
    axis(polSelAx,'image','off');
    polSelAx.Color = [0 0 0];
    colorbar
    polSelAx.CLim(1) = polSelThreshold;
    
    % pol mag image
    polMagAx = subplot(2,3,5);
    imagesc(fftMagImg,'Parent',polMagAx)
    colormap(polMagAx,flipud(gray))
    colormap(polMagAx,flipud(hot))
    
    title(['magnitude (FFT)[take top' num2str(polRespMagThreshold) '%]'])
    axis(polMagAx,'image','off');
    polMagAx.Color = [0 0 0];
    colorbar
    
    % median image
    oldAx = subplot(2,3,4);
    imagesc(f,'Parent',oldAx)
    title('median activity')
    colormap(oldAx,flipud(gray))
    colormap(oldAx,flipud(hot))
    axis(oldAx,'image','off');
    colorbar
    
    % Ang pref for each color mask, polar plot
    pAx = subplot(2,3,1,polaraxes);
    hold on
    set(pAx,'ThetaDir','clockwise','ThetaZeroLocation','right')
    pAx.RAxis.Visible = 'off';
    
    % Time-series for each color mask
    mAx = subplot(2,3,3);
    axis(mAx,'square');
    hold on
    mAx.XLim = [polAngs(1) polAngs(0.5*length(polAngs)+1)];
    mAx.XLabel.String = 'angle ({\circ})';
    mAx.YLabel.String = 'intensity (mean subtracted, a.u.)';
    set(mAx,'xtick',[polAngs(1:0.5*length(polAngs))] ,'xticklabel',[])
    
    u = unique(k(~isnan(k)));
    r = zeros( size(f,1),size(f,2),length(u) );
    % Invert f for displaying
    fInv =  (1-(f./max(f(:))));
    % Create ~mask, setting pixels that will become colored, polSel
    % regions to zero (since we will add values for color)
    fInvMask = fInv.*(isnan(k));
    fullMasks = cat(3,fInvMask,fInvMask,fInvMask);
    examplemasks = cat(3,fInvMask,fInvMask,fInvMask);
    %     masksbw = zeros( size(f,1),size(f,2));
    %     maskscale = masksbw;
    f=double(f);
    
    autoROI = struct('angle',[],'mask',[],'color',[],'position',[],'response',[]);
    ridx = 0;
    for l = 1:length(u)
        
        ridx = ridx+1;
        
        col = polColMap(u(l),:);
        
        rk = k==u(l);
        maskcol = squeeze(cat(3,col(1)*rk,col(2)*rk,col(3)*rk));
        fullMasks = fullMasks + maskcol;
        
        % Get the pixels included at this angle
        rk = k==u(l);
        rp = regionprops(rk,'PixelIdxList','Area');
        [~,rpSrtIdx] = sort([rp.Area],'descend');
        % Sort by area size and for the top X areas, get the idx of
        % all pixels contained within:
        rp = rp(rpSrtIdx);
        % Group the top numROIs areas by size:
        %                 %{
        if length(rp) < numROIs
            pixForColor = vertcat(rp(1:end).PixelIdxList);
        else
            pixForColor = vertcat(rp(1:numROIs).PixelIdxList);
        end
        %                 %}
        % -- changed to take the top area, only:
        %                 pixForColor = vertcat(rp(1).PixelIdxList);
        
        % Get the idx of all pixels currently set to zero in
        % 'masks', awaiting color (because not all will be colored
        % now)
        pixAtZero = find(k==u(l));
        % Find the pixels which won't be colored and should be
        % reset to original average image (inverted)
        pixToReset = setdiff(pixAtZero,pixForColor);
        
        % Update polAng array to discard all areas but the top X
        k(pixForColor) = u(l);
        k(pixToReset) = nan;
        
        % Update masks array to discard all areas but the top X
        tempMask = zeros( size(f,1),size(f,2));
        tempMask(pixToReset) = 1;
        toAdd = tempMask.*fInv;
        toAdd3 = cat(3,toAdd,toAdd,toAdd);
        examplemasks = examplemasks + toAdd3;
        
        rk = k==u(l);
        maskcol = squeeze(cat(3,col(1)*rk,col(2)*rk,col(3)*rk));
        examplemasks = examplemasks + maskcol;
        
        maskResponse = scanROI(obj,rk);
        
        autoROI(ridx).mask = rk;
        autoROI(ridx).color = col;
        autoROI(ridx).position = [1 1 1 1];
        autoROI(ridx).response = maskResponse;
        if u(l) ~= controlROIidx
            autoROI(ridx).angle = polAngs(u(l));
        else
            autoROI(ridx).angle = nan;
        end
        
        numHalfCycles = 2*length(find(obj.TrialPatNum(obj.TrialSeqNum == expNum) == 360));
        maskPolResponse = maskResponse(polStart: polStop);
        CycleLims = floor(linspace(1,length(maskPolResponse),numHalfCycles+1));
        maskTrialResp = nan(numHalfCycles,max(diff(CycleLims)));
        for cycIdx = 1:numHalfCycles
            thisResp = maskPolResponse(CycleLims(cycIdx):CycleLims(cycIdx+1)-1);
            % Subtract off mean:
            maskTrialResp(cycIdx,1:CycleLims(cycIdx+1)-CycleLims(cycIdx)) = thisResp - nanmean(thisResp);
        end
        angleVector = linspace(polAngs(1),polAngs(0.5*length(polAngs)+1),size(maskTrialResp,2));
        lineProps.col = {col};
        mAx;
        mseb(angleVector,nanmean(maskTrialResp),nanstd(maskTrialResp),lineProps,1);
        
        % add line to polar axes
        theta = deg2rad(polAngs(u(l)));
        rho = mean(polSelImg(k==u(l)));
        polarplot(pAx,[theta-pi theta],[rho rho],'-','color',col,'LineWidth',2)
        
        
    end
    image(examplemasks,'Parent',maskAx)
    axis(maskAx,'image','off');
    
    text(mAx,[polAngs(1:0.5*length(polAngs))+ 0.5*(mode(diff(polAngs)))].',(mAx.YLim(1)+20)*ones(0.5*length(polAngs),1), ...
        num2str([polAngs(1:0.5*length(polAngs))].','%d'),'horizontalalign','center');
    for pIdx = 2:0.5*length(polAngs)
        pL = line([polAngs(pIdx) polAngs(pIdx)],[mAx.YLim(1) mAx.YLim(2)],'LineStyle',':','Color',[0.8 0.8 0.8], 'LineWidth',2 );
        uistack(pL,'bottom')
    end
    
    
    
    
    if ~isempty( [autoROI.mask] )
        if nargout
            varargout{1} = autoROI;
        else
            
            % assign to object
            % objROIidx = length(obj.ROI); % appends existing ROIs
            objROIidx = 0; % Overwrites existing ROIs
            
            if objROIidx == 0
                obj.ROI = struct('angle',[],'mask',[],'color',[],'position',[],'response',[]);
            end
            for roiIdx = 1:length(autoROI)
                objROIidx = objROIidx + 1;
                obj.ROI(objROIidx) = autoROI(roiIdx);
                %disp(['ROI assigned to object as ROI index ' num2str(objROIidx) ])
            end
            
        end
        
    end
    
    obj.fullAngMaskImg = fullMasks;
    obj.exampleAngMaskImg = examplemasks;
 end
%}

%% Auto ROI stuff
%{
autoROI = struct('mask',[],'color',[],'position',[],'response',[]);

switch plotFlag
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
        for kS = 1:length(unique(k(~isnan(k))))
            plot(dotAx,kS,mean(polSelImg(k==kS)),'o','color',ROIcolor(kS),'MarkerFaceColor',ROIcolor(kS))
            hold on
        end
        dotAx.YLim(1) = 0;
        
        subplot(2,2,[2])
        mAx = gca;
        axis(mAx,'square');
        hold on
        
        % Add listener to keep images at max size
        % addlistener(mAx, 'Resize', @(obj,event)resetAxes(maskAx,avgAx));
        
        ridx = 0;
        u = unique(k(~isnan(k)));
        r = zeros( size(f,1),size(f,2),length(u) );
        masks = zeros( size(f,1),size(f,2));
        masksbw = zeros( size(f,1),size(f,2));
        maskscale = masksbw;
        f=double(f);
        for l = 1:length(u)
            rk_scale =  mean(mean(f(k==u(l))))./max(f(:));
            if rk_scale*max(f(:)) > 0%mean(f(:)) +std(f(:))
                ridx = ridx+1;
                rk = k==u(l);
                col = ROIcolor(ridx);
                maskcol = squeeze(cat(3,col(1)*rk,col(2)*rk,col(3)*rk));
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
            u = unique(k);
            r = zeros( size(f,1),size(f,2),length(u) );
            masks = zeros( size(f,1),size(f,2));
            masksbw = zeros( size(f,1),size(f,2));
            maskscale = masksbw;
            f=double(f);
            for l = 1:length(u)
                rk_scale =  mean(mean(f(k==u(l))))./max(f(:));
                if rk_scale*max(f(:)) > mean(f(:)) +std(f(:))
                    ridx = ridx+1;
                    rk = k==u(l);
                    col = ROIcolor(ridx);
                    maskcol = cat(3,col(1)*rk,col(2)*rk,col(3)*rk);
                    masks = masks + maskcol;
                    masksbw = u(l).*rk + masksbw;
                    maskscale = rk_scale.*rk + maskscale;
                    r(:,:,ridx) = rk;
                    plot(mAx,scanROI(obj,rk),'Color',col,'LineWidth',2)
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

delete(mFig)
%}

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

angleFrameArray = logical(angleFrameArray);
end

function angleFrameArray = getAngleLogicArray2(obj,polAng)

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

% Don't add orthogonal angle frames

% Construct logical array for frameAngles == polAng
angleFrameArray = zeros(1,size(obj.Frames,3));
for n = 1:length(angStarts)
    if (angStarts(n) >= polStart) && (angStops(n) <= polStop)
        angleFrameArray(angStarts(n):angStops(n)) = 1;
    end
end

angleFrameArray = logical(angleFrameArray);
end


function connMask = smoothdiscreteTuning(discreteTuning)

discreteTuningCopy = discreteTuning;
connMask = zeros(size(discreteTuning));

minArea = 16;

uniqueAng = unique(discreteTuning(discreteTuning>0));
for n = 1:length(uniqueAng)
    I = discreteTuningCopy == uniqueAng(n);
    connMask = connMask +  uniqueAng(n).*bwareaopen(I,minArea);
end

smallAreas = discreteTuningCopy.*(discreteTuningCopy & ~connMask);
plusminus1Ang(:,1) = circshift(uniqueAng,1);
plusminus1Ang(:,2) = circshift(uniqueAng,-1);
for n = 1:length(uniqueAng)
    I = smallAreas == uniqueAng(n);
    rp = regionprops(I,'PixelIdxList');
    
    % Examine each small area individually
    for rIdx = 1:length(rp)
        % Get the surrounding perimeter of the area
        areaPix = zeros(size(discreteTuning));
        areaPix(rp(rIdx).PixelIdxList) = 1;
        perimPix = find(imdilate(areaPix, true(3)) - areaPix);
        perimPixVals = discreteTuningCopy(perimPix);
        
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
            [arow,acol] = ind2sub(size(discreteTuning),perimPix( perimPixVals==mode(perimPixVals(perimPixVals>0)) ));
            outerPix = find(grayconnected(discreteTuningCopy,arow(1),acol(1)));
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
%}