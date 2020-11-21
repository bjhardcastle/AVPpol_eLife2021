function varargout = getCurvilinearDist(obj,I,J,refreshMask)
% I,J are pixel indices within the layerMask: can be a subset, in which
% case we just return the corresponding distances for those pixels. If
% refreshing, we must calculate distances for all pixels then look at
% subset.
if nargin < 2 || isempty(I) || isempty(J)
    useWholeMask = 1;
else
    useWholeMask = 0;
end



% special case for E-PG recordings in the bridge:
% map x,y pixel coordinates onto a curve (from left to right across
% the bridge) and use distance along the curve as position x (and
% ignore y)
assert(~isempty(obj.layerMask),['PB fitting requires layerMask, objID: ' obj.DateStr obj.TimeStr ])

if nargin < 4 || isempty(refreshMask)
    if ~isprop(obj,'PBcurvilinearDistMask') || isempty(obj.PBcurvilinearDistMask) ...
            || sum(~isnan(obj.PBcurvilinearDistMask(:))) ~= sum(obj.layerMask.mask(:)) % if mask changes shape we must recalculate distances to curve
        refreshMask = 1;
    else
        refreshMask = 0;
    end
end

if refreshMask
    
    
    objMask = obj.layerMask.mask;
    maskOutline = obj.layerMask.position; %
    
    % get the indices of all pixels within the layerMask
    [y,x] = ind2sub( size(objMask), find(objMask) ); % [ x,y] or [vertical,lateral] in brain
    
    
    if strcmp(obj.Area,'Me')
        
        % For recordings across a layer of the medulla the points fall
        % approximately in a line when viewed dorsally: fit a line to
        % all points in the MIP image then rotate pixel positions (old
        % superXY code)
        [j,i] = ind2sub(size(objMask),find(objMask)); % [i,j] or [vertical,lateral] in brain
        p = polyfit(j,i,1);
        rotateAngle = -atand(p(1)) ;
        
        % Rotate image and extract data
        jj = [x*cosd(rotateAngle) + y*sind(rotateAngle)];
        ii = [-x*sind(rotateAngle) + y*cosd(rotateAngle)];
        
        x = jj;
        y = ii + abs(min(ii));
        
    end
    
    
    % add some noise to their positions (since they're ordered in a
    % grid)
    xn = x+randn(size(x));
    yn = y+randn(size(y));
    
    
    % First we fit a line across the PB. Settled on piecewise
    % linear fit as it seems robust across range of recordings
    
    
    % Get some x/horizontal points from left to right which will form
    % the breakpoints of the piecewise fit. Here the ends of the PB,
    % where they turn to point ventrally are the challening part - the
    % fitting wants to continue horizontally. Taking slightly less than
    % the most extreme left/right positions as the limits works well
    % enough to end the curve roughly in the last glomeruli on either
    % end
    [right,rightYidx] = max(y);
    [left,leftYidx] = min(y);
    width = right - left;
    YI = linspace(left+0.02*width,right-0.02*width,20)';
    [ XI ] = lsq_lut_piecewise(yn, xn,YI );
    
% find y-coord of nearest mask outline points to left and right limits
[~,leftPt] = min(abs(maskOutline(:,2)-left));
[~,rightPt] = min(abs(maskOutline(:,2)-right));

%%
% find the first inflection point moving in either direction along x-coords
% from these limits. 
% first rearrange outline points so they start at the left/right limits,
% then find the sign of the gradient at each point, and changes in sign
dir1FromLeftPt = circshift(maskOutline(:,1),-leftPt+1);
dir2FromLeftPt = flipud(circshift(maskOutline(:,1),-leftPt+1));
dir1FromRightPt = circshift(maskOutline(:,1),-rightPt+1);
dir2FromRightPt = flipud(circshift(maskOutline(:,1),-rightPt+1));

dir1FromLeftPt2ndDeriv = diff( sign( diff( dir1FromLeftPt )));
dir2FromLeftPt2ndDeriv = diff( sign( diff( dir2FromLeftPt )));
dir1FromRightPt2ndDeriv = diff( sign( diff( dir1FromRightPt )));
dir2FromRightPt2ndDeriv = diff( sign( diff( dir2FromRightPt )));

% check that the limit we've found isn't an inflection point itself first -
% diff operation won't give results for the points in that neighborhood
if any( diff(sign(diff(dir1FromLeftPt([end-1,end,1,2])))) )
    dir1FromLeftInflectIdx = 1;
    dir2FromLeftInflectIdx = 1;
else
    dir1FromLeftInflectIdx = find(dir1FromLeftPt2ndDeriv,1,'first') +1;
    dir2FromLeftInflectIdx = find(dir2FromLeftPt2ndDeriv,1,'first') +1;
end
if any( diff(sign(diff(dir1FromRightPt([end-1,end,1,2])))) )
    dir1FromRightInflectIdx = 1;
    dir2FromRightInflectIdx = 1;
else
    dir1FromRightInflectIdx = find(dir1FromRightPt2ndDeriv,1,'first') +1;
    dir2FromRightInflectIdx = find(dir2FromRightPt2ndDeriv,1,'first') +1;
end

dir1FromLeftInflect = dir1FromLeftPt(dir1FromLeftInflectIdx);
dir2FromLeftInflect = dir2FromLeftPt(dir2FromLeftInflectIdx);
dir1FromRightInflect = dir1FromRightPt(dir1FromRightInflectIdx);
dir2FromRightInflect = dir2FromRightPt(dir2FromRightInflectIdx);

% whichever point is more ventral (should have a positive inflection but
% might not be reliable) is the one we want: rearrange points around that
% bottom left corner
if dir1FromLeftInflect < dir2FromLeftInflect 
    newYPtsFromBottomLeftDir1 = circshift(dir1FromLeftPt, -dir1FromLeftInflectIdx+1);
    % repeat what we did for the Y points above for the X points
    newXPtsFromBottomLeftDir1 =  circshift( circshift(maskOutline(:,2),-leftPt+1) ,-dir1FromLeftInflectIdx+1);
    
elseif dir2FromLeftInflect < dir1FromLeftInflect
    newYPtsFromBottomLeftDir1 = circshift(dir2FromLeftPt, -dir2FromLeftInflectIdx+1);
    newXPtsFromBottomLeftDir1 =  circshift( flipud(circshift(maskOutline(:,2),-leftPt+1)), -dir2FromLeftInflectIdx+1);
end
newYPtsFromBottomLeftDir2 = circshift(flipud(newYPtsFromBottomLeftDir1),1);
newXPtsFromBottomLeftDir2 = circshift(flipud(newXPtsFromBottomLeftDir1),1);

% now we need the index of the bottom right corner in each of these new
% vectors, so we can take all points between left corner and right corner,
% travelling in both directions
if dir1FromRightInflect < dir2FromRightInflect 
    rightBottomCornerYval = dir1FromRightInflect;
elseif dir2FromRightInflect < dir1FromRightInflect
    rightBottomCornerYval = dir2FromRightInflect;
end

% get the two curves, which both start from the bottom left corner/tip
rightIdxDir1 = find(newYPtsFromBottomLeftDir1 == rightBottomCornerYval,1,'last');
rightIdxDir2 = find(newYPtsFromBottomLeftDir2 == rightBottomCornerYval,1,'first');

semiOutline1 = [newYPtsFromBottomLeftDir1(1:rightIdxDir1) newXPtsFromBottomLeftDir1(1:rightIdxDir1)];
semiOutline2 = [newYPtsFromBottomLeftDir2(1:rightIdxDir2) newXPtsFromBottomLeftDir2(1:rightIdxDir2)];

%{
% To check curves
figure,
plot(semiOutline1(:,1),semiOutline1(:,2))
hold on
plot(semiOutline2(:,1),semiOutline2(:,2))
%}
% now resample the curves to a common sampling frequency
samplingPts = linspace(0,1,21); % must be the same
[intPts1] = interparc(samplingPts,semiOutline1(:,1),semiOutline1(:,2),'linear');
[intPts2] = interparc(samplingPts,semiOutline2(:,1),semiOutline2(:,2),'linear');

% construct a new curve from the midpoint between each pair of points,
% essentially running through the center of the two previos semi outlines
midCurveX = mean([intPts1(:,1),intPts2(:,1)],2);
midCurveY = mean([intPts1(:,2),intPts2(:,2)],2);
%{
figure
hold on
plot(intPts1(:,1),intPts1(:,2))
plot(intPts2(:,1),intPts2(:,2))
plot(midCurveX, midCurveY )
%}
XI = midCurveX;
YI = midCurveY;
%%
    % To check the fit:
    %{
        figure,
        plot(x,y,'b.')
        hold on
        plot(XI,YI,'r')
    %}
    
    % now find the closest point of the piecewise fit to each pixel in
    % the layer mask, and find the distance of that point along the fit
    % (left being 0, right being 1) - this is our new X position
    
    mapxy = [y,x];
    curvexy = [YI,XI];
    
    [closestCurvePt,~,dist] = distance2curve(curvexy,mapxy,'linear');
       % To check against regular projection of pixel position onto horizontal axis
    %{
        
        figure
        
        subplot(1,2,1)
        cm = flipud(hsv(100));
        cm = [cm; cm];
        Yup = round(((y-left)./(right-left)).*199);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
        scatter(y,x,100,C,'filled')
        axis image
        title('horizontal position')

        subplot(1,2,2)
        cm = flipud(hsv(100));
        cm = [cm; cm];
        Yup = round(dist.*199);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
        scatter(y,x,100,C,'filled')
        axis image
        title('position along curve')
    %}
 
    newMask = nan(size(obj.layerMask.mask));
    newMask(find(obj.layerMask.mask)) = dist;
    

    % We want to skew the points in the distance mask so that at the
    % midline the value is 0.5 
        
    % Get the coordinates of curve at the brain midline
    midPointX = round(mean([right,left]));
    [~,minIdx]= min(abs(curvexy(:,1)-midPointX));
    % Get the distance along the curve at the midline 
    midVal = newMask(floor(curvexy(minIdx,2)),floor(curvexy(minIdx,1)));
    midVal = nanmean(newMask(midPointX,:)) ;
    
    if ~isnan(midVal)
        gamma = log(midVal)/log(0.5); % gamma < 1 will bias towards left side (Zero)
    else % failsafe 
        gamma = 1;
    end
    obj.PBcurvilinearDistMask = newMask.^gamma;
    
    % To check against uncorrected mask 
    %{
        
        figure
        
        subplot(1,2,1)
        cm = flipud(hsv(101));
        cm = [cm; cm];
        Yup = round(dist.*199);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
        scatter(y,x,100,C,'filled')
        axis image
        title('uncorrected distance')

        subplot(1,2,2)
        cm = flipud(hsv(100));
        cm = [cm; cm];
        Yup = round((dist.^gamma).*199);
        Yup(isnan(Yup)) = 0;
        C = cm(Yup+1,:);
        scatter(y,x,100,C,'filled')
        axis image
        title('gamma corrected')
    %}
    
end

% Now we're sure we have the curvilinear distances, we can return
% pixel values if they were requested
if nargout >0
    wholeMask = obj.PBcurvilinearDistMask ;
    if useWholeMask
        varargout{1} = wholeMask;
    else
        varargout{1} = wholeMask(sub2ind(size(obj.layerMask.mask),I,J));
    end
end
if nargout >1
    varargout{2} =curvexy; 
end
if nargout > 2 
    closestCurvePtMapCell{1} = nan(size(obj.layerMask.mask));
    closestCurvePtMapCell{2} = nan(size(obj.layerMask.mask));

    closestCurvePtMapCell{2}(find(obj.layerMask.mask)) = closestCurvePt(:,2);
    closestCurvePtMapCell{1}(find(obj.layerMask.mask)) = closestCurvePt(:,1);

    varargout{3} = closestCurvePtMapCell;
end

end
