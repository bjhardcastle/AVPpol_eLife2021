function [xyRectangle,zLim] = getCropBox(objarray)
% Return a rectangle which can be used directly in 'imcrop' to crop image
% to the outer boundary of all ROI masks across planes. If a second output
% argument is requested zLimits, indicating planes with non-empty ROIs,
% are also returned.

if isempty([objarray.ROI])
    assignStackROIs(objarray)
end

nonemptyPlane = [1:length(objarray)];
sumROI = zeros(size(objarray(1).AverageFrame));
for oidx = 1:length(objarray)
    for ridx = 1:size(objarray(oidx).ROI,2)
        sumROI = [sumROI + objarray(oidx).ROI(ridx).mask];
    end
    if sum(sum([objarray(oidx).ROI.mask])) == 0
        nonemptyPlane(oidx) = [];
    end
end

activityBinaryImg = sumROI > 0;

% Dectect areas of activity
bw = regionprops(activityBinaryImg);
% Sort by size
[~,sortIdx] = sort([bw.Area],'descend');
% Find the bounding box of the largest detected region
xyRectangle = bw(sortIdx(1)).BoundingBox;

if nargout == 2
    zLim = [min(nonemptyPlane) max(nonemptyPlane)];
end

end