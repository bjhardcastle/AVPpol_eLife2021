function out = scanROI(mov, mask)
%SCANROI Find mean value within roi across frames
%   out = scanROI(obj, mask, scan_extents) finds the mean value of the
%   pixels within the passed ROI for each frame within the passed SlidebookObj. 
%   This function accepts the following arguments:
%
%       mask - A N-dimensional logical array the same size as the image
%           data with true values for pixels that should be included in the
%           ROI. For efficiency, all image data must have the dimensions.
%           If image data is encountered with a different dimension, then
%           that slice is given a value of NaN.
%
%       scan_extents - A [2x1] array that specifies how to scan through the
%           position data. The first element specifies the first frame to scan, 
%           and the second element specifies the last frame to scan to. 
%           If one or more entries is missing, the entire tiff will be
%           scanned.
%
%
%   This function has the following outputs:
%
%       out - A [3xN] vector where N is equal to the number of processed
%           slices. The first row contains the slice index from which the
%           ROI was calculated, and the second row contains the time
%           associated with the slice, and the third row contains the
%           calculated ROI value.
%
%        mean pixel intensity
%

% Check input arguments

%    if empty or only [1x1] then discard and use full extent

% Check ROI is same size as tiff frames
% Check if Frames have already been scanned into object



%   If no ROI specified then use entire image frame
if isempty(mask)
    disp('No ROI specified: finding mean of entire image')
    mask = true(size(mov,1),size(mov,2));
end

% If an ROI was specified it must be the same size as the tiff image frames
assert(isequal(size(mask),size(mov(:,:,1))) ,'ROI size does not match dimensions of tiff frames')
if ~islogical(mask)
    mask = logical(mask);
end

% for n = 1:size(mov,3)
%     % Grab the current frame
%     currentframe = mov(:,:,n);    
%     % Extract the mean pixel intensity within the ROI
%     out(n) = mean(currentframe(mask));
% end
% % Above code, vectorized:
out = sum(reshape(mov.*mask,length(mask(:)),size(mov,3)))./sum(mask(:));
  