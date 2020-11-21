function F0 = findPolF0(obj,ROI,expNum)
assert(nargin==3,'Must input expNum and ROI') 
assert(~isempty(expNum) && isnumeric(expNum), 'Must specify expNum as a numerical value')
assert(~isempty(ROI), 'Must provide an ROI strucutre or mask')
% Deal with ROI input which can be either a mask or a regular ROI struct
if isstruct(ROI)
    if ~isfield(ROI,'mask') || isempty(ROI.mask)
       disp('ROI mask is missing')
       return
    end
elseif isnumeric(ROI) && length(unique(ROI))==2 % should be 0s and 1s only, but may not be logical
    mask = ROI;
    clear ROI
    ROI = struct;
    ROI.mask = mask;
   
    if isempty(obj.MIP.Frames)
        getFrames(obj.MIP)
    end
    ROI.response = scanROI(obj.MIP,mask);
end

    if expNum ~= 4 || ( expNum == 4 && isempty(obj.MIP.polActivityImg) )
        if isempty(obj.MIP.Frames)
            getFrames(obj.MIP)
        end
        meanPol = mean(obj.MIP.Frames(:,:,obj.MIP.expStart(expNum): obj.MIP.expStop(expNum)),3);
    else
        meanPol = obj.MIP.polActivityImg;        
    end
    
    % Find mean of lowest 10% of pix:
    
    maskPix = sort(meanPol(ROI.mask));
    
    F0 = mean(maskPix(1:floor(length(maskPix)/10)));
end



