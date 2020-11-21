function updateROIs(objarray)
%UPDATEROIS Update the fixed response stored for each ROI and create a mask
%if missing

for oidx = 1:length(objarray)
    
    if isempty(objarray(oidx).Frames)
        % Get object's frames
        getFrames(objarray(oidx)) % Careful of the setting objarray(oidx).UseBackSubFrames
        CLEAR_FRAMES_FLAG = 1;
    else
        CLEAR_FRAMES_FLAG = 0;
    end
    
    % Store the current setting for using fixed resp, as we will have to
    % turn it off temporarily
    FIXED_RESP_FLAG = objarray(oidx).UseFixedResp;
    
    objarray(oidx).UseFixedResp = 0;
    for ridx = 1:length(objarray(oidx).ROI)
        
        % if no mask exists yet, create one
        if ~isfield(objarray(oidx).ROI(ridx),'mask') ...
                || isempty(objarray(oidx).ROI(ridx).mask)
            maskIm = zeros(size(objarray(oidx).AverageFrame));
            mfig = figure;
            imshow(maskIm);
            roiMask = impoly(gca,objarray(oidx).ROI(ridx).position);
            objarray(oidx).ROI(ridx).mask = createMask(roiMask);
            close(mfig)
        end
        
        scanROI(objarray(oidx),objarray(oidx).ROI(ridx).mask);
    end
    
    % Restore fixed resp setting
    objarray(oidx).UseFixedResp = FIXED_RESP_FLAG;
    
    if CLEAR_FRAMES_FLAG
        % Clear frames from memory
        objarray(oidx).Frames = [];
    end
end