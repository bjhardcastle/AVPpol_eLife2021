function [tuningAng,resultantLength] = tuningROI(obj,ROImask,convertToPubAngs)
if nargin<3 || isempty(convertToPubAngs)
    convertToPubAngs = 0;
end
if isempty(obj.polAngResp)
    getPolMaps(obj)
end
if ~isempty(obj.polExp)
    % Find circular mean for each pixel:
    ROIresp = mean(obj.polAngResp(ROImask,:),1);
    if length(ROIresp) == 6 
        % This was a half-cycle response
            polAngsOrig = [180/length(ROIresp):180/length(ROIresp):180];
    else
    polAngsOrig = [360/length(ROIresp):360/length(ROIresp):360];
    end
    [actualROIVec,axtheta] = circ_rangle(circ_axial(deg2rad(repmat(polAngsOrig(1:length(polAngsOrig)),size(ROIresp,1),1)),2), ROIresp, deg2rad(mode(diff(polAngsOrig))), 2);
    theta = axtheta/2;
    Ang =0.5.*(wrapTo360((2.*(rad2deg(theta)))));
    resultantLength = actualROIVec;
    tuningAng = Ang; % Expressed in original angles
    if convertToPubAngs
        tuningAng = wrapTo180(-tuningAng-270);
    end
else
    resultantLength = nan;
    tuningAng = nan;
end