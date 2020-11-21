function assignStackROIs(objarray)        
if isempty([objarray.TrialPatNum])
    runAcrossPlanes(objarray,'getParameters')
end

angles = [30 60 90];
angleStk = zeros(size(objarray(1).AverageFrame,1),size(objarray(1).AverageFrame,2),objarray(1).numZPlanes,length(angles));

for oidx = 1:length(objarray)
    objarray(oidx).Unattended = 1;
    getFrames(objarray(oidx))
    
    for aidx = 1:length(angles)
            % Get activity frame for each angle, a, equal to:
            % average frame during presentation of [a,a+180] minus
            % average frame during [a+90,a+270]
            aFrame = getPolAngleActivityFrame(objarray(oidx),angles(aidx));
        
            % Threshold activity frame to find pixels which respond to this
            % angle, and pixels which respond to angle+/-90
            maskStack1  = double((aFrame./max(aFrame(:)))>0.1);
            maskStack2  = double((aFrame./min(aFrame(:)))>0.1);
            
            aThresh = abs(aFrame) > (mean(abs(aFrame(:))) + 3*std(abs(aFrame(:))));
            aPos = aFrame.*(aFrame>0);
            aNeg = -aFrame.*(aFrame<0);
            mask1  = logical((aPos.*aThresh)>0);
            mask2  = logical((aNeg.*aThresh)>0);
            
            % Assign to object:
            objarray(oidx).ROI(aidx).mask = mask1;
            objarray(oidx).ROI(aidx).response = scanROI(objarray(oidx),mask1);

            objarray(oidx).ROI(aidx + length(angles)).mask = mask2;
            objarray(oidx).ROI(aidx+ length(angles)).response = scanROI(objarray(oidx),mask2);
            
            objarray(oidx).UseFixedResp = 1;
    end
    
    objarray(oidx).Frames = [];
    objarray(oidx).Unattended = 0;
end


end