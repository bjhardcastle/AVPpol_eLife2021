function assignPolPixRGB(objarray)
% Quick function to make color-mapped discrete tuning map in an
% activityframe, for drawing ROIs
for oidx= 1:length(objarray) 
    discCol = objarray(oidx).MIP.polPixDiscrete; % get discrete tuning map (With layerMask
    avg = objarray(oidx).MIP.InactivityFrame; % add static grayscale image
    avg = avg./max(avg(:));
    discCol(isnan(discCol)) = -180*avg(isnan(discCol));
    discCol = (discCol + 179);
    im =ind2rgb(round(discCol),([flipud(gray(180));flipud(hsv(180))]) );
    
    objarray(oidx).MIP.ActivityFrame2 = im;
    
    if ~isempty(objarray(oidx).staticMIP)

    % also add an average image from the static channel
im = objarray(oidx).staticMIP.AverageFrame;
if isempty(im)
getFrames(objarray(oidx).staticMIP)
im = objarray(oidx).staticMIP.AverageFrame;
objarray(oidx).staticMIP.Frames = [];
end
im = 255.*(im./max(im(:)));
objarray(oidx).MIP.ActivityFrame3 = ind2rgb(round(im),([1 0 0].*ones(255,3).*linspace(0,1,255)') );
end
end

  discCol = objarray(oidx).MIP.polTuningImg; % get discrete tuning map (With layerMask
    im =ind2rgb(round(discCol),([flipud(hsv(180))]) );
    
    objarray(oidx).MIP.ActivityFrame4 = im;


end