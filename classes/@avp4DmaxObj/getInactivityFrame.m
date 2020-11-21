function getInactivityFrame(objarray)
for oidx = 1:length(objarray)
    if isempty(objarray(oidx).TrialSeqNum)
        getParameters(objarray(oidx))
    end
    clearFrames = 0;
    if isempty(objarray(oidx).Frames)
        getFrames(objarray(oidx));
        clearFrames = 1;
    end
    exps = unique(objarray(oidx).TrialSeqNum);
    frameIdx = ones(1,length(objarray(oidx).Frametimes));
    
    framesec = ceil(1/(objarray(oidx).IFI./objarray(oidx).AIrate));
    
    for eIdx = exps
        frameIdx(objarray(oidx).expStart(eIdx)-framesec : objarray(oidx).expStop(eIdx)+framesec ) = 0;
    end
    if length(objarray(oidx).Frametimes) > size(objarray(oidx).Frames,3) ...
            || length(frameIdx) > size(objarray(oidx).Frames,3)
        frameIdx = frameIdx(1:size(objarray(oidx).Frames,3));
    end
    objarray(oidx).InactivityFrame = mean(objarray(oidx).Frames(:,:,logical(frameIdx)),3);
    
    if clearFrames
        objarray(oidx).Frames = [];
    end
end
end