function assignExpActivityFrames(objarray)
for oidx = 1:length(objarray)
    if isempty(objarray(oidx).TrialPatNum)
        getParameters(objarray(oidx))
    end
    
    if isempty(objarray(oidx).Frames)
        getFrames(objarray(oidx))
    end
    
    expnum = unique([objarray(oidx).TrialSeqNum]);
    
    if length(expnum)>25
        disp(['num Exps is ' num2str(length(expnum)) '. Only room for 25 activity frames.'])
    end
    
    frameName = 'ActivityFrame';
            
    for pidx = 1:length(expnum)
        if pidx > 25
            break
        end
        fields = [];
        fields.TrialSeqNum = expnum(pidx);
        objarray(oidx).([frameName num2str(expnum(pidx))]) = getActivityFrame(objarray(oidx),fields);

    end
    
    objarray(oidx).Frames = [];

end