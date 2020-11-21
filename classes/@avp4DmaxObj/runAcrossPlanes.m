function runAcrossPlanes(objarray,funcstring)

eval([funcstring '(objarray(1))']);
for oidx = 2:length(objarray)
    
    objarray(oidx).Daq = objarray(oidx-1).Daq;
    objarray(oidx).Frames = objarray(oidx-1).Frames;
    
    objarray(oidx-1).Daq = [];
    objarray(oidx-1).Frames = [];
    
    eval([funcstring '(objarray(oidx))']);

end
objarray(1).Daq = objarray(oidx).Daq;
objarray(1).Frames = objarray(oidx).Frames;
