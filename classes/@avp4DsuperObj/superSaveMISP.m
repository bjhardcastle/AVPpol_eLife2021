function superSaveMISP(objarray,overwrite)
%superSaveMISPFrames For super objects, make MISP frames and save to disk as a mat file
% (Currently used)
if nargin <2 || isempty(overwrite)
    overwrite = 0;
end
for oidx = 1:length(objarray)
        if ~isempty(objarray(oidx).MIP)
            filename = [objarray(oidx).MIP.DateStr objarray(oidx).MIP.TimeStr '_MSP.mat'];
            folder = [objarray(oidx).MIP.Folder];
            if ~exist(fullfile(folder,filename)) || overwrite
                superUseMSP(objarray(oidx),0)
                makeMISP(objarray(oidx))
                Frames = objarray(oidx).MIP.Frames;
                if ~isempty(Frames)
                    save( fullfile(folder,filename), 'Frames' ,'-v7.3');
                    objarray(oidx).MIP.UseMSP = 1;
                end
                objarray(oidx).MIP.Frames = [];
            end
        end
end
