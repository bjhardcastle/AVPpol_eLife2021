function superSaveMSP(objarray,overwrite)
%superSaveMSPFrames For super objects, make MSP frames and save to disk as a mat file
% (Not currently used)
if nargin <2 || isempty(overwrite)
    overwrite = 0;
end
for oidx = 1:length(objarray)
    if any(objarray(oidx).Exps==2|objarray(oidx).Exps==4|objarray(oidx).Exps==8|objarray(oidx).Exps==10)
        if ~isempty(objarray(oidx).MIP)
            filename = [objarray(oidx).MIP.DateStr objarray(oidx).MIP.TimeStr '_MSP.mat'];
            folder = [objarray(oidx).MIP.Folder];
            if ~exist(fullfile(folder,filename)) || overwrite
                makeMISP(objarray(oidx))
                Frames = objarray(oidx).MIP.Frames;
                if ~isempty(Frames)
                    save( fullfile(folder,filename), 'Frames' ,'-v7.3');
                    objarray(oidx).MIP.UseMSP = 1;
                end
            end
        end
    end
end
