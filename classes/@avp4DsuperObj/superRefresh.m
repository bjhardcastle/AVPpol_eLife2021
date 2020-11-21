function superRefresh(objarray)
% % get all current values for object
superSaveMISP(objarray,1) % Use inactivity frame, select max per layer
superUseMSP(objarray,1)
for oidx = 1:length(objarray)
   
    for lidx = 1:length(objarray(oidx).Layers)
        objarray(oidx).Layers(lidx).Unattended = 1;
        getFrames(objarray(oidx).Layers(lidx))
        getDaqData(objarray(oidx).Layers(lidx))
        
        getFrametimes(objarray(oidx).Layers(lidx))
        getTrialtimes(objarray(oidx).Layers(lidx))
        getParameters(objarray(oidx).Layers(lidx))
        
%         % %         getInactivityFrame(objarray(oidx).Layers(lidx)) % this is done in makeMISP
        
        loadLayerMasks(objarray(oidx).Layers(lidx))
        getPolMaps(objarray(oidx).Layers(lidx))
        objarray(oidx).Layers(lidx).polROI = [];
        objarray(oidx).Layers(lidx).UseFixedResp = 0;
        % %getPolROIs(objarray(oidx).Layers(lidx))
        % %savePolROIs(objarray(oidx).Layers(lidx))
        % objarray(oidx).Layers(lidx).polROI = [];
        
        objarray(oidx).Layers(lidx).ROI = [];
        loadROIs(objarray(oidx).Layers(lidx))
        updateROIs(objarray(oidx).Layers(lidx))
        saveROIs(objarray(oidx).Layers(lidx))
        objarray(oidx).Layers(lidx).UseFixedResp = 1;
        objarray(oidx).Layers(lidx).ROI = [];
        objarray(oidx).Layers(lidx).layerMask = [];
        objarray(oidx).Layers(lidx).Frames = [];
        objarray(oidx).Layers(lidx).Daq = [];
    end
    
    if ~isempty(objarray(oidx).MIP)
        objarray(oidx).MIP.Unattended = 1;
        getFrames(objarray(oidx).MIP)
        getDaqData(objarray(oidx).MIP)
        
        getFrametimes(objarray(oidx).MIP)
        getTrialtimes(objarray(oidx).MIP)
        getParameters(objarray(oidx).MIP)
                
        if  (objarray(oidx).MIP.numZPlanes > 1) && exist([objarray(oidx).MIP.Folder,objarray(oidx).MIP.DateStr,objarray(oidx).MIP.TimeStr,'_MSP.mat'],'file')
%             % ~isempty(objarray(oidx).MIP.polExp)
            objarray(oidx).MIP.Frames = [];
            objarray(oidx).MIP.UseMSP = 1;
            getFrames(objarray(oidx).MIP)
        else
                        objarray(oidx).MIP.Frames = [];
            objarray(oidx).MIP.UseMSP = 0;
                    getFrames(objarray(oidx).MIP)
        end
        getInactivityFrame(objarray(oidx).MIP) % this is also done in makeMISP

% %         superUseMSP(objarray(oidx),1)
        loadLayerMasks(objarray(oidx).MIP)
        objarray(oidx).MIP.polROI = [];
        
        objarray(oidx).MIP.UseFixedResp = 0;
        getPolMaps(objarray(oidx).MIP)
        getPolROIs(objarray(oidx).MIP)
        savePolROIs(objarray(oidx).MIP)
        objarray(oidx).MIP.polROI = [];
        
        objarray(oidx).MIP.ROI = [];
        loadROIs(objarray(oidx).MIP)
        updateROIs(objarray(oidx).MIP)
        saveROIs(objarray(oidx).MIP)
        
        objarray(oidx).MIP.UseFixedResp = 1;
        
        objarray(oidx).MIP.ROI = [];
        objarray(oidx).MIP.layerMask = [];
        objarray(oidx).MIP.Frames = [];
        objarray(oidx).MIP.Daq = [];
    end
end
superUseMSP(objarray,1)
superPolThreshold(objarray,-1)

end