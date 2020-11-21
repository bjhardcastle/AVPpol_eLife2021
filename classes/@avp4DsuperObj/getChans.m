function getChans(obj)
% Channel number labels out of Slidebook are inconsistent between different
% recording modes. Depending on the number of layers recorded we can deduce
% which type of recording the current object was (assumes that a single
% plane image is never produced by a 4D Capture with a single layer, but
% instead used the regular single-plane Stream to Disk method)
% Assign to object as numbers (0,1) and strings ('C0','C1')
GCAMP = [];
tdTOMATO = [];

if obj.numZPlanes == 1 && obj.Channels == 2
    % Single layer recording, channel numbers are different ("Stream"
    % output from Slidebook)
    GCAMP = 0;
    tdTOMATO = 1;
    
elseif obj.numZPlanes >2 && obj.Channels == 2
    % Regular 4D avp recording ("Capture" output from Slidebook)
    GCAMP = 1;
    tdTOMATO = 0;
    
elseif obj.Channels == 1
    % No static channel was recorded
    GCAMP = 0;
    
end

obj.activeChan = GCAMP;
obj.staticChan = tdTOMATO;
obj.activeChanStr = ['C' num2str(GCAMP)];
if ~isempty(obj.staticChan)
    obj.staticChanStr = ['C' num2str(tdTOMATO)];
else
    obj.staticChanStr = [];
end
end