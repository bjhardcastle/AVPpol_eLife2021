function getFrametimes(obj)
%GETFRAMETIMES Detect frame markers in DAQ data and store in object
% getFrametimes(obj)
%
% Reads the data stored in obj.Daq(:,1) and looks for peaks that correspond
% to frame markers output by 3i 2P system. The sample index of each frame 
% is stored in the object. Use the functions prevFrame and nextFrame to
% find the closest frame marker to a particular sample index.
%
% This functions writes to object:
%   obj.Frametimes  -   [1 x numFrames] array of sample indices. Each entry
%                       corresponds to a Tiff frame detected in obj.Daq(:,1)
%
%   obj.IFI         -   double, equal to the median number of samples
%                       between each frame detected. Generally used for
%                       finding a Tiff's framerate.
%
% See also getDaqData, nextFrame, prevFrame.
if isempty(obj.DaqFile)
    disp('No DAQ file exists - try running ''getDaqFile(obj)''')
    return
end

if isempty(obj.Daq)
    getDaqData(obj);    
end

% Some abbreviations for readability:
wstim = obj.Daq;
airate = obj.AIrate;

% For each tiff frame, DAQ AI0 records a pulse: 
wframes=find(diff(round(wstim(:,1)/5))==1)+1;
% wframes(1) first grabbed frame
% wframes(end) last grabbed frame

% Get rid of any duplicate frames resulting from double blips in signal
wframes(find(diff(wframes)<0.1*mean(diff(wframes))) + 1) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOD FOR AVP 4D TIFFs, WHERE A MARKER WAS RECORDED FOR EVERY PLANE
    if ~isempty(obj.numZPlanes) && ~isempty(obj.ZPlane)
    % multiple avp tiffs (xyt): 
    % Keep frame markers associated with the plane of the current object
    wframes = wframes(obj.ZPlane:obj.numZPlanes:end);    
    elseif ~isempty(obj.numZPlanes) && contains(obj.File,'_max_')
        % max int. projection object
        wframes = wframes(1:obj.numZPlanes:end);
    elseif ~isempty(obj.numZPlanes) && contains(obj.File,'_Zmax_')
        % max int. projection object
        wframes = wframes(1:obj.numZPlanes:end);
    elseif obj.numZPlanes == 1 
        disp('Object has 1 z-plane and does not appear to be a 4D-stack');      
    else
        disp('Error occurred in getFrametimes')
        return
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get inter-frame interval (in samples)
ifi = (median(diff(wframes)));

% Push to object
obj.IFI = ifi;
obj.Frametimes = wframes';
    
    
