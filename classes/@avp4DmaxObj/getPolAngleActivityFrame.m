function varargout = getPolAngleActivityFrame(obj,angle)
%GETACTIVITYFRAME Gets the average activity frame during trials
% [image] = getActivityFrame(obj,fields,detrend)
%
% This function finds the mean Tiff frame recorded during trials, subtracts
% from it the mean of all extra-trial frames, to enhance activity
% correlated with stimulus presentation, and saves the image to the object.
% If fields is provided as an input, specific trials will be searched for,
% and only their frames will be included ( see findTrials ). 
%
% Writes to object (if 'fields' is not specified): 
%   obj.ActivityFrame   - can be toggled on/off in 'play' to aid drawing ROIs
%
% Additional images can be also be assigned manually:
%   obj.ActivityFrame1
%   obj.ActivityFrame2
%   etc.
% and can be displayed in the 'play' GUI.
%
% This function accepts the following inputs:
% 
%   fields      - (optional) structure with trial parameters and value. 
%                 eg. fields.TrialPatNum = 2 
%                 If provided, the activity frame will only include the 
%                 trials specified in fields.
%
%   detrend     - (optional)  Set to 1 to improve signal to noise by
%                 subtracting constant linear trend or the mean from all 
%                 frames. Similar to background subtraction. Takes time. 
%                 Default is 0.
%   
% This function returns the following outputs:
%
%   image       - (optional) 2D image data of the activity frame for the
%                 trials requested in 'fields'
%
% See also findTrials.

if isempty(obj.TrialStartSample) 
    getTrialtimes(obj)
end
if isempty(obj.DaqFile)
    return
end
if isempty(obj.TrialPatNum)
	getParameters(obj)
end
if isempty(obj.Frames)
    if ~obj.Unattended 
        undoflag = 1;
    else
        undoflag = 0;
    end
    obj.Unattended = 1;
    getFrames(obj);
    if undoflag
        obj.Unattended = 0;
    end
end

assert( nargin > 1 && ~isempty(angle),'Must specify an angle')
    fields=[];
    fields.TrialSeqNum = 4;
    fields.TrialPatNum = angle;
    if fields.TrialPatNum == 0
         fields.TrialPatNum = 360;
    end
    PrefTrialIdx = findTrials(obj,fields)';
    fields.TrialPatNum = mod(angle+180,360);
    if fields.TrialPatNum == 0
         fields.TrialPatNum = 360;
    end
    PrefTrialIdx = sort([PrefTrialIdx findTrials(obj,fields)']);
    
    fields=[];
    fields.TrialSeqNum =4;
    fields.TrialPatNum = mod(angle+90,360);
    if fields.TrialPatNum == 0
         fields.TrialPatNum = 360;
    end
    AntiPrefTrialIdx = findTrials(obj,fields)';
    fields.TrialPatNum = mod(angle+270,360);
    if fields.TrialPatNum == 0
         fields.TrialPatNum = 360;
    end
    AntiPrefTrialIdx = sort([AntiPrefTrialIdx findTrials(obj,fields)']);

if nargin < 3 
    if isprop(obj,'Detrend') && obj.Detrend 
        detrend = 1;
    else
    detrend = 0;
    end
end
    
% Make a 1-D array indicating frames:
%  -1 : from anti-pref trials - will form 'background image'
%   0 : from non-selected trials - will be skipped
%   1 : from pref angle trials - will form 'activity image' 

t_on = zeros(1,size(obj.Frames,3));

for m = PrefTrialIdx
    t_on( obj.TrialStartFrame(m):obj.TrialEndFrame(m) ) = 1;
end
for m = AntiPrefTrialIdx
    t_on( obj.TrialStartFrame(m):obj.TrialEndFrame(m) ) = -1;
end
% % Set all frames in other trials to 0, these will be skipped.
% trialNull = setdiff( [1:length(obj.TrialStartFrame)] , trialIdx ) ;
% for n = trialNull
%     t_on( obj.TrialStartFrame(n):obj.TrialEndFrame(n) ) = 0;
% end

% % Find the F0 region before Trial1 of the pol rotate experiment
% % Find all trials with SeqNum == 4
% firstTrial = find([obj.TrialSeqNum == 4] ,1,'first');
% firstTrialFrame = obj.TrialStartFrame(firstTrial);
% % Now subtract find range { firstTrial - 6s : firstTrial - 2s }
% FPS = obj.AIrate/obj.IFI;
% 
%     
% Now set frames within this range to -1
% t_on( rangeStart:rangeEnd ) = -1;


if detrend
    Frames = detrendFrames(obj);
else
    Frames = obj.Frames;
end

% Get mean frame where this array is 1 (within selected trials)
avgTrialONframe = (mean(Frames(:,:,t_on==1),3));
% Get mean frame where the array is -1 (before LED turns on, before first pol rotate trial)
avgTrialOFFframe = (mean(Frames(:,:,t_on==-1),3));

%%% When no particular trial activity is sought, it's probably better to
%%% take the background image as the spontaneous activity before any
%%% stimuli run, as there can be activity between trials which distorts the
%%% subtraction 
% if  nargin < 2 || isempty(fields)
% avgTrialOFFframe = (mean(Frames(:,:,5:10),3));
% end

% Subtract one from the other and normalize
aFrame = (avgTrialONframe - avgTrialOFFframe);
% aFrame = avgTrialONframe ;
% And keep only positive values (negative value means less active
% during trials than before/after)
% aFrame(aFrame<0) = 0;

if nargin < 2
    % Push to object
    obj.ActivityFrame = aFrame;
else
    varargout{1} = aFrame;
end
