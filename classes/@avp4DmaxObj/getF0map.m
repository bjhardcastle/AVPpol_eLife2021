function F0map = getF0map(obj)
if any(obj.TrialSeqNum==4)
    expNum = 4; % experiments with polarizer
elseif any(obj.TrialSeqNum==2)
    expNum = 2; % experiments without polarizer
else
    disp('No pol tuning exp (2 or 4) available')
    return
end

% While we have the frames, also store the F0 values, from just before the
% experiment starts:
F0FrameArray = findF0LogicArray(obj,expNum);
% Extract mean response:
t0 = find(F0FrameArray,1,'first');
t1 = find(F0FrameArray,1,'last');
if isempty(obj.Frames)
    F0map = mean(getFrames(obj,[t0,t1]),3);
else
    F0map = mean(obj.Frames(:,:,t0:t1),3);
end
end


function F0FrameArray = findF0LogicArray(obj,expNum)
% 5sec recorded after each exp, 5sec recorded before each. Use 5s before. 
% We don't have a marker for where the 'set' starts, only the first trial.
% However, before each trial there is at least 1sec where the LED is on
% before the trial start marker, which we can use

lastFrame = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first')) - 1 - (floor(2.5*obj.AIrate/obj.IFI));
firstFrame = lastFrame - (floor(5*obj.AIrate/obj.IFI));
if firstFrame < 1
    firstFrame = 1;
end

% lastFrame = obj.TrialStartFrame(1);
% 
% firstFrame = lastFrame - (ceil(1*obj.IFI/obj.AIrate));%- (floor(1*obj.IFI/obj.AIrate));
% if firstFrame < 1
%     firstFrame = 1;
% end
% Construct logical array for F0 == 1
F0FrameArray = zeros(1,size(obj.Frames,3));
F0FrameArray(firstFrame:lastFrame)=1;

F0FrameArray = logical(F0FrameArray);
end
