function getPolMapsUnpooled(obj)
% For 4D avp Layer or MIP objects. Find each pixel's polarization tuning
% (exp2 w/o polarizer, exp4 with polarizer) and store as an array in the
% object, along with some average activity frames. Masks can be applied
% later, ROIs can be created later.
clearFramesFlag=0;
if isempty(obj.Frames)
    getFrames(obj)
    clearFramesFlag = 1;
end

detrendFlag = 0;
% Detrend Frames:
if detrendFlag
    detrendedFrames = detrendFrames(obj);
else
    detrendedFrames = obj.Frames;
end
if isempty(obj.TrialSeqNum)
    getParameters(obj)
end

if any(obj.TrialSeqNum==4)
    expNum = 4; % experiments with polarizer
elseif any(obj.TrialSeqNum==8)
    expNum = 8; % experiments with polarizer (fine steps)
elseif any(obj.TrialSeqNum==2)
    expNum = 2; % experiments without polarizer
else
    disp('No pol tuning exp (2,4 or 8) available')
    return
end

% For this experiment (whether pol attached or not) find start/stop frames:
polStart = obj.expStart(expNum);
polStop = obj.expStop(expNum);
polAngs = obj.expPolAngs(expNum);

f = median(detrendedFrames(:,:,polStart:polStop),3);


% Get indices for frame size
indLength = size(obj.Frames,1)*size(obj.Frames,2);
rows = reshape(repmat(1:size(obj.Frames,1),1,size(obj.Frames,2)),indLength,1);
cols = reshape(repmat(1:size(obj.Frames,2),size(obj.Frames,1),1),indLength,1);
inds = sub2ind(size(obj.Frames(:,:,1)), rows, cols);

% Extract every pixel's time-series (whole experiment)
resp = reshape(detrendedFrames, indLength, size(obj.Frames,3));
detrendedFrames = [];

% Extract every pixel's time-series (pol tuning experiment only)
% Detrend Frames:
if detrendFlag
    respPolExp = detrend(resp(inds,polStart:polStop)')';
else
    respPolExp = resp(inds,polStart:polStop);
end

%% FFT stuff
% Find FFT of every pixel's time-series. Old method used phase to determine
% tuning for each pixel, but doesn't work well with the discontinous
% stimulus. Now we just use the FFT magnitude as a way of thresholding
% responsiveness to the stimulus (capturing any pixels whose activity
% varies at 2x the polarizer rotation frequency) 
L = size(respPolExp,2);
n = 2^nextpow2(L); % sampling of FFT
Y = fft(respPolExp,n,2); % FFT array [pixel index , frequency index]
Fs = obj.AIrate/obj.IFI;
F = Fs*(0:(n/2))/n;
P = abs(Y/n);

% Now work out the frequency of responses caused by polarizer:
% (one cycle of polarizer rotation = two cycles of activity in
% pol-sensitive cells)
numHalfCycles = 2*length(find(obj.TrialPatNum(obj.TrialSeqNum == expNum) == 360));
halfCycFreq = 1/((1/numHalfCycles)*(polStop-polStart+1)/Fs);

% Find corresponding FFT frequency:
[~,fIdx]=min(abs(halfCycFreq - F));

% Array of fft magnitudes for each pixel (useful for eliminating pixels in
% padded frame edges which have artifically high pol selectivity (from
% registration), but their time-series don't follow a sinusoidal profile):
fftMagImg = reshape(sqrt(abs(Y(:,fIdx))),size(f));
obj.fftMagImg = fftMagImg./100;

% Find the pref. pol angle for each pixel:
%   +ve phase of FFT means shifted to 'earlier' angle, which is
%   counterintuitive when angle usually increases with time, so mulitply by -1
%   The pol stimulus might not start at zero, so add the first angle
%   (before division by 2, so must be doubled)
%   Then convert -90:90 to 0:180 with wrapTo360
%   Phase angle from FFT is in range -pi:pi, which maps to -90:90 of the
%   pol stimulus, so we must divide phase angle by 2
angPref = wrapTo360(polAngs(1)*2 + rad2deg(-angle(Y(:,fIdx))))./2;

fftAngImg = reshape(angPref,size(f));
obj.fftAngImg = fftAngImg;


%% Preferred angle map
% Here we get the average response of every pixel for each angle presented.
% From this we can find the pol selectivity: (max - min)/max
% which is a good parameter for excl./including pixels

% Construct array [numPix x numPolAngles/2], into which we'll enter the
% average response for each pol angle, for every pixel:
pixPolResp = nan(indLength,length(polAngs));
% For each angle:
for aIdx = 1:length(polAngs)
    % Get logic array that indicates every frame that was recorded at this
    % angle (without pooling with angle+180 frames)
    angleFrameArray = getAngleLogicArray2(obj,polAngs(aIdx));
    
    % Extract mean response for this angle:
        pixPolResp(:,aIdx) = mean(respPolExp(:,angleFrameArray(polStart:polStop)),2);
%         pixPolResp(:,aIdx) = max(respPolExp(:,angleFrameArray(polStart:polStop)),[],2);

end
% Assign to object:
obj.polAngResp = pixPolResp;

% Find circular mean for each pixel:
pixPolWeight = pixPolResp; 
pixPolWeight(isnan(pixPolWeight)) = 0;
[r,axtheta] = circ_rangle(circ_axial(deg2rad(repmat(polAngs(1:length(polAngs)),size(pixPolResp,1),1)),2), pixPolWeight, deg2rad(mode(diff(polAngs))), 2);
theta = axtheta;
% [r,axtheta] = circ_rangle( deg2rad( repmat(polAngs(1:length(polAngs)),size(pixPolResp,1),1) ) , pixPolWeight, deg2rad(mode(diff(polAngs))), 2);
% theta = axtheta;
Ang =(wrapTo360(((rad2deg(theta)))));
obj.polTuningImg = reshape(Ang,size(f));
obj.polVecLength = reshape(r,size(f));
%% Selectivity map
% Discretize the preferred angle values to find stimulus which gave mean
% pref resp, then find response at orthogonal angles:

% If polAngs = [30 60 90], make bin edges [15 45 75 105]
angEdges = polAngs(1:0.5*length(polAngs)+1)-0.5*(mode(diff(polAngs)))  ;
% deal with 0/180 and make bins centered on each angle :
shiftAng =0.5.*(wrapTo360((2.*(rad2deg(theta)-0.5*mode(diff(polAngs))))))  + 0.5*mode(diff(polAngs));

maxIdx = discretize(shiftAng,angEdges);  % fft mag threshold applied

% For each orientation:
pixPolResp = [];
for aIdx = 1:length(polAngs)
    % Get logic array that indicates every frame that was recorded at the
    % angle (don't pool angle+180 frames)
    angleFrameArray = getAngleLogicArray2(obj,polAngs(aIdx));
    
    % Extract mean response for this angle:
    pixPolResp(:,aIdx) = mean(respPolExp(:,angleFrameArray(polStart:polStop)),2);
end
% Assign to object:
obj.polAngRespPooled = pixPolResp;

respMax = pixPolResp( sub2ind(size(pixPolResp), [1:indLength]', maxIdx ) );

% From the max angle, find the orthogonal angle, which will be used as the
% minimum response:
maxIdx = maxIdx - 1; % To allow 'mod' function below to work. Add 1 again after using.
respMin = pixPolResp( sub2ind(size(pixPolResp), [1:indLength]', 1+mod(maxIdx + 0.25*length(polAngs),0.5*length(polAngs)) ) );
maxIdx = maxIdx + 1;

% Get the pol selectivity vector:
polSelVec = ( respMax - respMin )./(respMax + respMin);
% Raw image intensity values cannot be negative so selectivity cannot be >1
% and by definition (max - min)/max cannot be less than 0 (if max
% orientation was detected correctly), so clip to give an values
% between 0:1
polSelVec(polSelVec<0) = 0;
polSelVec(polSelVec>1) = 1;
polSelVec(isnan(polSelVec)) = 0;
% Reshape to form image:
polSelImg = reshape(polSelVec,size(f));
% obj.polSelImg = polSelImg;

obj.polAngImg = reshape(maxIdx,size(f)).*mode(diff(polAngs));

obj.avgPolImg = f;
obj.varImg = std(obj.Frames,[],3).*mean(obj.Frames,3);
% obj.polActivityImg = g; 

if clearFramesFlag
    obj.Frames = [];
    obj.Daq = [];
end
end

function angleFrameArray = getAngleLogicArray2(obj,polAng)
    % Return a logic array (1 x numFrames) indicating where the specified
    % angle was presented (alternative version below included angle+180)
    
% Extract the start:stop frame numbers for the pol tuning experiment:
if any(obj.TrialSeqNum==4)
    expNum = 4;
elseif any(obj.TrialSeqNum==8)
    expNum = 8; % experiments with polarizer (fine steps)
elseif any(obj.TrialSeqNum==2)
    expNum = 2;
else
    disp('No pol tuning exp (2 or 4) available')    
    return
end

polStart = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first'));
polStop = obj.TrialEndFrame(find(obj.TrialSeqNum == expNum,1,'last'));


% Find the frames recorded with polarizer at polAng:
angStarts = obj.TrialStartFrame(obj.TrialPatNum == polAng);
angStops = obj.TrialEndFrame(obj.TrialPatNum == polAng);

% Don't add orthogonal angle frames 

% Construct logical array for frameAngles == polAng
angleFrameArray = zeros(1,size(obj.Frames,3));
for n = 1:length(angStarts)
    if (angStarts(n) >= polStart) && (angStops(n) <= polStop)
        angleFrameArray(angStarts(n):angStops(n)) = 1;
    end
end

angleFrameArray = logical(angleFrameArray);
end

function angleFrameArray = getAngleLogicArray(obj,polAng)
  % Return a logic array (1 x numFrames) indicating where the specified
    % angle was presented AND where angle+180 was presented
% Extract the start:stop frame numbers for the pol tuning experiment:

if any(obj.TrialSeqNum==4)
    expNum = 4;
elseif any(obj.TrialSeqNum==8)
    expNum = 8; % experiments with polarizer (fine steps)
elseif any(obj.TrialSeqNum==2)
    expNum = 2;
else
    disp('No pol tuning exp (2 or 4) available')
    return
end

polStart = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first'));
polStop = obj.TrialEndFrame(find(obj.TrialSeqNum == expNum,1,'last'));


% Find the frames recorded with polarizer at polAng:
angStarts = obj.TrialStartFrame(obj.TrialPatNum == polAng);
angStops = obj.TrialEndFrame(obj.TrialPatNum == polAng);

% Add the equivalent +180deg angle
angStarts = [angStarts, obj.TrialStartFrame(obj.TrialPatNum == polAng + 180)];
angStops = [angStops, obj.TrialEndFrame(obj.TrialPatNum == polAng + 180)];

% Construct logical array for frameAngles == polAng
angleFrameArray = zeros(1,size(obj.Frames,3));
for n = 1:length(angStarts)
    if (angStarts(n) >= polStart) && (angStops(n) <= polStop)
        angleFrameArray(angStarts(n):angStops(n)) = 1;
    end
end

angleFrameArray = logical(angleFrameArray);
end


function F0FrameArray = findF0LogicArray(obj,expNum)
% 5sec recorded after each exp, 5sec recorded before each. Use 5s before. 
% We don't have a marker for where the 'set' starts, only the first trial.
% However, before each trial there is at least 1sec where the LED is on
% before the trial start marker, which we can use

lastFrame = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first')) - 1 - (floor(1*obj.IFI/obj.AIrate));
firstFrame = lastFrame - (floor(5*obj.IFI/obj.AIrate));
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

