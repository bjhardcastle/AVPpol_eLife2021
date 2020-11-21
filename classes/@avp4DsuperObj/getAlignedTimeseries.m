function [timeVector,fittedData] = getAlignedTimeseries(obj,ROI,expNum,inclPre,inclPost)
% obj    avp4dsuperObj, not a subobj
% ROI    any ROI structure containing at least .mask OR just a mask 
%        This will refer to one of the subobjects
% inclPre (0 default, or 1) to include frames at the beginning of the
% experiment, but  before the first trial 
% inclPost (0 default, or 1) to include frames at the end of the
% experiment, after the last trial 
assert(nargin>2,'Must input expNum and ROI') 
assert(~isempty(expNum) && isnumeric(expNum), 'Must specify expNum as a numerical value')
assert(~isempty(ROI), 'Must provide an ROI strucutre or mask')
if nargin < 4 || isempty(inclPre)
    inclPre = 0;
end
if nargin < 5 ||isempty(inclPost)
    inclPost = 0;
end
% Deal with ROI input which can be either a mask or a regular ROI struct
if isstruct(ROI)
    if ~isfield(ROI,'mask') || isempty(ROI.mask)
       disp('ROI mask is missing')
       return
    end
elseif isnumeric(ROI) && length(unique(ROI))==2 % should be 0s and 1s only, but may not be logical
    mask = ROI;
    clear ROI
    ROI = struct;
    ROI.mask = mask;
   
    if isempty(obj.MIP.Frames)
        getFrames(obj.MIP)
    end
    ROI.response = scanROI(obj.MIP,mask);
end
    

% Upsampling rate:
Fs = 20;

switch expNum
    
    case {2,4}
        % % Exp2/4
        %
        % % Before exp :
        % obj.pSet(expNum).recPreExpPauseLength % lower bound
        %
        % % Before first trial:
        % % LED ON (at zero deg)
        %
        % % % Rotate polarizer:
        % % obj.pSet(expNum).trialMotorPauseLength
        % %
        % % %Trial start
        % % obj.pSet(expNum).trialTestPauseLength
        % % % Trial end
        %
        % After exp:
        % pause(pSet.recPostExpPauseLength)

        numAngleTrials = obj.pSet(expNum).trialReps*360/obj.pSet(expNum).polAngleStep;
        %TrialStart(1): TrialStop(end) (in number of Frames)
        numSec = obj.pSet(expNum).trialMotorPauseLength*(numAngleTrials-1)  + obj.pSet(expNum).trialTestPauseLength*numAngleTrials;
        % Typically: 119s + delay incurred by motor rotation on each trial
        
        preSec =  obj.pSet(expNum).trialMotorPauseLength;
        postSec = obj.pSet(expNum).recPostExpPauseLength;
        
        
        R = ( obj.MIP.Frametimes(obj.MIP.expStop(expNum))  - obj.MIP.Frametimes(obj.MIP.expStart(expNum)) )/obj.MIP.AIrate/numSec;

    case {1,3}
        % Exp1/3
        %
        % Before exp:
        % pause(pSet.recPreExpPauseLength)  % lower bound
        %
        % Between polAngles (pSet.polAngleArray) (not before first)
        % pause(pSet.trialPreRotatePauseLength)
        %
        % pause(pSet.trialBaselinePauseLength)
        % Trial starts, LED ON
        %  pause(pSet.trialTestPauseLength)
        
        % After exp:
        % pause(pSet.recPostExpPauseLength)
        
        numTrials = length( obj.pSet(expNum).polAngleArray )*obj.pSet(expNum).trialReps;
        numInterTrialPauses = length( obj.pSet(expNum).polAngleArray ) -1;
        
        numSec = numTrials*( obj.pSet(expNum).trialTestPauseLength + obj.pSet(expNum).trialBaselinePauseLength ) ...
            + numInterTrialPauses*obj.pSet(expNum).trialPreRotatePauseLength;
        
        
        preSec = obj.pSet(expNum).recPreExpPauseLength;
        postSec = obj.pSet(expNum).recPostExpPauseLength;

             R = ( obj.MIP.Frametimes(obj.MIP.expStop(expNum))  - obj.MIP.Frametimes(obj.MIP.expStart(expNum)) )/obj.MIP.AIrate/numSec;
   
end



% Get the sample index of each frame within the exp
trial_samp_idx = [obj.MIP.Frametimes( obj.MIP.expStart(expNum) : obj.MIP.expStop(expNum) )];

% Convert into samples relative to start of first trial in exp
trial_samp_idx_rel = trial_samp_idx - obj.MIP.TrialStartSample(find(obj.MIP.TrialSeqNum==expNum,1,'first'));

% Convert into times, in seconds, relative to trial start
exp_time = trial_samp_idx_rel./(obj.MIP.AIrate*R);

expTimeVector = linspace(0,numSec- 1/Fs,Fs*numSec );

if inclPre
    % Get the sample index of each frame within the 'pre exp' section
    pre_samp_idx = [obj.MIP.Frametimes( floor(obj.MIP.expStart(expNum) -1 - preSec*R ) : obj.MIP.expStart(expNum)-1  )];
    % Convert into samples relative to start of first trial in exp
    pre_samp_idx_rel = pre_samp_idx - obj.MIP.TrialStartSample(find(obj.MIP.TrialSeqNum==expNum,1,'first'));

    % Convert into times, in seconds, relative to trial start
    pre_time = pre_samp_idx_rel./(obj.MIP.AIrate*R);
 
    preTimeVector = linspace(-preSec,-1/Fs,Fs*preSec );
    
end
if inclPost
    % Get the sample index of each frame within the 'pre exp' section
    post_samp_idx = [obj.MIP.Frametimes( ceil(obj.MIP.expStop(expNum) + 1 : obj.MIP.expStop(expNum) + 1 + postSec*R ) )];
    % Convert into samples relative to start of first trial in exp
    post_samp_idx_rel = post_samp_idx - obj.MIP.TrialStartSample(find(obj.MIP.TrialSeqNum==expNum,1,'first'));

    % Convert into times, in seconds, relative to trial start
    post_time = post_samp_idx_rel./(obj.MIP.AIrate*R);
 
    postTimeVector = linspace(numSec,numSec+postSec,Fs*postSec );
    
end


expData = ROI.response(obj.MIP.expStart(expNum): obj.MIP.expStop(expNum));
preData = ROI.response( floor(obj.MIP.expStart(expNum) -1 - preSec*R ) : obj.MIP.expStart(expNum)-1  );
postData = ROI.response(  ceil(obj.MIP.expStop(expNum) + 1 : obj.MIP.expStop(expNum) + 1 + postSec*R ) );

if ~inclPre && ~inclPost
    timeVector = expTimeVector;
    fittedData = interp1( exp_time , expData , timeVector);
elseif inclPre && ~inclPost
    timeVector = [preTimeVector expTimeVector];
    fittedData = interp1( [pre_time exp_time ], [ preData expData], timeVector);
elseif ~inclPre && inclPost
    timeVector =  [expTimeVector postTimeVector];
    fittedData = interp1( [ exp_time post_time], [  expData postData], timeVector);
elseif inclPre && inclPost
    timeVector = [preTimeVector expTimeVector postTimeVector];
    fittedData = interp1( [pre_time exp_time post_time], [ preData expData postData], timeVector);
end

end




