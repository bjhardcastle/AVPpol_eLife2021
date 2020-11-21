function [responseArray, timeVector, F0Array, ObjIdx,stimTimeArray] = findRespArrayExp4(objarray, ROImaskidx, fields)
%FINDRESPARRAY Extract the aligned, interpolated response time-series for
%all trials which match the specified parameters, across an object array
% [responseArray, timeVector, F0Array, ObjIdx] = findRespArray(objarray, ROImaskidx, fields)
%
% See documentation for findTrials. This is an extension of that function
% which accepts a [1xN] array of objects and returns all the trials, resampled
% to a common time series, timeVector, with correct duration in seconds.
% All trials are accumulated, so no information on which animal or Tiff
% they came from is retained, but the number of SlidebookObjs is counted
% and returned as numExps. If there are multiple tiffs per animal this
% could cause a problem for calculating SEM
%
% See also findTrials.

if nargin < 3 || isempty(fields)
    disp('''fields'' is empty. All trials will be returned.')
    fields = struct;
end
if nargin < 2 || isempty(ROImaskidx)
    ridx = 1;
    NOMASKFLAG = 1;
else
    ridx = ROImaskidx;
    NOMASKFLAG = 0;
end

% Get a common re-samplng rate and timeVector for plotting:
% [timeVector, normFPS, trialDuration] = findStandardTrial(objarray);
% Custom exp4 version of above:
t = [];
s = [];
inclObjIdx = [];
s0 = nan(length(objarray),1);
s1 = nan(length(objarray),1);
t0 = nan(length(objarray),1);
t1 = nan(length(objarray),1);
tStartFidx = [];sStartSidx = [];
tEndFidx = []; sEndSidx = [];

for oidx = 1:length(objarray)
    polExp = objarray(oidx).polExp;
if polExp >2 
     numPolTrials = 24;
else
     numPolTrials = objarray(oidx).pSet(polExp).trialReps*(360/objarray(oidx).pSet(polExp).polAngleStep);
     if numPolTrials > 24
         numPolTrials = 24;
     end
end
if  (ismember(polExp,objarray(oidx).TrialSeqNum)) && ...  % regular pol mapping exists
            (objarray(oidx).pSet(polExp).trialTestPauseLength == 4) && ...
            (objarray(oidx).pSet(polExp).trialMotorPauseLength == 1) && ...
            ( ~isfield(objarray(oidx).pSet(polExp),'StepDIR') || objarray(oidx).pSet(polExp).StepDIR ~= 5  ) && ...  % mapping was in the regular direction
        ( ~isfield(objarray(oidx).pSet(polExp),'polOffBetweenTrials') || objarray(oidx).pSet(polExp).polOffBetweenTrials == 0  )
        % (polarizer could not be reversed during early exps)
        
        inclObjIdx(end+1)= oidx;
        s0(oidx) =  objarray(oidx).TrialStartSample(find([objarray(oidx).TrialSeqNum == polExp],1,'first'));
        s1(oidx) =  objarray(oidx).TrialEndSample(find([objarray(oidx).TrialSeqNum == polExp],1,'first')+numPolTrials-1);
        
        t0(oidx) =  objarray(oidx).TrialStartFrame(find([objarray(oidx).TrialSeqNum == polExp],1,'first'));
        t1(oidx) = objarray(oidx).TrialEndFrame(find([objarray(oidx).TrialSeqNum == polExp],1,'first')+numPolTrials-1);
        
            sStartSidx(oidx,:) = objarray(oidx).TrialStartSample(find([objarray(oidx).TrialSeqNum == polExp],1,'first'):find([objarray(oidx).TrialSeqNum == polExp],1,'first')+numPolTrials-1);      
        sEndSidx(oidx,:)= objarray(oidx).TrialEndSample(find([objarray(oidx).TrialSeqNum == polExp],1,'first'):find([objarray(oidx).TrialSeqNum == polExp],1,'first')+numPolTrials-1);

        if ~isempty(s0)
            s(end+1) = ( s1(oidx) -s0(oidx) )/10000;
            t(end+1) = t1(oidx) -t0(oidx) ;        
            tStartFidx(end+1,:) = objarray(oidx).TrialStartFrame(find([objarray(oidx).TrialSeqNum == polExp],1,'first'):find([objarray(oidx).TrialSeqNum == polExp],1,'first')+numPolTrials-1);      
        tEndFidx(end+1,:)= objarray(oidx).TrialEndFrame(find([objarray(oidx).TrialSeqNum == polExp],1,'first'):find([objarray(oidx).TrialSeqNum == polExp],1,'first')+numPolTrials-1);

        end
    else 
        continue
    end
end
if isempty(inclObjIdx)
    responseArray = [];
F0Array = [];
stimTimeArray=[];
ObjIdx = [];
timeVector = [];
    return
end
trialDuration = round(nanmedian(s)); % seconds
maxTlength = nanmax(t); % frames
normFPS = 20*maxTlength / trialDuration;

% Get bounds for time vector from average frame-interval:
oneFrame = mean([objarray.IFI]./[objarray.AIrate]);

% Make standard trial time vector:
timeVector = [oneFrame:1/normFPS:trialDuration-oneFrame];


% Now cycle through objects in objarray, and store matching trials in
% 'trials_returned' (after interpolating):
tcount = 0;
ObjIdx = []; % Store the object index from which each trial was taken
trials_returned = [];
F0_returned = [];
for oidx = 1:length(objarray)
    
    if  NOMASKFLAG == 1 || ( ~isempty( objarray(oidx).ROI ) && ...
            length( objarray(oidx).ROI ) >= ridx && ...
            isfield( objarray(oidx).ROI(ridx), 'mask' ) && ...
            ~isempty( objarray(oidx).ROI(ridx).mask ) ) && ...
            ismember(oidx, inclObjIdx)
        
        assert( ~isempty( objarray(oidx).TrialStartFrame ) , [ 'No trial info stored. Run ''getTrialtimes(objarray(' num2str(oidx) ')) first'] );
        
        %trialidx = findTrials( objarray(oidx) , fields );
        %for tidx = trialidx'
        
        %Get ROI mask
        if NOMASKFLAG
            mask = [];
        else
            
            mask = objarray(oidx).ROI(ridx).mask;
        end
        
        % Get trial start / end frames
        %             secframes = ceil(objarray(oidx).AIrate ./ objarray(oidx).IFI);
        %             scan_extent = [objarray(oidx).TrialStartFrame(tidx) objarray(oidx).TrialEndFrame(tidx)];
        scan_extent = [t0(oidx)  t1(oidx)];
        
        % Get trial data within ROI
        this_trial = scanROI(objarray(oidx), mask, scan_extent);
        
        % Store trial, interpolated to common time vector
        %             trials_returned( tcount , : ) = findFittedTrial(objarray(oidx), this_trial, [], timeVector, normFPS, trialDuration);
        trial_samp_idx = [objarray(oidx).Frametimes( t0(oidx) : t1(oidx) )];
        
        % Convert into samples relative to trial start
        trial_samp_idx_rel = trial_samp_idx - s0(oidx);
        % Convert into times, in seconds, relative to trial start
        this_time = trial_samp_idx_rel./objarray(oidx).AIrate;
        
        tcount = tcount + 1;
        % Then interpolate the current trial data to fit the time points in the
        % standard trial (timeVector)
        trials_returned( tcount , : ) = interp1( this_time , this_trial , timeVector );
        
        % Also get the times of stimulus transitions
        stim_samp_idx = sort([ sStartSidx(oidx,:) sEndSidx(oidx,:) ],'ascend');
        % Convert into samples relative to trial start
        stim_samp_idx_rel = stim_samp_idx - stim_samp_idx(1);
        % Convert into times, in seconds, relative to trial start
        stim_times( tcount , : ) = stim_samp_idx_rel./objarray(oidx).AIrate;

    
    
    % Store F0 value for the trial
    F0_returned( tcount, : ) = findF0Exp4( objarray(oidx), mask, trials_returned( tcount , : ),[], timeVector, normFPS  );
    
    % Store the index of the object for the trial
    ObjIdx(tcount) = oidx;
    
    %end
    
    %if ~isempty(trialidx)
    
    if objarray(oidx).BackgroundSubtracted
        % Warn if background has been subtracted
        disp(['Object(' num2str(oidx) '): background subtraction has previously been applied']);
    end
    
    end
    

%     Otherwise, specified ROI mask does not exist for this object
end

responseArray = trials_returned;
F0Array = F0_returned;
stimTimeArray=stim_times;
if isempty( trials_returned )
    disp(['No ROI exists or no matching trials found for, ROI(' num2str(ridx) ')' ]);
end


