function getThrustMaps(obj)
% For 4D avp Layer or MIP objects. Find each pixel's bar response and position
% preference (exp5) and store as an array in the object, along with some
% average activity frames. Masks can be applied later, ROIs can be created
% later.
clearFramesFlag=0;
if isempty(obj.Frames)
    getFrames(obj)
    clearFramesFlag = 1;
end

detrendFlag = 0;
% Detrend Frames:
if detrendFlag
    Frames = detrendFrames(obj);
else
    Frames = obj.Frames;
end
if isempty(obj.TrialSeqNum)
    getParameters(obj)
end

expNum = 6;
if ~any(obj.TrialSeqNum==expNum)   
    disp('No thrust experiment (exp6) available')
    return
end

% For this experiment find start/stop frames:
barStart = obj.expStart(expNum);
barStop = obj.expStop(expNum);
barPos = unique(obj.TrialPatNum(obj.TrialSeqNum==expNum)); % This should always be the same
avgTrialLength = median(diff([obj.TrialStartFrame(obj.TrialSeqNum == expNum); obj.TrialEndFrame(obj.TrialSeqNum == expNum)]));

avgImg = mean(Frames(:,:,barStart:barStop),3);
stdImg = std(Frames(:,:,barStart:barStop),[],3);
%% Make avg & max activity maps for each bar position
% We make an activity map for each bar position, which includes the frames
% where the bar was presented plus all following frames until the next
% trial or the end of the experiment

% Find all trials which included bar presentation 
fields=[];
fields.TrialSeqNum = expNum;
allTrialIdx = findTrials(obj,fields)';

% For each bar position
for pIdx = 1:length(barPos)
    
    % Find trials which included bar presentation at this position
    fields=[];
    fields.TrialSeqNum = expNum;
    fields.TrialPatNum = barPos(pIdx);
    
    posTrialIdx = findTrials(obj,fields)';
    
    % Make a 1-D array indicating frames:
    %   0 : from all other frames - will be skipped
    %  >1 : from bar pos trial frames - will form 'within-trial activity image'
    %  -1 : from other bar trial frames - will form 'extra-trial activity image'
    t_on = zeros(1,size(obj.Frames,3));
    max_on = [];
    mIdx = 0;
    for m = posTrialIdx
        
        % Set marker for trial frames 
        t_on( obj.TrialStartFrame(m):obj.TrialEndFrame(m) ) = m;
        
        % Add post-trial sequence
        if m < posTrialIdx(end)
            t_on( obj.TrialEndFrame(m)+1 : obj.TrialStartFrame(m+1)-1 ) = m;
        else
            t_on( obj.TrialEndFrame(m)+1 : obj.TrialEndFrame(m) +avgTrialLength ) = m;
        end
        
        % Get the max activity for each trial at this bar position, so we
        % can average across repetitions
        mIdx= mIdx+1;
        max_on(:,:,mIdx) = max(Frames(:,:,[t_on==m]),[],3);
    end
    
    avgPosMax(:,:,pIdx) = mean(max_on,3);

    %{    
    %% Simplified sept2019: this is all unnecessary 
            
    otherTrialIdx = setdiff(allTrialIdx,posTrialIdx);

        for n = otherTrialIdx
            % Set marker for all other bar trials
            t_on( obj.TrialStartFrame(n):obj.TrialEndFrame(n) ) = -1;
            % Add post-trial sequence
            if otherTrialIdx < n
                t_on( obj.TrialEndFrame(n)+1 : obj.TrialStartFrame(n+1)-1 ) = -1;
            else
                t_on( obj.TrialEndFrame(n)+1 : obj.TrialEndFrame(n) + avgTrialLength ) = -1;
            end
        end
    
    % Get mean frame where this array is 1 (within selected trials)
    maxTrialONframe = max(Frames(:,:,t_on==1),[],3);
    % Get mean frame where the array is -1 (before LED turns on, before first pol rotate trial)
    maxTrialOFFframe = max(Frames(:,:,t_on==0),[],3);
    
    % Get mean frame where this array is 1 (within selected trials)
    avgPostONframe = mean(Frames(:,:,p_on==1),3);
    % Get mean frame where the array is -1 (before LED turns on, before first pol rotate trial)
    avgPostOFFframe = mean(Frames(:,:,p_on==-1),3);
    
    avgPostOFFframe2 = mean(Frames(:,:,p_on==1),3);
    avgPostONframe2 = mean(Frames(:,:,t_on==1),3);
    
    avgOtherTrialframe = mean(Frames(:,:,t_on==-1),3);
    maxOtherTrialframe  = max(Frames(:,:,t_on==-1),[],3);

    
    %%% When no particular trial activity is sought, it's probably better to
    %%% take the background image as the spontaneous activity before any
    %%% stimuli run, as there can be activity between trials which distorts the
    %%% subtraction
    % if  nargin < 2 || isempty(fields)
    % avgTrialOFFframe = (mean(Frames(:,:,5:10),3));
    % end
  
    % Subtract one from the other
%     aFrame{aIdx} = (maxTrialONframe - maxTrialOFFframe);
 aFrame{bIdx} =  avgPrefMax.*(avgPrefMax>( mean(Frames(:,:,t_on==-1),3) + 4*std(Frames(:,:,t_on==-1),[],3) ));
    pFrame{bIdx} = (avgPostONframe - avgPostOFFframe);
    
    bFrame{bIdx} = (avgPostONframe2 - avgPostOFFframe2);
   
    maxFrame(:,:,bIdx) = avgPrefMax;
    
    se = strel('square',2);
    binaryFrame{bIdx} = imfill(imclose(aFrame{bIdx}>0,se),'holes')  - imfill(imclose(pFrame{bIdx}>3*(std(std([pFrame{:}]))),se),'holes');
  %}
end

%% Push to object
obj.thrustMaxImg = avgPosMax;
obj.thrustAvgImg = avgImg;
obj.thrustStdImg = stdImg;

if clearFramesFlag
    obj.Frames = [];
    obj.Daq = [];
end
end

