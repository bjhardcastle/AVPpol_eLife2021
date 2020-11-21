function keepExpTrials(obj,expNum,specificTrials)
   % afterwards run      restoreAllTrials(obj,expNum)
        % cycleNum refers to exp4: cycle refers to half the conventional
        % polarizer cycle i.e. 30:180, not 30:360
        if nargin < 3 || isempty(specificTrials)
            specificTrials = 0; % ignore. Retain all trials specified by expNum
        end
        
        if isempty(expNum) && ~any(specificTrials==0)
           % in this mode, specificTrials are the absolute trial numbers.:
           % if expNum is specified, specificTrials will refer to the trial
           % position within that experiment only
           absTrialNum = 1;
        else
            assert(~isempty(expNum),'Must specify at least one input argument')
            absTrialNum = 0;
        end
        
      fields =[];
    fields.TrialSeqNum = expNum;
    trialKeepIdx = findTrials(obj,fields);
    if ~any(specificTrials==0)
        assert(length(specificTrials)<=length(trialKeepIdx),'specified trials exceed trials in this exp')
        if absTrialNum        
            assert(all(ismember(specificTrials,trialKeepIdx)),'specified trials not contained in this exp')
            % ignore trials just found with findTrials:
            trialKeepIdx = specificTrials;
        else
            assert(max(specificTrials)<=length(trialKeepIdx),'specified trials exceed trials in this exp')
            trialKeepIdx = trialKeepIdx(specificTrials);
        end
    end
    
    varStr{1} = 'TrialStartFrame';
    varStr{2} = 'TrialStartSample';
    varStr{3} = 'TrialEndFrame';
    varStr{4} = 'TrialEndSample';
    varStr{5} = 'TrialPatNum';
    varStr{6} = 'TrialSeqNum';
    
    % Store original data and keep only the requested exp trials
    for vidx = 1:length(varStr)
        origVar.(varStr{vidx}) = obj.(varStr{vidx});
        obj.(varStr{vidx}) = obj.(varStr{vidx})(trialKeepIdx);
        if expNum == 9 % also change the TrialPatNum for the blue flash trials, which is ambiguous
            obj.TrialPatNum(1:3) = 2;
        end
    end     
    
    obj.storedTrialVariables = origVar;
    
end