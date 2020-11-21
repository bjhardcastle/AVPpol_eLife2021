function snapStruct = superCycleSnapshots(objarray)
% For checking cycle-by-cycle variation in selectivity index and
% tuning angle in all pixels. For each MIP object, we discard all but the
% trials in a certain cycle, then recalculate maps with getPolMaps.
% Arrays are stored in a structure for storage since it takes a while to
% run.
% (recycles a lot of code from superPBMIP)

useFirstFourCyclesOnly = 0; % PB recordings may have more - option to limit here for comparison to other areas
includeMultipleExpsPerRecording = 1;
if useFirstFourCyclesOnly
    includeMultipleExpsPerRecording = 0;
end

snapStruct = struct();

for oidx = 1:length(objarray)
    polSelImg = [];
    polTuningImg = [];
    if ( ~isempty(objarray(oidx).MIP.polExp) && objarray(oidx).MIP.polExp == 4 ) ... % only look at regular mapping exps (so tuning shifts are compared fairly across exps)
            && ( ~isfield(objarray(oidx).MIP.pSet(4),'polOffBetweenTrials') || objarray(oidx).MIP.pSet(4).polOffBetweenTrials == 0 ) ...
            && ( ~isfield(objarray(oidx).MIP.pSet(4),'trialRandomizeOrder') || objarray(oidx).MIP.pSet(4).trialRandomizeOrder == 0 )
        % also can't use randomized presentations since
        % we'll break everything down to half-cycles
        
        %  get number of half cycles in experiment
        trialsPerCycle = 180/objarray(oidx).MIP.pSet(4).polAngleStep;
        if trialsPerCycle ~= 6
            disp(['different angle step size to regular experiment (' num2str(oidx) ')'])
        end
        
        if includeMultipleExpsPerRecording
            numTrialsOrig = length(find(objarray(oidx).MIP.TrialSeqNum==4));
            numCycles = numTrialsOrig/trialsPerCycle;
        elseif useFirstFourCyclesOnly
            % only look at first four cycles in first exp4 in recording
            numTrialsOrig = 4*trialsPerCycle;
            numCycles = 4;
        else
            % only look at first exp4 in recording
            numTrialsOrig = objarray(oidx).MIP.pSet(4).trialReps*trialsPerCycle*2;
            numCycles = objarray(oidx).MIP.pSet(4).trialReps*2;
        end
        
        getFrames(objarray(oidx).MIP)
        
        % Now get MIPs with subsets of trials and remake pol maps
        % Progress from 1:numCycles and store snapshots at each stage
        %
        % We'll make four different types of subset:
        %
        % 1) trials in each individual cycle on their own, in
        % chronological order
        % 2) a cumulative set of trials, adding each new cycle to the
        % previous ones, progressing from start to finish in time
        % 3) cumulative set as above, but in reverse chronological order
        % 4) cumulative set as above, but in a shuffled order to eliminate
        % time effects
        %
        % NB:
        % fundamental object info is modified and must be restored with
        % restoreAllTrials()
        %
        % We overwrite the existing MIP object's trial info each time, and
        % must restore the original or we'll run into problems (Daq limits
        % will be incorrect for one). Try/catch handles function crashes,
        % make sure to run it manually if debugging
        caseStr = {'individual';'cumul_fwd_chrono';'cumul_rev_chrono';'cumul_shuffled'};
        try
            
            for caseIdx = 1:size(caseStr,1)
                
                trialArray = reshape( 1:numTrialsOrig, trialsPerCycle, numCycles)';
                
                for cycIdx = 1:size(trialArray,1) % number of cycles
                    
                    switch caseStr{caseIdx}
                        
                        case 'individual' % individual cycles, chronological order
                            
                            trialsToKeep = trialArray(cycIdx,:);
                            
                        case 'cumul_fwd_chrono' % cumulative, chronological
                            trialsToKeep = trialArray(1:cycIdx,:);
                            
                        case 'cumul_rev_chrono' % cumulative, reverse chronological
                            flipTrialArray = flipud(trialArray);
                            trialsToKeep = flipTrialArray(1:cycIdx,:);
                            
                        case 'cumul_shuffled'
                            % cumulative, shuffled
                            % seed with object identifier for repeatability
                            rng(str2double([objarray(oidx).DateStr objarray(oidx).TimeStr]))
                            shuffTrialArray = trialArray(randperm(numCycles),:);
                            trialsToKeep = shuffTrialArray(1:cycIdx,:);
                            
                    end
                    
                    % get subset
                    keepExpTrials(objarray(oidx).MIP,4,sort(trialsToKeep(:)))
                    objarray(oidx).MIP.limitCycles = 0;
                    getPolMapsHalfCycle(objarray(oidx).MIP)
                    
                    % store pol maps
                    snapStruct(oidx).(caseStr{caseIdx}).polSelImg(:,:,cycIdx) = objarray(oidx).MIP.polSelImg;
                    snapStruct(oidx).(caseStr{caseIdx}).polTuningImg(:,:,cycIdx) = objarray(oidx).MIP.polTuningImg;
                    
                    restoreAllTrials(objarray(oidx).MIP)
                    
                end
            end
            
            objarray(oidx).MIP.Frames=[];
            objarray(oidx).MIP.snapStruct = snapStruct(oidx);
        catch ME
            restoreAllTrials(objarray(oidx).MIP)
            objarray(oidx).MIP.Frames=[];
            rethrow(ME)
        end
    end
end
end