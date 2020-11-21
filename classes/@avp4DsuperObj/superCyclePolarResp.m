function snapStruct = superCyclePolarResp(objarray,roiIdx)
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
        
        % get all ROIs if not specified
        if nargin < 2 || isempty(roiIdx)
           roiIdx = 1:length(objarray(oidx).MIP.ROI); 
        end

        %  get number of half cycles in experiment
        trialsPerCycle = 360/objarray(oidx).MIP.pSet(4).polAngleStep;
        if trialsPerCycle ~= 12
            disp(['different angle step size to regular experiment (' num2str(oidx) ')'])
        end
        
                    
        if includeMultipleExpsPerRecording
            inclTrials = (find( objarray(oidx).MIP.TrialSeqNum==4 | objarray(oidx).MIP.TrialSeqNum==10 ));
            numTrialsOrig = length(inclTrials);
            numCycles = numTrialsOrig/trialsPerCycle;
            inclExps = inclTrials(1:trialsPerCycle:end);

        elseif useFirstFourCyclesOnly
            % only look at first four cycles in first exp4 in recording
            numTrialsOrig = 2*trialsPerCycle;
            numCycles = 2;
            inclExps = 4*ones(1,numCycles);
            inclTrials = (find( objarray(oidx).MIP.TrialSeqNum==4 , numTrialsOrig, 'first'));

        else
            % only look at first exp4 in recording
            numTrialsOrig = objarray(oidx).MIP.pSet(4).trialReps*trialsPerCycle;
            numCycles = objarray(oidx).MIP.pSet(4).trialReps;
            inclExps = 4*ones(1,numCycles);
            inclTrials = (find( objarray(oidx).MIP.TrialSeqNum==4 ));
        end
        
        trialArray = reshape( inclTrials, trialsPerCycle, numCycles)';

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
            
            for caseIdx = 1
                
                
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
                    keepExpTrials(objarray(oidx).MIP,[],sort(trialsToKeep(:)))
                    objarray(oidx).MIP.limitCycles = 1;
                    getPolMaps(objarray(oidx).MIP)
                    objarray(oidx).MIP.UseFixedResp = 0;
                    [ROImeanTuningAng,ROIresultantLength] = plotSinglePolarResp(objarray(oidx),roiIdx,1);
                    close gcf
                    
                    % store angle tunings
                    snapStruct(oidx).(caseStr{caseIdx}).ROImeanTuningAng(cycIdx,:) = ROImeanTuningAng;
                    snapStruct(oidx).(caseStr{caseIdx}).ROIresultantLength(cycIdx,:) = ROIresultantLength;
                    for rIdx = roiIdx
                       snapStruct(oidx).(caseStr{caseIdx}).ROIpsi(cycIdx,rIdx) = mean(objarray(oidx).MIP.polSelImg(objarray(oidx).MIP.ROI(rIdx).mask)); 
                    end
                    snapStruct(oidx).(caseStr{caseIdx}).bkgpsi(cycIdx) = mean(objarray(oidx).MIP.polSelImg(objarray(oidx).MIP.bkgMask));
                    snapStruct(oidx).(caseStr{caseIdx}).bkgpsiStd(cycIdx) = std(objarray(oidx).MIP.polSelImg(objarray(oidx).MIP.bkgMask));

                    snapStruct(oidx).(caseStr{caseIdx}).Fly = objarray(oidx).Fly;
                    snapStruct(oidx).(caseStr{caseIdx}).ID = str2double([ objarray(oidx).DateStr objarray(oidx).TimeStr]);

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