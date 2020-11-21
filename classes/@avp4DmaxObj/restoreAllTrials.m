function restoreAllTrials(obj)
if ~isempty(obj.storedTrialVariables)
    varStr{1} = 'TrialStartFrame';
    varStr{2} = 'TrialStartSample';
    varStr{3} = 'TrialEndFrame';
    varStr{4} = 'TrialEndSample';
    varStr{5} = 'TrialPatNum';
    varStr{6} = 'TrialSeqNum';
    
    if isprop(obj,'storedTrialVariables')
        origVar = obj.storedTrialVariables;
        
        
        % Restore original data
        for vidx = 1:length(varStr)
            obj.(varStr{vidx}) = origVar.(varStr{vidx});
        end
        obj.storedTrialVariables = [];
    end
    
end
end