% Get snapshots of tuning/selectivity data for ROIs in E-PGs (protocerebral bridge)
% for every cycle or every half cycle of the stimulus, then save as .mat
% files for looking at tuning over multiple cycles and variable polarotopy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Half cycle tuning snapshots 
clear x 
loadSS00096_PB
snapStruct = superCyclePolarHalfResp(x,1:8);
save(snapshotStruct_HalfCycle_path(savename),'snapStruct')
clear x % clear in case we have partial trial info remaining (See superCyclePolarHalfResp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single cycle tuning snapshots 
%{
clear x
loadSS00096_PB
snapStruct = superCyclePolarResp(x,1:8);
save(snapshotStruct_path(savename),'snapStruct')
clear x % clear in case we have partial trial info remaining (See superCyclePolarHalfResp)
%}