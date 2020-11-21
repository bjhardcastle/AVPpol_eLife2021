function refresh_psi()
%%% Previously 'refreshPSIstructs.m' in 18_data/plotting folder

% load super object arrays one by one and get their PSI struct containing
% values of pol sel index in cells and background, with and without the
% polarizer attached. 
% save all of these structures in a .mat file for quick access (loading
% each superobj array takes enough time that this is tedious to do whenever
% we want to compare PSI values in different celltypes)

pathsAVP 
PSI = struct();
for fieldIdx = 1:length(objnames)

    
    % Load superobjarray x
    clear x
    run(fullfile(load_objarray_func_path,['load' objnames{fieldIdx} '.m']))
    
    superUseMSP(x,1)
    superPolThreshold(x,-1)
    loadLayerMasks([x.MIP])
    
    % Get PSI data
    PSI.(objnames{fieldIdx}) = getPSIstruct(x,1);

end

save(psi_mat_path,'PSI','-v7.3')