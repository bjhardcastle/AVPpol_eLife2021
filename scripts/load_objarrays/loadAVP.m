
switch orig_path_mode
    
    case 0 % Use backed up minimal data versions
        
        load(objpath.(savename))
        
        % switch to standalone data directory
        superChangeRoot(x,original_data_path,aux_backup_path)
        
    case 1 % Use raw data versions
        
        disp('orig_path_mode is On: loading objarrays from raw data storage')
        load(objpath_orig.(savename))
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Common commands to run after loading objarray

% print array name s
fprintf(['\n' savename ,' object array loaded as variable ''x''\n\n' ])

loadLayerMasks([x.MIP])
loadROIs([x.MIP])

% clear leftover plotting-specific variables from a previous object
clear prefix suffix selectObj selectLayer

% find index of recordings made on right or left side
inclR = find(cellfun('isempty',strfind({x.Name},' L')));
inclL = find(~cellfun('isempty',strfind({x.Name},' L')));

