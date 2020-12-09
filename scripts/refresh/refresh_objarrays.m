% Load original object arrays
pathsAVP
assert(exist(original_data_path,'dir')==7,'Requires original data folder (1 TB size): if folder is downloaded, update ''original_data_path'' in ''pathsAVP.m''')
for oidx = 1:length(objnames)
    thisobj = objnames{oidx};
    switch orig_path_mode
        case 0
            % copies objarrays from raw data path to minimal data storage
            copyfile(objpath_orig.(thisobj),objpath.(thisobj));
        case 1
            % refresh all values, autoROIs etc. and save back to raw data
            % path - run this first, then copy to standalone dir by
            % re-running with orig_path_mode = 0
            load(objpath_orig.(thisobj))
            superRefresh(x)    
            save(objpath_orig.(thisobj),'x','-v7.3')
    end
end
