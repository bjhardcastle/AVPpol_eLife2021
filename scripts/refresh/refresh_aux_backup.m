% Copy small .mat files including ROIs, layer masks etc. from the original
% tiff data folders, to include in the standalone data package
pathsAVP

assert(exist(original_data_path,'dir')==7,'Requires original data folder (1 TB size): if folder is downloaded, update ''original_data_path'' in ''pathsAVP.m''')

if exist(aux_backup_path,'dir')
    rmdir(aux_backup_path,'s')
end

for ridx = 1:length(objnames)
    objname = objnames{ridx}; % object array name
    disp(['Copying ' objname ' data'])
    
%     if strcmp(objname,'R88A06_Bu_ant') % corresponds to the same folder as 'R88A06_Bu'
%         continue
%     else
        % load object array (without using loadAVP func that potentially
        % changes object folder root)
       load(objpath.(objname))

        for oidx = 1:length(x) % for each object included in analysis
            
            thisID = [x(oidx).DateStr,x(oidx).TimeStr]; % recording identifier used in some filenames
            
            objroot_orig = fullfile(original_data_path, rel_data_paths.(objname), x(oidx).Name);
            if exist(objroot_orig,'dir')
                [~,objfolder] = fileparts(objroot_orig);
            else
                objfolder = ls([fullfile(original_data_path, rel_data_paths.(objname), x(oidx).Name) '*']);            
                objroot_orig = fullfile(original_data_path, rel_data_paths.(objname),objfolder);
            end
            objroot_copy = fullfile(aux_backup_path, rel_data_paths.(objname), objfolder);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % copy log file
            % copyfile(fullfile(objroot_orig,x(oidx).LogFile),fullfile(objroot_copy,x(oidx).LogFile))
            copyfile(fullfile(objroot_orig,[x(oidx).LogFile,'*']),fullfile(objroot_copy))

            folderStr = {'layers';'projection_max'};
            subfolderStr = {'C0';'C1'};
            for fidx = 1:2
                for sidx = 1:2
                    thisdir = fullfile(objroot_orig,folderStr{fidx},subfolderStr{sidx});
                    copydir = fullfile(objroot_copy,folderStr{fidx},subfolderStr{sidx});
                    if exist(thisdir,'dir')
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % copy MSP file
                        thismsp = fullfile(thisdir,[thisID,'_MSP.mat']);
                        if ~strcmp(folderStr{fidx},'layers') && exist(thismsp,'file')
                            copyfile([thismsp '*'],copydir)
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % copy ROIs
                        thisROI = fullfile(thisdir,[thisID,'ObjROIs.mat']);
                        if exist(thisROI,'file')
                            copyfile([thisROI '*'],copydir)
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % copy polROIs
                        thisPolROI = fullfile(thisdir,[thisID,'PolROIs.mat']);
                        if exist(thisPolROI,'file')
                            copyfile([thisPolROI '*'],copydir)
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % copy layer masks
                        switch fidx
                            case 1
                                layermaskfiles = ls(fullfile(thisdir,[thisID '*_Z*_layerMask.mat']));
                            case 2
                                layermaskfiles = ls(fullfile(thisdir,[thisID '_layerMask.mat']));
                        end
                        for midx = 1:size(layermaskfiles,1)
                            copyfile( fullfile(thisdir, [strtrim(layermaskfiles(midx,:)) '*']), fullfile(copydir) )
                        end
                     
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % copy alt files
                        
                        if fidx == 2 && isprop(x(oidx).MIP,'altMaskStr') && ~isempty(x(oidx).MIP.altMaskStr)
                            altfiles = ls(fullfile(thisdir,[thisID '_' x(oidx).MIP.altMaskStr '_*.mat']));
                            for aidx = 1:size(altfiles,1)
                                copyfile( fullfile(thisdir, [strtrim(altfiles(aidx,:)) '*']), fullfile(copydir) )
                            end
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % copy PB glomerulus limits (auto-generated, then
                        % manually checked/adjusted)
                        
                        thisPBlim = fullfile(thisdir,[thisID,'PBglomLims.mat']);
                        if exist(thisPBlim,'file')
                            copyfile([thisPBlim '*'],copydir)
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    end
                end
            end
            
        end
%     end
end
disp('Done')