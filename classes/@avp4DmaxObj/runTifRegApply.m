function runTifRegApply(obj,regDispXYfullpath)
%RUNTIFREG Runs ca_RegSeries_v4 image alignment on an object's Tiff file
% runTifRegApply(obj,regDispXYfullpath)
%
% This function uses the displacement data in the text file stored at
% 'regDispXYfullpath' to apply an existing registration to another tiff
% file. Useful for applying registration of activity-independent images
% with a more stable image to activity-dependent images which are more
% variable over time. 
%
% See also ca_RegSeries_v4_Apply_regDispXY, runTifRegStatic, getFrames, runBackSub, Tiff.

NETWORK_FILE = 0;

% If file is not on C:\ or D:\ it's likely on a portable or network drive
if ~strcmpi(obj.Folder(1), 'c') && ~strcmpi(obj.Folder(1), 'd')
    
    NETWORK_FILE = 1;
    
    % .. so copy tiff file to a new local folder
    local_root = 'C:\tempMatlab\';
    local_dir = [local_root obj.File '\'];
    network_filepath = [obj.Folder obj.TiffFile];
    
    % If the temp folder already exists, remove it
    if exist(local_root, 'dir')
        
        % Matlab can't remove the folder if it's the current
        % working directory, so move one folder up
        if strcmp(local_dir,cd)
            cd('\')
        end
        
        % Remove non-empty folder
        rmdir(local_root,'s')
        
    end
    
    % Make local temp folder and subfolder for this file
    mkdir(local_root)
    mkdir(local_dir)
    
    % Copy into the subfolder the object's Tiff file
    disp(['Copying files to ' local_root ''])
    copyfile(network_filepath, local_dir)
    
else
    % If the Tiff file is on a local drive, work on it directly
    local_dir = [obj.Folder];
end


% Run registration:
ca_RegSeries_v4_Apply_regDispXY(obj.TiffFile, local_dir, regDispXYfullpath);

% Copy reg_tif back to original folder on network drive
if NETWORK_FILE
    disp('Copying files back to network drive..')
    copyfile([local_dir obj.TifRegFile], obj.Folder)
    rmdir(local_root, 's')
end