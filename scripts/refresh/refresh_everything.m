% Refresh everything, in a particular order
pathsAVP
thisfilepath = mfilename('fullpath');
thisfolderpath = fullfile(fileparts(thisfilepath));
mat_refresh_folderpath = thisfolderpath; %fullfile(thisfolderpath,'mat');
if exist(mat_refresh_folderpath,'dir')
    
    % refresh objects first:
    % run(fullfile(mat_refresh_folderpath,'refresh_objarrays.m'))
    
    % % change 'orig_path_mode' back to 0
    
    % copy all raw data: 
    % run(fullfile(mat_refresh_folderpath,'refresh_aux_backup.m'))
    
    run(fullfile(mat_refresh_folderpath,'refresh_psi.m'))
    
    run(fullfile(mat_refresh_folderpath,'refresh_misc.m'))
else
    warning('\scripts\refresh\mat not found - check relative file paths')
end

make_all_figs

% thisfilepath = mfilename('fullpath');
% thisfolderpath = fullfile(fileparts(thisfilepath));
% mat_refresh_folderpath = fullfile(thisfolderpath,'mat');
% if exist(mat_refresh_folderpath,'dir')
%     refresh_scripts = cellstr(ls(fullfile(mat_refresh_folderpath,'*.m')));
%     for ridx = 1:length(refresh_scripts)
%         run(refresh_scripts{ridx})
%     end
% else
%     warning('\scripts\refresh\mat not found - check relative file paths')
% end
