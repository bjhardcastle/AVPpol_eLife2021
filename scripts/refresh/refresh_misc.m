% Refresh various stored mat files containing derived values used for plotting
thisfilepath = mfilename('fullpath');
thisfolderpath = fullfile(fileparts(thisfilepath));
misc_folderpath = fullfile(thisfolderpath,'misc');
if exist(misc_folderpath,'dir')
    refresh_scripts = cellstr(ls(fullfile(misc_folderpath,'*.m')));
    for ridx = 1:length(refresh_scripts)
        run(refresh_scripts{ridx})
    end
else
    warning('\scripts\refresh\mat\misc not found - check relative file paths')
end