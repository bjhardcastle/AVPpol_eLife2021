thisfilepath = mfilename('fullpath'); 
rootpathAVP = fullfile(fileparts(thisfilepath));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orig_path_mode = 0; % for testing updated objarrays before copying to minimal data storage
original_data_path = 'Y:\ben\avp_pol\18_data\';
raw_data_path = fullfile(fileparts(rootpathAVP), '18_data');

publicationXYcorrStats = 1; % plots will appear the same either way, but setting to 1 will run more time-consuming stats (total run time to make all figs becomes hours rather than minutes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add all required subfolders to search path (temporary)

% check required plotting functions are on search path:
plot_script_path = fullfile(rootpathAVP,'scripts');
if isempty(which('pbgridplot_distsigned')) % function with a unique name which should be on path
    
    % if not, try to add subfolder
    if exist(plot_script_path,'dir')
        subfolderlist = genpath(plot_script_path);
        addpath(subfolderlist);
    else
        error('Folder required for plotting not found')
    end
end
    
% check required object classes are on search path:
if isempty(which('avp4DsuperObj')) || isempty(which('avp4DmaxObj')) || isempty(which('SlidebookObj'))
    
    % if not, try to add subfolder
    obj_classes_path = fullfile(rootpathAVP,'classes');    
    if exist(obj_classes_path,'dir')
        subfolderlist = genpath(obj_classes_path);
        addpath(subfolderlist);
    else
        error('Required object classes not found')
    end  
end

% finally, add the root folder, where 'pathsAVP' is located
addpath(rootpathAVP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures:
% plotpathAVP = fullfile(rootpathAVP, ['plots_' datestr(today,'yymmdd')] );
plotpathAVP = fullfile(rootpathAVP, 'plots');

fig1path = fullfile(plotpathAVP, 'fig1' );
fig2path = fullfile(plotpathAVP, 'fig2' );
fig3path = fullfile(plotpathAVP, 'fig3' );
fig6path = fullfile(plotpathAVP, 'fig6' );
fig7path = fullfile(plotpathAVP, 'fig7' );
fig8path = fullfile(plotpathAVP, 'fig8' );
fig9path = fullfile(plotpathAVP, 'fig9' );
fig10path = fullfile(plotpathAVP, 'fig10' );

fig1s1path = fullfile(plotpathAVP, 'fig1s1' );
fig1s2path = fullfile(plotpathAVP, 'fig1s2' );
figX1path = fullfile(plotpathAVP, 'figX1' );
fig3s1path = fullfile(plotpathAVP, 'fig3s1' );
fig6s1path = fullfile(plotpathAVP, 'fig6s1' );
fig8s1path = fullfile(plotpathAVP, 'fig8s1' );
fig10s1path = fullfile(plotpathAVP, 'fig10s1' );

fig_command_window_output = fullfile(plotpathAVP, 'commandWindow.log' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % compiled mat files, derived data, ROI masks etc.
compiled_mat_path =  fullfile(rootpathAVP, 'mat' );


% compiled mat files for plotting:
plot_mat_path = fullfile(compiled_mat_path, 'plotting'); % additional data stored for some plots

% psi values extracted from all objects for convenient plotting
psi_mat_path = fullfile(compiled_mat_path, 'psi', 'PSIstruct_exclSinglePlanes_std.mat'); 

% other resp data extracted from objects
misc_mat_path = fullfile(compiled_mat_path, 'misc' ); 

% folder structure of original 18_data folder, with manually generated
% ROIs etc. and MSIP projections, no large files. Backed up in case we ever make new ROIs or
% layers for the same data, which would overwrite these:
aux_backup_path = fullfile(compiled_mat_path, 'aux_data_backup' ); 

% TuBu reponse data for blue LED flashes in AOTu and bulb 
TuBu_flash_data_path = fullfile(plot_mat_path,'TuBu_flash_resp.mat');

%  ROI tuning/selectivity data for ROIs in:
%   TuBu_a   (anterior bulb)            lineStr = 'R34H10_Bu'
%   R4m      (anterior bulb)            lineStr = 'R34D03_Bu'
%   E-PG     (protocerebral bridge)     lineStr = 'SS00096_PB'
ROIstruct_path = @(lineStr) fullfile(plot_mat_path, ['ROIstruct_' lineStr '.mat']);

% 
snapshotStruct_HalfCycle_path = @(savename) fullfile(plot_mat_path, ['snapshotStruct_HalfCycle_' savename '.mat']);
snapshotStruct_path = @(savename) fullfile(plot_mat_path, ['snapshotStruct_' savename '.mat']);

autocorrStruct_path = fullfile(plot_mat_path, 'autocorrStruct.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
objnames = {
    %%%
    'R7R8'; ...
    'DmDRA';  ...
    %%%
    'R17F12_AOTU';  ...
    %%%
    'R56F07_AOTU';  ...
    'R73C04_AOTU';  ...
    %%%
    'R49E09_AOTU';  ...
    'R88A06_AOTU';  ...
    'R34H10_AOTU';  ...
    %%%
    'R49E09_Bu';  ...
    'R88A06_Bu';  ...
    'R88A06_Bu_ant';  ...
    'R34H10_Bu';  ...
    %%%
    'R19C08_Bu';  ...
    'R34D03_Bu';  ...
    %%%
    'R34D03_EB';  ...
    %%%
    'SS00096_PB'; ...
    };
 
% relative paths to subfolders in raw_data_path and aux_backup_path:
rel_data_paths.R7R8 =           '\1_Me\R7R8\panR7-gal4\';
rel_data_paths.DmDRA =          '\1_Me\DmDRA\VT059871,R13E04-split-gal4\';
rel_data_paths.R17F12_AOTU =    '\2_AOTU\TuTu\R17F12-gal4\';
rel_data_paths.R56F07_AOTU =    '\2_AOTU\MeTu\R56F07-gal4\';
rel_data_paths.R73C04_AOTU =    '\2_AOTU\MeTu\R73C04-gal4\';
rel_data_paths.R49E09_AOTU =    '\2_AOTU\TuBu\R49E09-gal4\';
rel_data_paths.R88A06_AOTU =    '\2_AOTU\TuBu\R88A06-gal4\';
rel_data_paths.R34H10_AOTU =    '\2_AOTU\TuBu\R34H10-gal4\';
rel_data_paths.R49E09_Bu =      '\3_Bu\TuBu\R49E09-gal4\';
rel_data_paths.R88A06_Bu =      '\3_Bu\TuBu\R88A06-gal4\';
rel_data_paths.R88A06_Bu_ant =  '\3_Bu\TuBu\R88A06-gal4\';
rel_data_paths.R34H10_Bu =      '\3_Bu\TuBu\R34H10-gal4\';
rel_data_paths.R19C08_Bu =      '\3_Bu\Ring\R19C08-gal4\';
rel_data_paths.R34D03_Bu =      '\3_Bu\Ring\R34D03-gal4\';
rel_data_paths.R34D03_EB =      '\4_EB\Ring\R34D03-gal4\';
rel_data_paths.SS00096_PB =     '\5_PB\E-PG\SS00096-gal4\';

objarray_path =  fullfile(compiled_mat_path, 'objarrays' ); % [1xN] arrays of avp4DsuperObj objects for each driver line
          
for oidx = 1:length(objnames)
    thisobj = objnames{oidx};
    
    % standalone object array paths:
    objpath.(thisobj) = fullfile(objarray_path, [thisobj,'.mat']);
    
    % original object array paths:
    objpath_orig.(thisobj) = fullfile(original_data_path, rel_data_paths.(thisobj),'superObjArray_ALL.mat');
    if strcmp(thisobj,'R88A06_Bu_ant')
        objpath_orig.(thisobj) = fullfile(original_data_path, rel_data_paths.(thisobj),'superObjArray_AnteriorAlt.mat');
    end
end
load_objarray_func_path = fullfile(rootpathAVP, 'scripts', 'load_objarrays' ); % customized load functions for each object array

