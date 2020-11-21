function getParameterFile(obj)
daqName = obj.DaqFile;
pName = strrep(daqName,'MASTER','PARAMETERS');
pName = strrep(pName,'.daq','.mat');
pFolder = fullfile('Y:\ben\avp_pol\18_data\dump\',obj.DateStr);
pPath = fullfile(pFolder,pName);


%Check the new file exists before storing it
if exist(pPath, 'file') ~= 2
    disp(['Copy parameter file not found: <a href="matlab:winopen(''' pFolder ''')">dump folder</a> '])
else
    params = load(pPath,'parameterSet');
    obj.pSet = params.parameterSet;
end
