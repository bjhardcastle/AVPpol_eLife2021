function loadPolROIs(objarray)
%SAVEPOLROIS Load ROI structures from disk 
% loadPolROIs(objarray)

for oidx = 1:length(objarray)
    if isprop(objarray(oidx),'altMaskStr') && ~isempty(objarray(oidx).altMaskStr)
        maskStr = ['_' objarray(oidx).altMaskStr '_'];
    else
        maskStr = '';
    end
    filename = [objarray(oidx).DateStr objarray(oidx).TimeStr maskStr 'PolROIs.mat'];
   folder = [objarray(oidx).Folder];  
   if exist(fullfile(folder,filename),'file')
   load( fullfile(folder,filename), 'polROI' );
   objarray(oidx).polROI = polROI;
   end
end