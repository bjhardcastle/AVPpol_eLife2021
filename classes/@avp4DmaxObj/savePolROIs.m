function savePolROIs(objarray)
%SAVEPOLROIS Save ROI structures to disk in their own .mat file
% savePolROIs(objarray)

for oidx = 1:length(objarray)
    if isprop(objarray(oidx),'altMaskStr') && ~isempty(objarray(oidx).altMaskStr)
        maskStr = ['_' objarray(oidx).altMaskStr '_'];
    else
        maskStr = '';
    end
    filename = [objarray(oidx).DateStr objarray(oidx).TimeStr maskStr 'PolROIs.mat'];
   folder = [objarray(oidx).Folder];
   polROI = objarray(oidx).polROI;
   save( fullfile(folder,filename), 'polROI' );
end