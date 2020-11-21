function savePBGlomeruliLimits(objarray)
%SAVEPOLROIS Save ROI structures to disk in their own .mat file
% savePolROIs(objarray)

for oidx = 1:length(objarray)
    if isprop(objarray(oidx),'altMaskStr') && ~isempty(objarray(oidx).altMaskStr)
        maskStr = ['_' objarray(oidx).altMaskStr '_'];
    else
        maskStr = '';
    end
    filename = [objarray(oidx).DateStr objarray(oidx).TimeStr maskStr 'PBglomLims.mat'];
   folder = [objarray(oidx).Folder];
   glomLims = objarray(oidx).PBglomLims;
   save( fullfile(folder,filename), 'glomLims' );
end