function saveLayerMasks(objarray)
%saveLayerMasks save layer mask structures to disk
% saveLayerMasks(objarray)

for oidx = 1:length(objarray)
    if isempty(objarray(oidx).ZPlane)
        layerStr = '_';
    else
        layerStr = ['_Z' num2str(objarray(oidx).ZPlane) '_'];
    end
    if isprop(objarray(oidx),'altMaskStr') && ~isempty(objarray(oidx).altMaskStr)
        maskStr = [objarray(oidx).altMaskStr '_'];
    else
        maskStr = '';
    end
    filename = [objarray(oidx).DateStr objarray(oidx).TimeStr layerStr maskStr 'layerMask.mat'];
    folder = [objarray(oidx).Folder];
    layerMask = objarray(oidx).layerMask;   
    save( fullfile(folder,filename), 'layerMask' );
    if isempty(layerMask)
        disp(['obj(' num2str(oidx) ').layerMask is empty'])
    end
end
