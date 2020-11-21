function loadLayerMasks(objarray)
%loadLayerMasks Load layer mask structures from disk
% loadLayerMasks(objarray)

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
    if exist(fullfile(folder,filename),'file')
        load( fullfile(folder,filename), 'layerMask' );
        objarray(oidx).layerMask = layerMask;
    else
        % try to deal with an old filename format
        layerStr = '_';
        filename = [objarray(oidx).DateStr objarray(oidx).TimeStr layerStr maskStr 'layerMask.mat'];
        folder = [objarray(oidx).Folder];
        if exist(fullfile(folder,filename),'file')
            load( fullfile(folder,filename), 'layerMask' );
            objarray(oidx).layerMask = layerMask;
        end
    end
end