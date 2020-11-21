function getSubObjArrays(obj)

%% Find individual layers tiffs
activeLayersFolder = fullfile(obj.Folder,'layers',obj.activeChanStr);
tifList1 = ls(fullfile(activeLayersFolder,['*_' obj.activeChanStr '._reg.tif']));
if isempty(tifList1)
    tifList1 = ls(fullfile(activeLayersFolder,['*_' obj.activeChanStr '.tiff']));
    assert(~isempty(tifList1),['No tiff files in ''layers\' obj.activeChanStr '\'': ' obj.Link])
    disp(['Tiff files in ''layers\' obj.activeChanStr '\'' have not been registered: ' obj.Link])
end
assert(size(tifList1,1) == obj.numZPlanes, ['Number of .tiff files does not equal obj.numZPlanes in ''layers\' obj.activeChanStr '\'': ' obj.Link])

% Go through all or selection of layers, construct objects and store in object array
if isempty(obj.ZSelect)
    ZPlanes = 1:obj.numZPlanes;
else
    ZPlanes = obj.ZSelect;
end
fCount = 0;
for fIdx = ZPlanes
    fCount = fCount+1;
    subObjarray(fCount) = avp4DmaxObj(fullfile(activeLayersFolder,tifList1(fIdx,:)));
end
obj.Layers = subObjarray;


%% Find MIP tiff
if obj.numZPlanes > 1
    
    activeMIPFolder = fullfile(obj.Folder,'projection_max',obj.activeChanStr);
    tifList2 = ls(fullfile(activeMIPFolder,['*_' obj.activeChanStr '._reg.tif']));
    if isempty(tifList2)
        tifList2 = ls(fullfile(activeMIPFolder,['*_' obj.activeChanStr '.tiff']));
        assert(~isempty(tifList2),['No tiff files in ''projection_max\' obj.activeChanStr '\'': ' obj.Link])
        disp(['Tiff file in ''projection_max\' obj.activeChanStr '\'' have not been registered: ' obj.Link])
    end
    assert(size(tifList2,1) == 1, ['More than one .tiff file in ''projection_max\' obj.activeChanStr '\'': ' obj.Link])

    % Construct object and store in object array
    mipObj = avp4DmaxObj(fullfile(activeMIPFolder,tifList2));
    obj.MIP = mipObj;
    
end

%% Find individual layers static tiffs
if  obj.Channels > 1
staticLayersFolder = fullfile(obj.Folder,'layers',obj.staticChanStr);
tifList3 = ls(fullfile(staticLayersFolder,['*_' obj.staticChanStr '._reg.tif']));
if isempty(tifList3)
    tifList3 = ls(fullfile(staticLayersFolder,['*_' obj.staticChanStr '.tiff']));
    assert(~isempty(tifList3),['No tiff files in ''layers\' obj.staticChanStr '\'': ' obj.Link])
    disp(['Tiff files in ''layers\' obj.staticChanStr '\'' have not been registered: ' obj.Link])
end
assert(size(tifList3,1) == obj.numZPlanes, ['Number of .tiff files does not equal obj.numZPlanes in ''layers\' obj.activeChanStr '\'': ' obj.Link])

% Go through all or selection of layers, construct objects and store in object array
if isempty(obj.ZSelect)
    ZPlanes = 1:obj.numZPlanes;
else
    ZPlanes = obj.ZSelect;
end
fCount = 0;
for fIdx = ZPlanes
    fCount = fCount+1;
    staticsubObjarray(fCount) = avp4DmaxObj(fullfile(staticLayersFolder,tifList3(fIdx,:)));
end
obj.staticLayers = staticsubObjarray;
end

%% Find static MIP tiff
if obj.numZPlanes > 1 && obj.Channels > 1
    
    staticMIPFolder = fullfile(obj.Folder,'projection_max',obj.staticChanStr);
    tifList4 = ls(fullfile(staticMIPFolder,['*_' obj.staticChanStr '._reg.tif']));
    if isempty(tifList4)
        tifList4 = ls(fullfile(staticMIPFolder,['*_' obj.staticChanStr '.tiff']));
        assert(~isempty(tifList4),['No tiff files in ''projection_max\' obj.staticChanStr '\'': ' obj.Link])
        disp(['Tiff file in ''projection_max\' obj.staticChanStr '\'' have not been registered: ' obj.Link])
    end
    assert(size(tifList4,1) == 1, ['More than one .tiff file in ''projection_max\' obj.staticChanStr '\'': ' obj.Link])

    % Construct object and store in object array
    staticmipObj = avp4DmaxObj(fullfile(staticMIPFolder,tifList4));
    obj.staticMIP = staticmipObj;
    
end

end