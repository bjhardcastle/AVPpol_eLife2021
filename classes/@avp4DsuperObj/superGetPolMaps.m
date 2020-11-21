function superGetPolMaps(objarray)
for oidx = 1:length(objarray)
    for lidx = 1:length(objarray(oidx).Layers)
        getPolMaps(objarray(oidx).Layers(lidx));
    end
    getPolMaps(objarray(oidx).MIP);
end