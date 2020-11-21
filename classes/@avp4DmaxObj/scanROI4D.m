function [Ftrace] = scanROI4D(objarray,maskIdx,scan_extent)

assert(length(objarray) == objarray(1).numZPlanes, ...
    'Objarray must be one complete stack of z-planes')

if nargin<3
    scan_extent = [];
end

for oidx = 1:length(objarray)
    objarray(oidx).Unattended = 1;
    
    trace_by_layer(:,oidx) = scanROI(objarray(oidx),objarray(oidx).ROI(maskIdx).mask,scan_extent);
    mask_size_by_layer(oidx) = length(find(objarray(oidx).ROI(maskIdx).mask));
    
    objarray(oidx).Unattended = 0;
end

Ftrace = nansum(trace_by_layer./mask_size_by_layer,2);