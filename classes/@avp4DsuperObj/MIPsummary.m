function MIPsummary(obj)
if ~isempty(obj.MIP)
    if isempty(obj.MIP.fullAngMaskImg)
        getPolROIs(obj.MIP)
    end
    plotPolData(obj.MIP)
else
    disp('No MIP exists')
end
end