function getRecProps(obj)
% For reference, get the parent folder names:
parts = strsplit(obj.Folder,filesep);
obj.Line = parts{end-1};
obj.Cell = parts{end-2};
areaStr = parts{end-3};
obj.Area = areaStr(3:end);
end