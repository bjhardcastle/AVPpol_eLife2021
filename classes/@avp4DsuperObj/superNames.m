function superNames(obj)

for oidx = 1:length(obj)
    idxName{oidx} = [ num2str(oidx) ' - ' obj(oidx).Name];
end
  disp(cellstr(idxName)')
end