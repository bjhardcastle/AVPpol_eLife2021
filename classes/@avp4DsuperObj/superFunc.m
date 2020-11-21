function superVal = superFunc(obj,funcString,propString)
assert(nargin == 3,'Enter string for function and property')
assert(ischar(funcString),'Enter character strings only')
assert(ischar(propString),'Enter character strings only')
assert(isprop(obj(1).Layers(1),propString),'String does not match a valid property of the object')

propStore = [];
for oidx = 1:length(obj)
for lidx = 1:length(obj(oidx).Layers)
    if isempty(obj(oidx).Layers(lidx).(propString))
        disp([propString ' is empty in obj.(' num2str(oidx) ').Layers(' num2str(lidx) ')'])
    else
        propStore = [propStore;obj(oidx).Layers(lidx).(propString)(:)];
    end
end
end

eval(['superVal = ' funcString '(propStore);'])