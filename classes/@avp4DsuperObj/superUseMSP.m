function superUseMSP(objarray,value)
if nargin < 2 || isempty(value)
    disp('Setting UseMSP = 1 (use input argument to set to 0)')
    value = 1;
end
for oidx = 1:length(objarray)
   if ~isempty(objarray(oidx).MIP) && (length(objarray(oidx).Layers) > 1) && exist([objarray(oidx).MIP.Folder,objarray(oidx).MIP.DateStr,objarray(oidx).MIP.TimeStr,'_MSP.mat'])
       objarray(oidx).MIP.UseMSP = value;
   elseif objarray(oidx).numZPlanes == 1
       objarray(oidx).MIP = objarray(oidx).Layers(1);
        objarray(oidx).MIP.UseMSP = 0;
   else 
       objarray(oidx).MIP.UseMSP = 0;
   end
   for lidx = 1:objarray(oidx).numZPlanes
       objarray(oidx).Layers(lidx).UseMSP = 0;
   end
end
