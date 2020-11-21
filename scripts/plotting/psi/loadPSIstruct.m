function [PSI, objNames] = loadPSIstruct()

% load PSI struct
pathsAVP
% saveFileDir = 'Y:/ben/avp_pol/18_data/plotting/pol_sel_funcs/PSIstruct_exclSinglePlanes_std.mat';
load(psi_mat_path)

getAVPplotParams

objNames = fieldnames(PSI);

for fIdx = 1:length(objNames)
   PSI.(objNames{fIdx}) = orderfields(PSI.(objNames{fIdx}));
end

end
