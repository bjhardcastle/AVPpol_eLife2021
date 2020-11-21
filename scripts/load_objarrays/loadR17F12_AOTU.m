pathsAVP
savename = 'R17F12_AOTU';
if ~exist('x') || ~strcmp(x(1).Line,'R17F12-gal4')
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190622'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1222' )));
selectLayer(1) = 5;
% Right hand side, different fly
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'190622'))  &  ~cellfun('isempty',strfind({x.TimeStr},'0944' )));
selectLayer(2) = 5;
