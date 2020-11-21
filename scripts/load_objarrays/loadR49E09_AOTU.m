pathsAVP
savename = 'R49E09_AOTU';
if ~exist('x') || ~ ( strcmp(x(1).Line,'R49E09-gal4') && strcmp(x(1).Area,'AOTu') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190625'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1726' )));
selectLayer(1) = 5;

