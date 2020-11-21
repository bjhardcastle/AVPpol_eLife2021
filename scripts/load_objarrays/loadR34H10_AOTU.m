pathsAVP
savename = 'R34H10_AOTU';
if ~exist('x') || ~ ( strcmp(x(1).Line,'R34H10-gal4') && strcmp(x(1).Area,'AOTu') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190611'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1627' )));
selectLayer(1) = 5;
% Different fly
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'190612'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1423' )));
selectLayer(2) = 6;

