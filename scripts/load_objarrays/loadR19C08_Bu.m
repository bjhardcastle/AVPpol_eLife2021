pathsAVP
savename = 'R19C08_Bu';
if ~exist('x') || ~ ( strcmp(x(1).Line,'R19C08-gal4') && strcmp(x(1).Area,'Bu') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray



% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'191108'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1425' )));
selectLayer(1) = 4;
% % Right hand side, different fly
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'191108'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1512' )));
selectLayer(2) = 3;
