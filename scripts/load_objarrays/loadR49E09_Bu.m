pathsAVP
savename = 'R49E09_Bu';
if ~exist('x') || ~ ( strcmp(x(1).Line,'R49E09-gal4') && strcmp(x(1).Area,'Bu') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'191025'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1645' )));
selectLayer(1) = 7;
% % Right hand side, different fly
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'191025'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1710' )));
selectLayer(2) = 8;
