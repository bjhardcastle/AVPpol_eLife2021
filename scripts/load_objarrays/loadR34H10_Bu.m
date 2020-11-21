pathsAVP
savename = 'R34H10_Bu';
if ~exist('x') || ~ ( strcmp(x(1).Line,'R34H10-gal4') && strcmp(x(1).Area,'Bu') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% Bad drift in this recording:
x(find(~cellfun('isempty',strfind({x.DateStr},'181009'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1422' )))) = [];

% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190622'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1432' )));
selectLayer(1) = 6;
% % Right hand side, different fly
% selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'190622'))  &  ~cellfun('isempty',strfind({x.TimeStr},'0944' )));
% Left hand side, same fly
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'190622'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1414' )));
selectLayer(2) = 8;
