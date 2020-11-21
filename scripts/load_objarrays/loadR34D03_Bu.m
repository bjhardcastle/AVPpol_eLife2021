pathsAVP
savename = 'R34D03_Bu';
if ~exist('x') || ~ ( strcmp(x(1).Line,'R34D03-gal4') && strcmp(x(1).Area,'Bu') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray



% % % Right hand side, different fly
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190704'))  &  ~cellfun('isempty',strfind({x.TimeStr},'0900' )));
selectLayer(1) = 4;
% Right hand side example
% selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'191108'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1908' )));
% selectLayer(2) = 7;

% selectObj(1) = 29