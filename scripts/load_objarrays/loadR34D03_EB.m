pathsAVP
savename = 'R34D03_EB';
if ~exist('x') || ~ ( strcmp(x(1).Line,'R34D03-gal4') && strcmp(x(1).Area,'EB') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray



% single tuning: whole EB
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190521'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1810' )));
% responses concentrated: upper left
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'181008'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1312' )));
selectObj(3) = find(~cellfun('isempty',strfind({x.DateStr},'190522'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1132' )));
selectObj(4) = find(~cellfun('isempty',strfind({x.DateStr},'190522'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1203' )));
% responses concentrated: upper right
selectObj(5) = find(~cellfun('isempty',strfind({x.DateStr},'191231'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1639' )));
% single tuning: upper right
selectObj(6) = find(~cellfun('isempty',strfind({x.DateStr},'190530'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1946' )));
% full rec which zoomed upper right above comes from 
selectObj(7) = find(~cellfun('isempty',strfind({x.DateStr},'190530'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1942' )));

selectObj(8) = find(~cellfun('isempty',strfind({x.DateStr},'190522'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1139' )));

