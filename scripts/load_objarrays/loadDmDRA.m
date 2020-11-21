pathsAVP
savename = 'DmDRA';
if ~exist('x') || ~strcmp(x(1).Cell,'DmDRA')
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% Recording with streaking artefacts in Layers and MIP
x(find(~cellfun('isempty',strfind({x.DateStr},'190614')) & ~cellfun('isempty',strfind({x.Name},'fly2'))  )) = [];

% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'200103'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1900' )));
% Left hand side, same fly
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'200103'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1920' )));
% Right hand side example, no pol
selectObj(3) = find(~cellfun('isempty',strfind({x.DateStr},'200103'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1906' )));

% pol example
selectObj(4) = find(~cellfun('isempty',strfind({x.DateStr},'200104'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1401' )));
% no pol, same fly
selectObj(5) = find(~cellfun('isempty',strfind({x.DateStr},'200104'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1408' )));

% pol example
selectObj(6) = find(~cellfun('isempty',strfind({x.DateStr},'200103'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1836' )));
% no pol, same fly
selectObj(7) = find(~cellfun('isempty',strfind({x.DateStr},'200103'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1842' )));



