pathsAVP
savename = 'R88A06_Bu_ant';
if ~exist('x') || ~( strcmp(x(1).Line,'R88A06-gal4') && strcmp(x(1).Area,'Bu') && strcmp(x(1).Layers(1).altMaskStr,'AnteriorAlt') ) 
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% No anterior bulb in this recording
x(find(~cellfun('isempty',strfind({x.DateStr},'190624')) & ~cellfun('isempty',strfind({x.TimeStr},'1535'))  )) = [];
% Massive drift in this recording
x(find(~cellfun('isempty',strfind({x.DateStr},'190625')) & ~cellfun('isempty',strfind({x.TimeStr},'2126'))  )) = [];
% Bright explosion in this recording
x(find(~cellfun('isempty',strfind({x.DateStr},'190625')) & ~cellfun('isempty',strfind({x.TimeStr},'2159'))  )) = [];

% % Right hand side example
% selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190628'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1437' )));
% selectLayer(1) = 4;
% % % Right hand side, different fly
% selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'190628'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1425' )));
% selectLayer(2) = 5;
% % Right hand side example
% selectObj(3) = find(~cellfun('isempty',strfind({x.DateStr},'190628'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1522' )));
% selectLayer(3) = 6;

% same as in superior bulb version:
% % Right hand side, imaging layer in superior bulb
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190628'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1511' )));
selectLayer(1) = 4;
% Same fly, deeper (more anterior) layer in anterior bulb
selectObj(2) = selectObj(1);
selectLayer(2) = 8;

selectObj(3) = find(~cellfun('isempty',strfind({x.DateStr},'190625'))  &  ~cellfun('isempty',strfind({x.TimeStr},'2040' )));
selectLayer(3) = 8;
% same fly as above, left side 
selectObj(4) = find(~cellfun('isempty',strfind({x.DateStr},'190625'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1959' )));
selectLayer(4) = 8;