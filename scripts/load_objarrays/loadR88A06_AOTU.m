pathsAVP
savename = 'R88A06_AOTU';
if ~exist('x') || ~( strcmp(x(1).Line,'R88A06-gal4') && strcmp(x(1).Area,'AOTu') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% one of two recordings in one fly which has a lot of dendrites in the
% MIP image making it difficult to interpret.
x(find(~cellfun('isempty',strfind({x.DateStr},'190628')) & ~cellfun('isempty',strfind({x.TimeStr},'1053')))  ) = [];

% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190625'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1544' )));
selectLayer(1) = 4;

% Right hand side example, ventral polarotopy
selectObj(2) = nan;
selectLayer(2) = nan;

% Right hand side example
selectObj(3) = find(~cellfun('isempty',strfind({x.DateStr},'190625'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1432' )));
selectLayer(3) = 4;

% Right hand side example
selectObj(4) = find(~cellfun('isempty',strfind({x.DateStr},'190624'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1526' )));
selectLayer(4) = 4;