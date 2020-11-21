pathsAVP
savename = 'R7R8';
if ~exist('x') || ~strcmp(x(1).Cell,'R7R8')
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% Recording not in MeDRA
x(find(~cellfun('isempty',strfind({x.DateStr},'190517')) & ~cellfun('isempty',strfind({x.TimeStr},'1741'))  )) = [];

% Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'190514'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1425' )));
% Hi res right hand example
selectObj(4) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1235' )));

% mid-medra no polarizer
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'190517'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1648' )));
% mid-medra with polarizer (same fly)
selectObj(3) = find(~cellfun('isempty',strfind({x.DateStr},'190517'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1642' )));

% mid-medra no polarizer
selectObj(5) = find(~cellfun('isempty',strfind({x.DateStr},'190514'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1412' )));
% mid-medra with polarizer (same fly)
selectObj(6) = find(~cellfun('isempty',strfind({x.DateStr},'190514'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1425' )));

% ant-medra with polarizer
selectObj(7) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1347' )));
% mid-medra with polarizer (same fly)
selectObj(8) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1342' )));
% post-medra with polarizer (same fly)
selectObj(9) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1352' )));


% ant-medra with polarizer
selectObj(10) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1329' )));
% mid-medra with polarizer (same fly)
selectObj(11) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1318' )));
% post-medra with polarizer (same fly)
selectObj(12) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1313' )));

% ant-medra with polarizer
selectObj(13) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1211' )));
% mid-medra with polarizer (same fly)
selectObj(14) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1201' )));
% post-medra with polarizer (same fly)
selectObj(15) = find(~cellfun('isempty',strfind({x.DateStr},'190516'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1155' )));
