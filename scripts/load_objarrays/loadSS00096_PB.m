pathsAVP
savename = 'SS00096_PB';
if ~exist('x','var') || ~( strcmp(x(1).Line,'SS00096-gal4') && strcmp(x(1).Area,'PB') )
	loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


x(find(~cellfun('isempty',strfind({x.DateStr},'191223'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1320' ))))= [];

% clear some unresponsive recordings or recordings of only one side of the PB
excl = [];
excl = [excl, find([x.Fly]==1911211)];
excl = [excl, find([x.Fly]==2001025)];
excl = [excl, find([x.Fly]==2001021)];
% Two recordings with short angle presentations
excl = [excl, find(~cellfun('isempty',strfind({x.DateStr},'191212'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1334' )))];
excl = [excl, find(~cellfun('isempty',strfind({x.DateStr},'191212'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1337' )))];
for xIdx = 1:length(x)
 if x(xIdx).pSet(4).trialTestPauseLength ~= 4 ||...
         x(xIdx).pSet(4).trialRandomizeOrder ~= 0 ||...
         isfield(x(xIdx).pSet(4),'polOffBetweenTrials') && x(xIdx).pSet(4).polOffBetweenTrials ~=0
    excl = [excl,xIdx]; 
 end
end
x(excl) = [];

incl = 1:length(x);
stepdir = nan(length(incl),1);
for xIdx = 1:length(x(incl))
   if x(incl(xIdx)).MIP.pSet(4).StepDIR == 5
       stepdir(xIdx) = 0;
   elseif x(incl(xIdx)).MIP.pSet(4).StepDIR == 1
       stepdir(xIdx) = 1;
   end
end
inclFwd = incl(stepdir==1);
inclRev = incl(stepdir==0);

selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'191231'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1428' )));
selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'191231'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1420' )));
selectObj(3) = find(~cellfun('isempty',strfind({x.DateStr},'191114'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1810' )));
selectObj(4) = find(~cellfun('isempty',strfind({x.DateStr},'191114'))  &  ~cellfun('isempty',strfind({x.TimeStr},'2000' )));
selectObj(5) = find(~cellfun('isempty',strfind({x.DateStr},'191115'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1926' )));
selectObj(6) = find(~cellfun('isempty',strfind({x.DateStr},'191224'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1726' )));
selectObj(7) = find(~cellfun('isempty',strfind({x.DateStr},'191224'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1642' )));
selectObj(8) = find(~cellfun('isempty',strfind({x.DateStr},'191224'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1636' )));
selectObj(9) = find(~cellfun('isempty',strfind({x.DateStr},'200102'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1304' )));
selectObj(10) = find(~cellfun('isempty',strfind({x.DateStr},'200102'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1256' )));
selectObj(11) = find(~cellfun('isempty',strfind({x.DateStr},'200102'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1730' )));

