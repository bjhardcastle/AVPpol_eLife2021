pathsAVP
savename = 'R73C04_AOTU';
if ~exist('x') || ~( strcmp(x(1).Line,'R73C04-gal4') && strcmp(x(1).Area,'AOTu') )
    loadAVP
else, disp([savename ' objarray ''x'' already in workspace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Specific commands to run after loading objarray


% Some early experiments were run with the LED on throughout all experiments
% by mistake: make sure they're not included here:
x(find(~cellfun('isempty',strfind({x.Name},'LED ON')))) = [];


% % Right hand side example
selectObj(1) = find(~cellfun('isempty',strfind({x.DateStr},'181003'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1052' )));
selectLayer(1) = 7;
% % Right hand side, different fly
% selectObj(2) = find(~cellfun('isempty',strfind({x.DateStr},'190622'))  &  ~cellfun('isempty',strfind({x.TimeStr},'0944' )));
% 
% single fly with opposite polarotopy
selectObj(3) = find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1826' )));
 
 


% These rhs objects had negative vertical pos / tuning correlations (less
% than-0.15,p<0.05)or an organization of tunings that doesn't fit the
% common positive one
selectSetOutliers = [ ...
    find(~cellfun('isempty',strfind({x.DateStr},'181003'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1142' ))), ...
    find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1826' ))) ...
    find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1850' ))) ... 
    ]; 


% This object has an organization of tunings that resembles those of the 17F12 TuTu
% recordings
selectSetTuTuCompare = [ ...
find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1445' ))), ... % This one resembles ipsilateral R17F12 (left side)
];

% These objects had positive vertical pos / tuning correlations (over 0.15,p<0.05)
% selectSetPos = [ ...
selectSetPos = setdiff( find(cellfun('isempty',strfind({x.Name},' L'))), [selectSetOutliers,selectSetTuTuCompare] );

selectPBTobjects_L = [ ...
find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1231' ))), ...
find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1445' ))), ...
find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1850' ))), ...
];
selectPBTobjects_R = [ ...
find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1139' ))), ...
find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1420' ))), ...
find(~cellfun('isempty',strfind({x.DateStr},'190531'))  &  ~cellfun('isempty',strfind({x.TimeStr},'1826' ))), ...
];

getAVPplotParams
loadROIs([x([selectPBTobjects_L,selectPBTobjects_R]).MIP])
for oidx = [selectPBTobjects_L,selectPBTobjects_R]
   for ridx = 1:6
       x(oidx).MIP.ROI(ridx).color = fig2s3ROIcols(ridx,:);
   end
   x(oidx).MIP.UseFixedResp = 1;
end
