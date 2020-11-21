function getLogData(obj)
% Read the log file for a particular tiff file and extract date and time of
% experiment.

% Log filenames are slightly different to the tiff files they're associated
% with. Fix the filename first:

% Find "T0" followed by additional "0"s in Tif filename
% logexpr = 'T0*0+_'; 
% logFile = [regexprep(obj.File, logexpr, 'T0_') '.log'];
filenamePrefix = obj.Name(8:end);
logFiles = cellstr( ls( fullfile(obj.Folder,[filenamePrefix '*.log']) ) );
logFile = logFiles{end};
% if exist([obj.Folder logFile]) == 2
% Read from .log text file
log_entry = textread(fullfile(obj.Folder, logFile), '%s', 'delimiter', '\n', 'whitespace', '');

% Extract info from 2nd line (Capture time)
% Use some markers as guides for extracting date + time
whitespaces = regexp(log_entry{2}, '\s');
colons = regexp(log_entry{2}, '\:');
fwdslashes = regexp(log_entry{2}, '/');

% Store capture time (used to find corresp. daq/param files)
hrstr = log_entry{2}(whitespaces(3) + 1 : colons(2) -1 );
minstr = log_entry{2}(colons(2) +1 : colons(3) -1 );

obj.TimeStr = strcat([hrstr minstr]);

% Store capture date
yrstr = log_entry{2}(fwdslashes(2) + 3 : fwdslashes(2) + 4);
daystr = log_entry{2}(fwdslashes(1) + 1 : fwdslashes(2) - 1);
if length(daystr)<2,daystr=['0' daystr];end
mthstr = log_entry{2}(whitespaces(2) + 1  : fwdslashes(1) - 1);
if length(mthstr)<2,mthstr=['0' mthstr];end

obj.DateStr = strcat([yrstr mthstr daystr]);

% Store number of ZPlanes & Channels 
obj.numZPlanes = str2double( log_entry{3}(11:end) );
obj.Channels = str2double( log_entry{5}(11:end) );
obj.micronsPerPixel = str2double( log_entry{6}(20:end) );
obj.micronsPerZStep = str2double( log_entry{7}(22:end) );

obj.LogFile = logFile;

%     error(['Could not find ' logFile ' in folder <a href="matlab:winopen(''' obj.Folder ''')">(Open in Explorer)</a> Log file should have the same name as the tiff file it belongs to.'''])
end

