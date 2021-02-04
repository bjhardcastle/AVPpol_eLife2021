if ~exist('suffix') || isempty(suffix)
    suffix = '';
elseif ~strcmp('_',suffix(1))
    suffix = ['_' suffix];
end
if ~exist('prefix') || isempty(prefix)
    prefix = '';
elseif ~strcmp('_',prefix(end))
    prefix = [prefix '_'];
end

% v = ver;
% if str2double( v(1).Date(end-3:end) ) >= 2020
% 	useNewMatlabExport = 1;
% else
%     useNewMatlabExport = 0;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PNG
if ~exist(fullfile(printpath,'png'),'dir')
    mkdir(fullfile(printpath,'png'))
end
fullSavePath = fullfile(printpath,'png',[prefix savename suffix]);
set(gcf, 'InvertHardCopy', 'off');
% set(gcf,'color','none');
if exist('export_fig','file')
    export_fig(fullSavePath , '-dpng','-r400','-q100','-opengl',gcf);
else
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf,fullSavePath,'-dpng','-r400')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PDF
if ~exist(fullfile(printpath,'pdf'),'dir')
    mkdir(fullfile(printpath,'pdf'))
end
fullSavePath = fullfile(printpath,'pdf',[prefix savename suffix]);
set(gcf, 'InvertHardCopy', 'off');

% first call to print is bad, but changes something (?) about the figure
% that allows a good second call ( spent too long trying to understand
% 'print': it is just broken and export_fig isn't a complete solution
% since it cannot currently export transparent patches in vector formats)
print(gcf, fullSavePath , '-dpdf','-painters');
set(gcf,'Units','centimeters')
set(gcf,'PaperUnits','centimeters')
figPos = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 figPos(3:4)])
set(gcf,'color','none');
% second (good) call overwriting the first:
print(gcf, fullSavePath , '-dpdf','-painters');
% export_fig(fullSavePath , '-pdf','-painters',gcf); % no transparent patches

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SVG
%{
if ~exist(fullfile(printpath,'svg'),'dir')
    mkdir(fullfile(printpath,'svg'))
end
fullSavePath = fullfile(printpath,'svg',[prefix savename suffix]);
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'color','none');
print(gcf, fullSavePath , '-dsvg','-painters');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % After printing is done
close(gcf)
