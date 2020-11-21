function addScalebar(obj,ax,micronLength,invertCol)
if nargin < 2 || isempty(ax)
    ax = gca;
end
if nargin < 3 || isempty(micronLength)
    micronLength = 5;
end
if nargin < 4 || isempty(invertCol)
    invertCol = 0;
end


reverseHold =ishold(ax);

barLength = round(micronLength/obj.micronsPerPixel);
if str2double(obj.DateStr) > 190502
    % after this date 20x objective was replaced with 40x, but never
    % updated in Slidebook software, so compensate for it here:
    barLength = 2*barLength;
end

% get axes limits
xlims = xlim;
ylims = ylim;

% place scalebar in bottom left
padding = 0.05;

xpos = xlims(2)*padding + [0 barLength];
ypos = ylims(2)*(1-padding)*[1 1];

hold(ax,'on')

% remove any existing scalebar
for cidx = length(ax.Children):-1:1
    if ~isempty(ax.Children)
        if strcmp('scalebar',ax.Children(cidx).Tag)
            delete(ax.Children(cidx));
        end
    end
end

% make new scalebar
L = line(ax,xpos,ypos,'LineWidth',2);
T = text(ax,xpos(1),ypos(1)-10,'','fontsize',7);

T.String = [ num2str(micronLength) ' ' char(956) 'm'];
if invertCol
    L.Color = 'w';
    T.Color = 'w';
else
    L.Color = 'k';
    T.Color = 'k';
end
T.FontWeight = 'bold';

L.Tag = 'scalebar';
T.Tag = 'scalebar';

if reverseHold
    hold(ax,'off');
end

end