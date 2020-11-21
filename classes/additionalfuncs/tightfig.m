function varargout = tightfig(hfig)
% tightfig: Alters a figure so that it has the minimum size necessary to
% enclose all axes in the figure without excess space around them.
% 
% Note that tightfig will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfig will produce an error, as these cannot easily be dealt with.
% 
% Input
%
% hfig - handle to figure, if not supplied, the current figure will be used
%   instead.
%
%

    if nargin == 0
        hfig = gcf;
    end

    % There can be an issue with tightfig when the user has been modifying
    % the contnts manually, the code below is an attempt to resolve this,
    % but it has not yet been satisfactorily fixed
%     origwindowstyle = get(hfig, 'WindowStyle');
    set(hfig, 'WindowStyle', 'normal');
    
    % 1 point is 0.3528 mm for future use

    % get all the axes handles note this will also fetch legends and
    % colorbars as well
    hax = [findall(hfig, 'type', 'axes') findall(hfig, 'type', 'polaraxes') ];
    % TODO: fix for modern matlab, colorbars and legends are no longer axes
    hcbar = findall(hfig, 'type', 'colorbar');
    hleg = findall(hfig, 'type', 'legend');
    htxt = findall(hfig,'Type','text'); %text position is in data units
    
    % get the original axes units, so we can change and reset these again
    % later
    origaxunits = get(hax, 'Units');
    
    % change the axes units to cm
    set(hax, 'Units', 'centimeters');
    
    pos = [];
    ti = [];
    
    % get various position parameters of the axes
    if numel(hax) > 1
%         fsize = cell2mat(get(hax, 'FontSize'));
        ti = cell2mat(get(hax,'TightInset'));
        pos = [pos; cell2mat(get(hax, 'Position')) ];
    else
%         fsize = get(hax, 'FontSize');
        ti = get(hax,'TightInset');
        pos = [pos; get(hax, 'Position') ];
    end
    
    if ~isempty (hcbar)
        for hidx = 1:length(hcbar)
        set(hcbar(hidx), 'Units', 'centimeters');
        
        % colorbars do not have tightinset property
        for cbind = 1:numel(hcbar(hidx))
            %         fsize = cell2mat(get(hax, 'FontSize'));
            [cbarpos, cbarti] = colorbarpos (hcbar(hidx));

            pos = [pos; cbarpos];
            ti = [ti; cbarti];
        end
        end
    end
        
    if ~isempty (hleg)
        
        set(hleg, 'Units', 'centimeters');
        
        % legends do not have tightinset property
        if numel(hleg) > 1
            %         fsize = cell2mat(get(hax, 'FontSize'));
            pos = [pos; cell2mat(get(hleg, 'Position')) ];
        else
            %         fsize = get(hax, 'FontSize');
            pos = [pos; get(hleg, 'Position') ];
        end
        ti = [ti; repmat([0,0,0,0], numel(hleg), 1); ];
    end
    
    % ensure very tiny border so outer box always appears
    ti(ti < 0.1) = 0.1;
    
    % we will check if any 3d axes are zoomed, to do this we will check if
    % they are not being viewed in any of the 2d directions
    views2d = [0,90; 0,0; 90,0];
    
    for i = 1:numel(hax)
        
        set(hax(i), 'LooseInset', ti(i,:));
%         set(hax(i), 'LooseInset', [0,0,0,0]);
        
        % get the current viewing angle of the axes
        [az,el] = view(hax(i));
        
        % determine if the axes are zoomed
        if strcmp(hax(i).Type,'polaraxes')
            iszoomed = 0;
        else
            iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
        end
        
        % test if we are viewing in 2d mode or a 3d view
        is2d = all(bsxfun(@eq, [az,el], views2d), 2);
               
        if iszoomed && ~any(is2d)
           error('TIGHTFIG:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.') 
        end
        
    end
    
    % we will move all the axes down and to the left by the amount
    % necessary to just show the bottom and leftmost axes and labels etc.
    moveleft = min(pos(:,1) - ti(:,1));
    
    movedown = min(pos(:,2) - ti(:,2));
    
    % we will also alter the height and width of the figure to just
    % encompass the topmost and rightmost axes and lables
    figwidth = max(pos(:,1) + pos(:,3) + ti(:,3) - moveleft);
    
    figheight = max(pos(:,2) + pos(:,4) + ti(:,4) - movedown);
    
    % move all the axes
    for i = 1:numel(hax)
        
        set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
        
    end
    
    for i = 1:numel(hcbar)
        
        set(hcbar(i), 'Position', [pos(i+numel(hax),1:2) - [moveleft,movedown], pos(i+numel(hax),3:4)]);
        
    end
    
    for i = 1:numel(hleg)
        
        set(hleg(i), 'Position', [pos(i+numel(hax)+numel(hcbar),1:2) - [moveleft,movedown], pos(i+numel(hax)+numel(hcbar),3:4)]);
        
    end
    
    origfigunits = get(hfig, 'Units');
    
    set(hfig, 'Units', 'centimeters');
    
    % change the size of the figure
    figpos = get(hfig, 'Position');
    
    set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);
    
    % change the size of the paper
    set(hfig, 'PaperUnits','centimeters');
    set(hfig, 'PaperSize', [figwidth, figheight]);
         set(hfig, 'PaperPositionMode', 'manual');
         set(hfig, 'PaperPosition',[0 0 figwidth figheight]);
   
    % reset to original units for axes and figure 
    if ~iscell(origaxunits)
        origaxunits = {origaxunits};
    end

    for i = 1:numel(hax)
        set(hax(i), 'Units', origaxunits{i});
    end

    set(hfig, 'Units', origfigunits);
    
%      set(hfig, 'WindowStyle', origwindowstyle);
     if nargout == 1 
         varargout{1} = hfig;
     end
end


function [pos, ti] = colorbarpos (hcbar)

    % 1 point is 0.3528 mm
    for hidx = 1:length(hcbar)
    pos = hcbar(hidx).Position;
    ti = [0,0,0,0];
    
    if ~isempty (strfind (hcbar(hidx).Location, 'outside'))

        if strcmp (hcbar(hidx).AxisLocation, 'out')
            
            tlabels = hcbar(hidx).TickLabels;
            
            fsize = hcbar(hidx).FontSize;
            
            switch hcbar(hidx).Location
                
                case 'northoutside'
                    
                    % make exta space a little more than the font size/height
                    ticklablespace_cm = 1.1 * (0.3528/10) * fsize;
                    
                    ti(4) = ti(4) + ticklablespace_cm;
                    
                case 'eastoutside'
                    
                    maxlabellen = max ( cellfun (@numel, tlabels, 'UniformOutput', true) );
            
                    % 0.62 factor is arbitrary and added because we don't
                    % know the width of every character in the label, the
                    % fsize refers to the height of the font
                    ticklablespace_cm = (0.3528/10) * fsize * maxlabellen * 0.62;

                    ti(3) = ti(3) + ticklablespace_cm;
                    
                case 'southoutside'
                    
                    % make exta space a little more than the font size/height
                    ticklablespace_cm = 1.1 * (0.3528/10) * fsize;

                    ti(2) = ti(2) + ticklablespace_cm;
                    
                case 'westoutside'
                    
                    maxlabellen = max ( cellfun (@numel, tlabels, 'UniformOutput', true) );
            
                    % 0.62 factor is arbitrary and added because we don't
                    % know the width of every character in the label, the
                    % fsize refers to the height of the font
                    ticklablespace_cm = (0.3528/10) * fsize * maxlabellen * 0.62;

                    ti(1) = ti(1) + ticklablespace_cm;
                    
            end
            
        end
    end
    end

end

%{
function hfig = tightfigadv(hfig)
% tightfigadv (Tight Figure Advanced): Alters a figure so that it has the
% minimum size necessary to enclose all axes in the figure without excess
% space around them.
%
% tightfigadv is inspired ny tightfig and contains improvements and bug
% fixes for HG2, colorbars, legends and manually positioned labels.
%
% Note that tightfigadv will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfigadv will produce an error, as these cannot easily be dealt with.
%
% hfig - handle to figure, if not supplied, the current figure will be used
% instead.

% The following code is an extension of tightfig
% Copyright (c) 2011, Richard Crozier (BSD 2-clause license)
% See: https://au.mathworks.com/matlabcentral/fileexchange/34055-tightfig

% Modifications by: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Date: 09 November 2017
% Version: 0.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    hfig = gcf;
end

% There can be an issue with tightfigadv when the user has been modifying
% the contnts manually, the code below is an attempt to resolve this,
% but it has not yet been satisfactorily fixed
set(hfig, 'WindowStyle', 'normal');

% 1 point is 0.3528 mm for future use

% get all the axes handles note this will also fetch legends and
% colorbars as well
hax = [findall(hfig, 'type', 'axes') findall(hfig, 'type', 'polaraxes') ];
% Concatenate HG2 objects into hax as graphics array and sort
hax = [hax; findall(hfig, 'type', 'legend')];

% get all the text handles so we can determine the tightinset manually
% because manual placement of labels and text removes them from
% tightinset calculations.
htx = findall(hfig, 'type', 'text');
hcb = findall(hfig, 'type', 'colorbar');
if ~isempty(hcb)
    htx = [htx; [hcb.Label]'];
end
if ~isempty(htx)
    htx(contains({htx.Visible}, 'off')) = [];
end
% Get parents of text objects because the extent is based on the parent
% position
htxPar = get(htx, 'Parent');

% get the original units, so we can change and reset these again
% later
origaxunits = get(hax, 'Units');
origtxunits = get(htx, 'Units');
origtxparunits = cellfun(@(x) x.Units, htxPar, 'un', 0);

% change the units to cm
set(hax, 'Units', 'centimeters');
set(htx, 'Units', 'centimeters');
cellfun(@set,htxPar,...
    repmat({'Units'},numel(htxPar),1),...
    repmat({'centimeters'},numel(htxPar),1));

% get various position parameters of the axes
hax_ti_ind = arrayfun(@(x) (isa(x,'matlab.graphics.axis.Axes')),...
    hax,'UniformOutput',false);
hax_ti=hax([hax_ti_ind{:}]); % Returns only Axes handles
if numel(hax) > 1
    pos = cell2mat(get(hax, 'Position'));
    ti_ = get(hax_ti,'TightInset');
    if iscell(ti_)
        ti_ = cell2mat(ti_);
    end
    ti = [ti_; ...
        zeros(size(pos)-[sum([hax_ti_ind{:}]) 0])];
else
    pos = get(hax, 'Position');
    ti = get(hax_ti,'TightInset');
end

% get the global extents of the text objects
if numel(htx) > 1 && ~isempty(htx)
    ext = cell2mat(get(htx, 'Extent'));
    extParPos = cell2mat(get([htxPar{:}],'Position'));
    ext(:,1:2) = ext(:,1:2) + extParPos(:,1:2);
elseif ~isempty(htx)
    ext = get(htx, 'Extent');
    ext(1:2) = ext(1:2) + htxPar.Position(1:2);
end

% ensure very tiny border so outer box always appears
ti(ti < 0.1) = 0.1;

% we will check if any 3d axes are zoomed, to do this we will check if
% they are not being viewed in any of the 2d directions
views2d = [0,90; 0,0; 90,0];

for i = 1:numel(hax_ti)
    set(hax(i), 'LooseInset', ti(i,:));
    % get the current viewing angle of the axes
    [az,el] = view(hax(i));
    % determine if the axes are zoomed
     if strcmp(hax(i).Type,'polaraxes')
         iszoomed = 0;
     else
         iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
     end
    % test if we are viewing in 2d mode or a 3d view
    is2d = all(bsxfun(@eq, [az,el], views2d), 2);
    
    if iszoomed && ~any(is2d)
        error('TIGHTFIGADV:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.')
    end
end

% we will move all the axes down and to the left by the amount
% necessary to just show the bottom and leftmost axes and labels etc.
if isempty(htx), ext = pos - ti; end
moveleft = min( min(pos(:,1) - ti(:,1)), min(ext(:,1)) );
movedown = min(min(pos(:,2) - ti(:,2)), min(ext(:,2)));

% we will also alter the height and width of the figure to just
% encompass the topmost and rightmost axes and lables
if isempty(htx), ext(:,1:2) = pos(:,1:2);ext(:,3:4) = pos(:,3:4) + ti(:,3:4); end
figwidth = max(max(pos(:,1) + pos(:,3) + ti(:,3)), max(sum(ext(:,[1 3]),2))) - moveleft;
figheight = max(max(pos(:,2) + pos(:,4) + ti(:,4)), max(sum(ext(:,[2 4]),2))) - movedown;

% Resets temporary changes made to colorbar pos
if numel(hax) > 1
    pos = cell2mat(get(hax, 'Position'));
else
    pos = get(hax, 'Position');
end

% move all the axes but dont move legends whose position is dependent on
% the parent axis
for i = 1:numel(hax)
    if strcmpi(hax(i).Tag,'legend') && ~strcmpi(hax(i).Location,'none')
        continue;
    end
    set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
end

% Move any colorbars that have manual positioning set because they will not
% move with their parent axis
hcbML = hcb(strcmpi(get(hcb,'Location'),'manual')); %ML: Manual Location
for i = 1:numel(hcbML)
    if numel(hcbML) > 1
        pos = cell2mat(get(hcbML, 'Position')); % Set pos to colorbar pos
    else
        pos = get(hcbML, 'Position'); % Set pos to colorbar pos
    end
    set(hcbML(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
end

origfigunits = get(hfig, 'Units');

set(hfig, 'Units', 'centimeters');

% change the size of the figure
figpos = get(hfig, 'Position');

set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);

% change the size of the paper
set(hfig, 'PaperUnits','centimeters');
set(hfig, 'PaperSize', [figwidth, figheight]);
set(hfig, 'PaperPositionMode', 'manual');
set(hfig, 'PaperPosition',[0 0 figwidth figheight]);

% reset to original units for axes and figure
if ~iscell(origaxunits)
    origaxunits = {origaxunits};
end
if ~iscell(origtxunits)
    origtxunits = {origtxunits};
end
if ~iscell(origtxparunits)
    origtxparunits = {origtxparunits};
end

drawnow; % Force drawing of objects before units are reset.
for i = 1:numel(hax)
    set(hax(i), 'Units', origaxunits{i});
end
for i = 1:numel(htx)
    set(htx(i), 'Units', origtxunits{i});
    set(htxPar{i}, 'Units', origtxparunits{i});
end

set(hfig, 'Units', origfigunits);

end

%}