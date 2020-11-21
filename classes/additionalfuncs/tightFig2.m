function tightFig2(varargin)
% tightFig resizes figure so that its children fit tightly
%
% <SYNTAX>
%   tightFig(); % equivalent to tightFig(gcf);
%   tightFig(fig);
%   tightFig(_, margin);
%
% <DESCRIPTION>
% tightFig(fig) resizes given figure so that its children fit tightly.
% It is recommanded to use tightFig right before saving current figure
% since it manipulate children position and location.
% Operations on figure after tightFig may not behave the same as default.
% 
% tightFig(_, margin) sets custom margin. margin can be either ratio or
% pixel.
% 
%
% 
% <INPUT>
%     - fig (matlab.ui.Figure)
%          Target figure. The default value is current figure, gcf.
%     - margin (double, ratio) or (integer, pixel)
%          [horizontal, vertical] margin or a scalar for both.
%
% See also figure, gcf, gca
% Copyright 2019 Dohyun Kim / CC BY-NC
% Contact: kim92n@gmail.com
% Developed using MATLAB.ver 9.7 (R2019b) on Microsoft Windows 10 Enterprise
%%
args = varargin;
fig = [];
margin = 0.05;
if ~isempty(args) && isa(args{1}, 'matlab.ui.Figure')
    fig = args{1};
    args(1) = [];
end
if ~isempty(args)
    margin = args{1};
    margin = margin(:).';
end
if length(fig) > 1 % if figure is more than one.
    arrayfun(@(f) tightFig(f, margin), fig);
    return;
end
if isempty(fig)
    fig = gcf;
end
        
    
drawnow;
% setup new color
figcolor = get(fig, 'Color');
set(fig, 'Color', [255, 255, 254]/255); % Almost white
img = sum(frame2im(getframe(fig)), 3) ~= (255+255+254); % find background
% find lines containing any pixel other than background
n = size(img, 1) - find(any(img, 2), 1);
e = find(any(img, 1), 1, 'last');
w = find(any(img, 1), 1);
s = size(img, 1) - find(any(img, 2), 1, 'last');
% update pixels with some margin
if all(margin < 1) % if margin is given as percentage
    margin = [(e-w),(n-s)].*margin;
elseif length(margin) == 1
    margin = [margin, margin];
end
n = n + margin(2);
e = e + margin(1);
w = w - margin(1);
s = s - margin(2);
% Get current position setting
pos = zeros(length(fig.Children), 4);
for i = 1 : length(fig.Children)
    fig.Children(i).Units = 'pixels';
    
    if isfield(fig.Children(i), 'OuterPosition')
        pos(i,:) = fig.Children(i).OuterPosition;
    else
        pos(i,:) = fig.Children(i).Position;
    end
end
% Remove dependency between objects
% This prevents abnormal behavior of colorbar
for i = 1 : length(fig.Children)
    if isfield(fig.Children(i), 'Location')
        fig.Children(i).Location = 'manual';
    end
end
% Update child position
for i = 1 : length(fig.Children)
    if isfield(fig.Children(i), 'OuterPosition')
        fig.Children(i).OuterPosition(1) = pos(i,1) - w;
        fig.Children(i).OuterPosition(2) = pos(i,2) - s;
        fig.Children(i).OuterPosition(3:4) = pos(i,3:4);
    else
        fig.Children(i).Position(1) = pos(i,1) - w;
        fig.Children(i).Position(2) = pos(i,2) - s;
        fig.Children(i).Position(3:4) = pos(i,3:4);
    end
end
% Update Figure Position
fig.Position(3:4) = [e-w, n-s];
% Restore original color of figure.
set(fig, 'Color', figcolor);
end