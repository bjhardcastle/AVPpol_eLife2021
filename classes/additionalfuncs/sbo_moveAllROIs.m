function sbo_moveAllROIs(fig, callback, live_info)
%NIA_SELECTBACKSUB Open GUI to select parameters for background subtraction
%   Provides the user with graphical controls that can be used to select
%   values for two parameters. These are returned for use in another function
%   for subtracting the background intensity in an object's image frames.
%   If a callback function is provided then it is assumed that this
%   function is managed, as such it will hide when it
%   receives a close request rather delete itself. The managing figure must
%   delete it when appropriate.
%
%   This function accepts the following arguments:
%
%       fig - Figure to use for axes.  If empty a new figure
%           will be created
%
%       ch_ids - A [1xN] array containing the channel identifier
%           for each channel in hist_info
%
%       hist_info - A structure containing information about
%           the image histogram, more information can be found
%           below.
%
%       callback - An optional function handle that is invoked
%           and passed a colormap specification whenever the
%           current colormap is altered.  The colormap
%           specification is described below.
%
%   The argument hist_info specifies the histogram information
%   that informs selection of the colormap.  It must be a structure
%   array with an element for each channel and the fields:
%
%       min - Minimum intensity value
%       max - Maximum intensity value
%       dist - A 2xN array specifying the histogram counts. The
%           first row contains the count location, and the 
%           second row contains the count value. This
%           field may also be left empty, in which case the
%           histogram plot is left empty
%
%   The colormap specification passed to the callback function
%   is a structure array with the fields:
%
%       min - A scalar value indicating minimum value for colormap
%       max - A scalar value indicating maximum value for colormap
%       colormap - A Nx3 matrix for use a colormap
%       nan_color - A 1x3 matrix used for NaN
%

default_pixel_step = 1;

% Fill in optional arguments
if nargin < 2
    callback = [];
end

% Check input arguments
if ~ishghandle(fig)
    error 'The argument ''fig'' must be a figure handle';
end

if ~isempty(callback)
    if ~isa(callback, 'function_handle')
        error 'The argument ''callback'' has an invalid type';
    end
end

if nargin < 3
    live_info = [];
else
    fnames = {'creator', 'uuids'};
    
    if ~isstruct(live_info) || length(live_info) ~= 1
        error 'The argument ''live_info'' has an invalid type';
    end
    
    [live_info_ok, live_info_msg] = nia_sbo_hasValidFieldNames(...
        live_info, fnames, fnames);
   	if ~live_info_ok
        error(live_info_msg, 'live_info');
    end
    
    if ~ishghandle(live_info.creator)
        error 'The argument ''live_info.creator'' has an invalid type';
    end
    
    if ~iscell(live_info.uuids) || ~isvector(live_info.uuids)
        error 'The argument ''live_info.uuids'' has an invalid type';
    end
    
    for idx=1:length(live_info.uuids)
        if ~nia_sbo_isString(live_info.uuids{idx})
            error 'The argument ''live_info.uuids'' has an invalid entry';
        end
    end
end

% Create the figure if necessary
if isempty(fig)
    fig = figure;
end

% Get figure dimensions to guide button placement
figpos = get(fig, 'Position');
fw = figpos(3); % Figure width
fh = figpos(4); % Figure height
bsize = round(fw/5); % Button size 
sh = round(fw/8); % Slider height
sw = round(2*fw/3); %Slider width

% Create push buttons
 right_btn_image = imread('nia_play_icon.png');
 right_ctrl = uicontrol(...
     'CData', right_btn_image, ...
     'Style', 'pushbutton', ...
     'TooltipString', 'Move all ROIs RIGHT', ...
     'HandleVisibility', 'off', ...
     'Interruptible', 'off', ...
     'BusyAction', 'cancel', ...
     'Position', [round((fw+bsize)/2) round((fh+bsize-sh)/2) bsize bsize]);
 
 left_btn_image = imrotate(right_btn_image,180);
 left_ctrl = uicontrol(...
     'CData', left_btn_image, ...
     'Style', 'pushbutton', ...
     'TooltipString', 'Move all ROIs LEFT', ...
     'HandleVisibility', 'off', ...
     'Interruptible', 'off', ...
     'BusyAction', 'cancel', ...
     'Position', [round((fw-bsize)/2)-bsize round((fh+bsize-sh)/2) bsize bsize]);
 
up_btn_image = imrotate(right_btn_image,90);
 up_ctrl = uicontrol(...
     'CData', up_btn_image, ...
     'Style', 'pushbutton', ...
     'TooltipString', 'Move all ROIs UP', ...
     'HandleVisibility', 'off', ...
     'Interruptible', 'off', ...
     'BusyAction', 'cancel', ...
     'Position', [round((fw-bsize)/2) round((fh+bsize-sh)/2)+bsize bsize bsize]);
 
down_btn_image = imrotate(right_btn_image,270);
 down_ctrl = uicontrol(...
     'CData', down_btn_image, ...
     'Style', 'pushbutton', ...
     'TooltipString', 'Move all ROIs DOWN', ...
     'HandleVisibility', 'off', ...
     'Interruptible', 'off', ...
     'BusyAction', 'cancel', ...
     'Position', [round((fw-bsize)/2) round((fh+bsize-sh)/2)-bsize bsize bsize]);

 % Create slider for pixel step control
slider_ctrl = uicontrol(...
    'Style','slider',...
    'Min',0.05,'Max',10, ...
    'Value',default_pixel_step,...
    'SliderStep',[0.05 0.2],...
    'Position',[round((fw-sw)/2) round(sh/4) sw sh]);

% Create a display/input for the current pixel step
value_ctrl = uicontrol(...
    'Style','edit',...
    'String', default_pixel_step ,...
    'Position',[round((fw-sw)/2) round(sh*1.5) 40 20]);
% Create a label next to it
value_text = uicontrol(...
    'Style','text',...
    'String', 'pixels per step' ,...
    'Position',[round((fw-sw)/2)+45 round(sh*1.5) round(fw/2) 15]);

% Push live info and callbacks to figure data
udata.live_info = live_info;
udata.callback = callback;
set(fig, 'UserData', udata);

% Add callbacks for buttons
set(left_ctrl, 'Callback', @(src, evt) leftCallback(fig, src, evt));
set(right_ctrl, 'Callback', @(src, evt) rightCallback(fig, src, evt));
set(up_ctrl, 'Callback', @(src, evt) upCallback(fig, src, evt));
set(down_ctrl, 'Callback', @(src, evt) downCallback(fig, src, evt));
set(slider_ctrl, 'Callback', @(src, evt) sliderChangeCallback(fig, src, evt));
set(value_ctrl, 'Callback', @(src, evt) valueChangeCallback(fig, src, evt));
set(fig, ...
    'CloseRequestFcn', @(src, evt) figCloseRequestCallback(fig, src, evt));

% Push controls to figure data
udata.slider_ctrl = slider_ctrl;
udata.value_ctrl = value_ctrl;
udata.pixel_step = default_pixel_step;

set(fig, 'UserData', udata);

end

function figCloseRequestCallback(fig, ~, ~)
% This function is invoked when the figure receives a close request.

udata = get(fig, 'UserData');

if ~isempty(udata.callback)
    % Hide rather than close the window
    set(fig, 'Visible', 'off');
end

end

function rightCallback(fig, ~, ~)
% This function is invoked when the colormap selector is altered.

udata = get(fig, 'UserData');
pixel_step = udata.pixel_step;

if isempty(udata.live_info) 
    return;
end

if ~ishghandle(udata.live_info.creator)
    return;
end

creator_udata = get(udata.live_info.creator, 'UserData');

if ~isempty(creator_udata.roi_list)
    
    for idx = 1:length(creator_udata.roi_list)
        
        roi_pos = getPosition(creator_udata.roi_list(idx).handle);
        
        if isequal( size(roi_pos), [1 4]) % Ellipse
            
            % Modify ROI position by current step size
            roi_pos(1) = roi_pos(1) + pixel_step;
        
        else % Polygon
            % Modify ROI position by current step size
            roi_pos(:,1) = roi_pos(:,1) + pixel_step;
        end
        
        % Set in parent figure
        setPosition(creator_udata.roi_list(idx).handle, roi_pos);
        udata.callback(creator_udata.roi_list(idx).handle);

    end
    
end


end

function leftCallback(fig, ~, ~)
% This function is invoked when the colormap selector is altered.

udata = get(fig, 'UserData');
pixel_step = udata.pixel_step;

if isempty(udata.live_info) 
    return;
end

if ~ishghandle(udata.live_info.creator)
    return;
end

creator_udata = get(udata.live_info.creator, 'UserData');

if ~isempty(creator_udata.roi_list)
    
    for idx = 1:length(creator_udata.roi_list)
        
        roi_pos = getPosition(creator_udata.roi_list(idx).handle);
        
        if isequal( size(roi_pos), [1 4]) % Ellipse
            
            % Modify ROI position by current step size
            roi_pos(1) = roi_pos(1) - pixel_step;
        
        else % Polygon
            % Modify ROI position by current step size
            roi_pos(:,1) = roi_pos(:,1) - pixel_step;
        end
        
        % Set in parent figure
        setPosition(creator_udata.roi_list(idx).handle, roi_pos);
        udata.callback(creator_udata.roi_list(idx).handle);
    end
    
end


end

function upCallback(fig, ~, ~)
% This function is invoked when the colormap selector is altered.

udata = get(fig, 'UserData');
pixel_step = udata.pixel_step;

if isempty(udata.live_info) 
    return;
end

if ~ishghandle(udata.live_info.creator)
    return;
end

creator_udata = get(udata.live_info.creator, 'UserData');

if ~isempty(creator_udata.roi_list)
    
    for idx = 1:length(creator_udata.roi_list)
        
        roi_pos = getPosition(creator_udata.roi_list(idx).handle);
        
        if isequal( size(roi_pos), [1 4]) % Ellipse
            
            % Modify ROI position by current step size
            roi_pos(2) = roi_pos(2) - pixel_step;
        
        else % Polygon
            % Modify ROI position by current step size
            roi_pos(:,2) = roi_pos(:,2) - pixel_step;
        end
        
        % Set in parent figure
        setPosition(creator_udata.roi_list(idx).handle, roi_pos);
        udata.callback(creator_udata.roi_list(idx).handle);

    end
    
end


end

function downCallback(fig, ~, ~)
% This function is invoked when the colormap selector is altered.

udata = get(fig, 'UserData');
pixel_step = udata.pixel_step;

if isempty(udata.live_info) 
    return;
end

if ~ishghandle(udata.live_info.creator)
    return;
end

creator_udata = get(udata.live_info.creator, 'UserData');

if ~isempty(creator_udata.roi_list)
    
    for idx = 1:length(creator_udata.roi_list)
        
        roi_pos = getPosition(creator_udata.roi_list(idx).handle);
        
        if isequal( size(roi_pos), [1 4]) % Ellipse
            
            % Modify ROI position by current step size
            roi_pos(2) = roi_pos(2) + pixel_step;
        
        else % Polygon
            % Modify ROI position by current step size
            roi_pos(:,2) = roi_pos(:,2) + pixel_step;
        end
        
        % Set in parent figure
        setPosition(creator_udata.roi_list(idx).handle, roi_pos);
        udata.callback(creator_udata.roi_list(idx).handle);

    end
    
end


end

function valueChangeCallback(fig,~,~)
udata = get(fig, 'UserData');

% Get the value entered in the textbox
str = udata.value_ctrl.String;
val = str2num(str);

if isempty(val)
    % If it isn't a number, pop up an advisory dialog
    warndlg(['Input must be numerical, in range ' num2str(udata.slider_ctrl.Min) '-' num2str(udata.slider_ctrl.Max) '']);
    % and restore current value
    udata.value_ctrl.String = udata.slider_ctrl.Value;

else
    % Otherwise, enforce the range of allowable values set for the slider
    if val > udata.slider_ctrl.Max
        pixel_step = udata.slider_ctrl.Max;
    elseif val < udata.slider_ctrl.Min
        pixel_step = udata.slider_ctrl.Min;
    else
        pixel_step = val;
    end
    
    % Update the text box and slider position
    udata.value_ctrl.String = pixel_step;
    udata.slider_ctrl.Value = pixel_step;
    set(fig, 'UserData', udata);
    
    % Then update the pixel step size for l/r/u/d buttons
    stepValueChanged(fig, pixel_step);
    
end
end

function sliderChangeCallback(fig,~,~)
udata = get(fig, 'UserData');

% Get the slider value
val = udata.slider_ctrl.Value;
% Change the numerical display next to the slider
udata.value_ctrl.String = val;

set(fig, 'UserData', udata);

% Update the pixel step size for l/r/u/d buttons
stepValueChanged(fig, val);

end

function stepValueChanged(fig,pixel_step) 
udata = get(fig, 'UserData');
udata.pixel_step = pixel_step;
set(fig, 'UserData', udata);
end
