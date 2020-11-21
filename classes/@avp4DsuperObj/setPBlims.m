function setPBlims(superobj)
% GUI for checking or creating masks for MIP objects and Layers

% Create the figure
fig = figure('CloseRequestFcn',@figDelete);
addExportFigToolbar(fig)
set(fig, 'MenuBar', 'none','Name',['Layer masking tool  |  ' superobj.Name],'NumberTitle','off');
set(fig, 'Color', [0.2 0.2 0.2]);
fig_pos = fig.Position;
fig.Position(3) = 2*fig_pos(4);
fig_bkg = fig.Color;

tbar = uitoolbar(fig);

activity_MIP_image = imread('sbo_activity_icon.png');
activity_MIP_image(~(activity_MIP_image(:,:,1)==activity_MIP_image(:,:,2) & activity_MIP_image(:,:,2)==activity_MIP_image(:,:,3))) ...
    = 1-activity_MIP_image(~(activity_MIP_image(:,:,1)==activity_MIP_image(:,:,2) & activity_MIP_image(:,:,2)==activity_MIP_image(:,:,3)));
activity_MIP_ctrl = uitoggletool(tbar, ...
    'CData', activity_MIP_image, ...
    'TooltipString', 'Show MIP tuning image', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Enable', 'on', ...
    'State', 'on');
set(activity_MIP_ctrl,'ClickedCallback', @activitiy_MIP_cb);

selectivity_MIP_image = flip(imread('sbo_activity_icon.png'),3);
selectivity_MIP_image(~(selectivity_MIP_image(:,:,1)==selectivity_MIP_image(:,:,2) & selectivity_MIP_image(:,:,2)==selectivity_MIP_image(:,:,3))) ...
    = 1-selectivity_MIP_image(~(selectivity_MIP_image(:,:,1)==selectivity_MIP_image(:,:,2) & selectivity_MIP_image(:,:,2)==selectivity_MIP_image(:,:,3)));
selectivity_MIP_ctrl = uitoggletool(tbar, ...
    'CData', flip(activity_MIP_image,3), ...
    'TooltipString', 'Show MIP selectivity image', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Enable', 'on', ...
    'State', 'off');
set(selectivity_MIP_ctrl,'ClickedCallback', @selectivity_MIP_cb);

poly_image = imread('nia_poly_icon.png');
poly_MIP_ctrl = uipushtool(tbar, ...
    'CData', poly_image, ...
    'TooltipString', 'Draw ROI on MIP', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'on', ...
    'Tag','MIP',...
    'BusyAction', 'cancel');
set(poly_MIP_ctrl,'ClickedCallback', @(src,evt) poly_MIP_cb(src,evt));

poly_layer_ctrl = uipushtool(tbar, ...
    'CData', poly_image, ...
    'TooltipString', 'Draw ROI on layer', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'on', ...
    'Tag','layer',...
    'Separator','on',...
    'BusyAction', 'cancel');
set(poly_layer_ctrl,'ClickedCallback', @(src,evt) poly_layer_cb(src,evt));

activity_image = imread('sbo_activity_icon.png');
activity_ctrl = uitoggletool(tbar, ...
    'CData', activity_image, ...
    'TooltipString', 'Show tuning image', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Enable', 'on', ...
    'State', 'off');
set(activity_ctrl,'ClickedCallback', @activitiy_image_cb);

selectivity_image = flip(imread('sbo_activity_icon.png'),3);
selectivity_ctrl = uitoggletool(tbar, ...
    'CData', selectivity_image, ...
    'TooltipString', 'Show selectivity image', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Enable', 'on', ...
    'State', 'off');
set(selectivity_ctrl,'ClickedCallback', @selectivity_image_cb);


up_image = imread('sbo_up_icon.png');
up_ctrl=uipushtool('parent',tbar, ...
    'CData', up_image, ...
    'TooltipString', 'Layer above', ...
    'Interruptible', 'off', ...
    'Separator','on',...
    'Enable', 'off');
set(up_ctrl,'ClickedCallback', @up_ctrl_cb);

down_image = imread('sbo_down_icon.png');
down_ctrl=uipushtool('parent',tbar, ...
    'CData', down_image, ...
    'TooltipString', 'Layer below', ...
    'Interruptible', 'off', ...
    'Enable', 'on');
set(down_ctrl,'ClickedCallback', @down_ctrl_cb);

save_image = imread('sbo_save_icon.png');
save_ctrl = uipushtool(tbar, ...
    'CData', save_image, ...
    'TooltipString', 'Save mask', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Separator','on',...
    'Enable', 'on');
set(save_ctrl,'ClickedCallback', @save_ctrl_cb);

xycorr_image = imread('sbo_XYcorr.png');

xycorrMIP_ctrl = uipushtool(tbar, ...
    'CData', flip(xycorr_image,3), ...
    'TooltipString', 'Plot MIP X/Y correlation', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Separator','on',...
    'Enable', 'on');
set(xycorrMIP_ctrl,'ClickedCallback', @xycorrMIP_ctrl_cb);
xycorr_ctrl = uipushtool(tbar, ...
    'CData', xycorr_image, ...
    'TooltipString', 'Plot X/Y correlation', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Separator','on',...
    'Enable', 'on');
set(xycorr_ctrl,'ClickedCallback', @xycorr_ctrl_cb);

thetacorr_ctrl = uipushtool(tbar, ...
    'CData', flip(xycorr_image,2), ...
    'TooltipString', 'Plot circ correlation', ...
    'HandleVisibility', 'off', ...
    'Interruptible', 'off', ...
    'BusyAction', 'cancel', ...
    'Separator','off',...
    'Enable', 'on');
set(thetacorr_ctrl,'ClickedCallback', @thetacorr_ctrl_cb);

%         poly_layer_ctrl = uipushtool(tbar, ...
%             'CData', poly_image, ...
%             'TooltipString', 'Draw ROI on layer', ...
%             'HandleVisibility', 'off', ...
%             'Interruptible', 'on', ...
%             'Tag','layer',...
%             'Separator','on',...
%             'BusyAction', 'cancel');
%         set(poly_layer_ctrl,'ClickedCallback', @(src,evt) poly_layer_cb(src,evt));

% Setup MIP context menu
MIPmenu = uicontextmenu('Parent', fig);
% Create context menu for clearing current MIP ROI
clear_MIP_menu = uimenu(MIPmenu, ...
    'Label', 'Clear MIP','Callback',@clear_MIP_cb);
% Create context menu for reloading current MIP ROI
reload_MIP_menu = uimenu(MIPmenu, ...
    'Label', 'Reload MIP','Callback',@reload_MIP_cb);
% Create context menu for applying MIP ROI to current Layer
apply_curr_menu = uimenu(MIPmenu, ...
    'Label', 'Apply MIP ROI to current Layer','Callback',@apply_MIP2layer_cb);
% Create context menu for applying MIP ROI to all Layers
apply_all_menu = uimenu(MIPmenu, ...
    'Label', 'Apply MIP ROI to ALL Layers','Callback',@apply_MIP2ALL_cb);
% Create context menu for saving current MIP ROI
save_MIP_menu = uimenu(MIPmenu, ...
    'Label', 'Save MIP only','Callback',@save_MIP_cb);


% Setup Layer context menu
layermenu = uicontextmenu('Parent', fig);
% Create context menu for clearing current Layer
clear_Layer_menu = uimenu(layermenu, ...
    'Label', 'Clear current Layer','Callback',@clear_layer_cb);
clear_ALL_menu = uimenu(layermenu, ...
    'Label', 'Clear ALL Layers and save','Callback',@clear_ALL_cb);
% Create context menu for reloading current Layer
reload_Layer_menu = uimenu(layermenu, ...
    'Label', 'Reload current Layer','Callback',@reload_layer_cb);
% Create context menu for applying MIP ROI to current Layer
apply_MIP2layer_menu = uimenu(layermenu, ...
    'Label', 'Apply MIP ROI to current Layer','Callback',@apply_MIP2layer_cb);
% Create context menu for applying current Layer ROI to MIP
apply_layer2MIP_menu = uimenu(layermenu, ...
    'Label', 'Apply current Layer ROI to MIP','Callback',@apply_layer2MIP_cb);

% Create context menu for saving current Layer
save_Layer_menu = uimenu(layermenu, ...
    'Label', 'Save current Layer only','Callback',@save_layer_cb);


% Create the main image axes
axWidth = 0.49;
axPad = 0.01;
axHeight = 0.92;

MIP_ax = axes;
set(MIP_ax, ...
    'Units','normalized',...
    'Position', [0.66*axPad axPad axWidth axHeight], ...
    'XTick', [], ...
    'YTick', [], ...
    'color', [0 0 0],...
    'Box', 'on');
layer_ax = axes;
set(layer_ax, ...
    'Units','normalized',...
    'Position', [0.5+0.33*axPad axPad axWidth axHeight], ...
    'XTick', [], ...
    'YTick', [], ...
    'color', [0 0 0],...
    'Box', 'on');

% Make cell array of layer names:

for n = 1:superobj.numZPlanes
    layerList{n} = ['Layer ' num2str(n,'%.2d')];
    layerVals{n} = n;
end

layer_ctrl = uicontrol(fig,...
    'Units','normalized',...
    'style','popupmenu',...
    'position',[0.5*axWidth+axPad axHeight+3*axPad 13*axPad 4*axPad],...
    'String',layerList,...
    'Interruptible', 'on', ...
    'Callback',@layer_ctrl_cb, ...
    'Enable','on');
layer_label = uicontrol(...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'String', 'Current: ', ...
    'position',[0.5*axWidth-14*axPad axHeight+3*axPad 15*axPad 4*axPad],...
    'FontSize', 8, ...
    'HorizontalAlignment', 'right', ...
    'Foreground',1-0.5*fig_bkg,...
    'Background', fig_bkg);


duplicate_label = uicontrol(...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'String', 'Copy mask from: ', ...
    'position',[1.3*axWidth-15*axPad axHeight+3*axPad 15*axPad 4*axPad],...
    'FontSize', 8, ...
    'HorizontalAlignment', 'right', ...
    'Foreground',1-0.5*fig_bkg,...
    'Background', fig_bkg);

duplicate_ctrl = uicontrol(fig,...
    'Units','normalized',...
    'style','popupmenu',...
    'position',[1.3*axWidth+axPad axHeight+3*axPad 13*axPad 4*axPad],...
    'String',layerList,...
    'Interruptible', 'off', ...
    'Enable','on');

apply_duplicate_ctrl = uicontrol(fig,...
    'Units','normalized',...
    'style','pushbutton',...
    'position',[1.3*axWidth+15*axPad axHeight+1*axPad 15*axPad 6*axPad],...
    'Interruptible', 'off', ...
    'string','Apply mask',...
    'Callback',@copy_Layer2Curr_cb,...
    'Enable','on');
scan_layers_ctrl = uicontrol(fig,...
    'Units','normalized',...
    'style','pushbutton',...
    'position',[5*axPad axHeight+1*axPad 15*axPad 6*axPad],...
    'Interruptible', 'off', ...
    'string','Scan all',...
    'Callback',@scanAllLayers, ...
    'Enable','on');

pol_sel_edit = uicontrol(fig,...
    'Units','normalized',...
    'style','edit',...
    'position',[1*axWidth-2*axPad axHeight+2*axPad 4*axPad 5*axPad],...
    'string','0',...
    'Callback',@pol_sel_cb, ...
    'Enable','on');

layer_label = uicontrol(...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'String', 'selectivity: ', ...
    'position',[1*axWidth-10*axPad axHeight+3*axPad 8*axPad 4*axPad],...
    'FontSize', 8, ...
    'HorizontalAlignment', 'right', ...
    'Foreground',1-0.5*fig_bkg,...
    'Background', fig_bkg);

%%
DELETE_MIP_FLAG = 0;
if isempty(superobj.MIP)
    superobj.MIP = superobj.Layers(1);
    DELETE_MIP_FLAG = 1;
    fig.Name = 'Layer masking tool         No MIP available! Using duplicate of Layer(1)';
end
% color_submenus = zeros(1, numel(cfg_color_choices));
% for submenu_idx=1:numel(cfg_color_choices)
% 	color_submenus(submenu_idx) = uimenu(color_menu,...
% 		'Label', cfg_color_choices{submenu_idx}{1}, ...
% 		'Callback', @(src, evt) roiColorCallback(fig, h, ...
%         cfg_color_choices{submenu_idx}{2}, src, evt, obj));
% end
%
loadLayerMasks(superobj.MIP)
[distMask,curveXY,closestCurvePtMap] = getCurvilinearDist(superobj.MIP,[],[],1);
loadPBGlomeruliLimits([superobj.MIP])
newLayer = 1;
currLayer = 0;
currROI = [];
MIPROI = [];
defaultLims = [];
currentLims= [];
currentPts = impoint.empty(0,15);
gLimIdx = [];
pol_sel_edit.String = num2str(superobj.Layers(1).polSelThreshold);
getDefaultLims
tuningFlag =0;
updateMIPFrame
updateAll

    function getDefaultLims(varargin)
        % setup estimated limits - can revert back to these if needed
        % There are 16 e-pg glomeruli (8 per side), but we'll give G1 (medial pair)
        % half the width of all others
        gLimits = linspace(0,1,16); % 16 boundaries to give 15 even widths.
        gLimsL = [ gLimits(1:8) 0.5]; % G1 fills remaining space up to centre
        gLimsR = 1-fliplr(gLimsL); % right side is mirror symmetry
        defaultLims = [gLimsL gLimsR];
    end
    function [ypos,xpos] = lims2curve(varargin)
        % gLimIdx = 1:16 % each glomerulus
        % Get boundaries of glomerulus, express as point on curve
        % through mask
        inclPix = find(distMask<=currentLims(gLimIdx+1+floor(gLimIdx/8)));
        [~,maxIdx] = nanmax(distMask(inclPix));
        LmaxIdx = inclPix(maxIdx);
        ypos = round(closestCurvePtMap{1}(LmaxIdx));
        xpos = round(closestCurvePtMap{2}(LmaxIdx));
    end

    function setDefaultLims(varargin)
        superobj.MIP.PBglomLims = defaultLims;
        updateAll
    end

    function getCurrentLims(varargin)
        if ~isempty(currentPts)
            % get closest point on curve to current points in gui, assign to
            % object, then find that point's distance along the curve from
            % the distance map
            gLimits = zeros(size(defaultLims));
            gLimits(end) = 1;
            for gLimIdx = 1:15
                ptpos=getPosition(currentPts(gLimIdx));
                yCurve = round(closestCurvePtMap{1}(round(ptpos(2)),round(ptpos(1))));
                xCurve = round(closestCurvePtMap{2}(round(ptpos(2)),round(ptpos(1))));
                limit = distMask(yCurve,xCurve);
                if gLimIdx < 8 % LHS
                    gLimits(gLimIdx+1) = limit;
                elseif gLimIdx == 8
                    % central: by old convention where we have left and right 
                    % limits in seperate vectors, we have the mid point at the 
                    % end of the left, and the start of the right. just continue with that here
                    gLimits(gLimIdx+1) = limit;
                    gLimits(gLimIdx+2) = limit;
                elseif gLimIdx > 8
                    gLimits(gLimIdx+2) = limit;
                end
            end
            currentLims = gLimits;
        else
            if isempty(superobj.MIP.PBglomLims)
                setDefaultLims
            end
            currentLims = superobj.MIP.PBglomLims;
        end

    end
    function pushLims2Obj(varargin)

        getCurrentLims
        superobj.MIP.PBglomLims = currentLims;

    end
    function plotCurrentPoints(varargin)
                    
        ROI_ax = MIP_ax;
        getCurrentLims
        
        delete(currentPts)
        currentPts = impoint.empty(0,15);
        updateMIPFrame

        for gLimIdx = 1:15 % each glomerulus
            
            % Get boundaries of glomerulus, express as point on curve
            % through mask
            [ypos,xpos] = lims2curve;
            currentPts(gLimIdx) =  impoint(ROI_ax,[xpos,ypos]);
            
            
            
            fcn = makeConstrainToRectFcn('impoint', ...
                get(ROI_ax,'XLim'), get(ROI_ax,'YLim'));
            setPositionConstraintFcn(currentPts(gLimIdx), fcn);
            if gLimIdx == 8
                currentPts(gLimIdx).setString(['mid'] );
                currentPts(gLimIdx).setColor([0 1 0]);
            elseif gLimIdx < 8
                currentPts(gLimIdx).setString([ num2str((8-gLimIdx+1)) '/' num2str(8-gLimIdx)] );
                currentPts(gLimIdx).setColor([1 0 0]);
            else
                currentPts(gLimIdx).setString([ num2str((gLimIdx-8)) '/' num2str(gLimIdx-7)] );
                currentPts(gLimIdx).setColor([0 0 1]);
            end
            poly_MIP_ctrl.Enable = 'off';
            
            
        end
        
    end

    function poly_layer_cb(varargin)
        ROI_ax = layer_ax;
        try
            
            currROI = impoly(ROI_ax);
            if ~isempty(currROI)
                fcn = makeConstrainToRectFcn('impoly', ...
                    get(ROI_ax,'XLim'), get(ROI_ax,'YLim'));
                setPositionConstraintFcn(currROI, fcn);
                
                currROI.setColor([1 0 0]);
                poly_layer_ctrl.Enable = 'off';
            end
            updateAll
        catch
            try delete(currROI);currROI = [];catch,end
        end
        
    end
    function poly_MIP_cb(varargin)
        ROI_ax = MIP_ax;
        try
            
            MIPROI = impoint(ROI_ax);
            if ~isempty(MIPROI)
                fcn = makeConstrainToRectFcn('impoint', ...
                    get(ROI_ax,'XLim'), get(ROI_ax,'YLim'));
                setPositionConstraintFcn(MIPROI, fcn);
                MIPROI.setString('1');
                MIPROI.setColor([1 0 0]);
                poly_MIP_ctrl.Enable = 'off';
            end
            updateAll
        catch
            try delete(MIPROI);MIPROI = [];catch,end
        end
        
    end
% Save callbacks (write to disk):
    function save_ctrl_cb(varargin)
        pushLims2Obj
        savePBGlomeruliLimits([superobj.MIP])
        updateAll
        disp('Saved.')
    end
    function xycorr_ctrl_cb(varargin)
        pushMIP2Obj
        pushCurr2Obj
        savePBGlomeruliLimits([superobj.MIP])
        updateAll
        superXY(superobj)
    end
    function xycorrMIP_ctrl_cb(varargin)
        save_MIP_cb
        pushMIP2Obj
        pushCurr2Obj
        savePBGlomeruliLimits([superobj.MIP])
        updateAll
        superMIPXY(superobj)
    end
    function thetacorr_ctrl_cb(varargin)
        pushMIP2Obj
        pushCurr2Obj
        savePBGlomeruliLimits([superobj.MIP])
        updateAll
        superRTheta(superobj)
    end
    function save_layer_cb(varargin)
        pushCurr2Obj
        savePBGlomeruliLimits(superobj.Layers(currLayer))
        updateAll
    end
    function save_MIP_cb(varargin)
        pushMIP2Obj
        savePBGlomeruliLimits(superobj.MIP)
        updateAll
    end
% Clear callbacks (remove ROI only):
    function clear_layer_cb(varargin)
        delete(currROI);
        currROI = [];
        updateUI
    end
    function clear_MIP_cb(varargin)
        delete(MIPROI);
        MIPROI = [];
        updateUI
    end
    function pol_sel_cb(varargin)
        val = str2double(pol_sel_edit.String);
        if val >=0 && val <=1
            superobj.superPolThreshold(val);
            updateLayerFrame
            updateMIPFrame
            updateAll
        else
            d = msgbox('Enter pol selectivity threshold in range [0,1]','Bad value','help');
            pol_sel_edit.String = num2str(superobj.Layers(1).polSelThreshold);
        end
    end
% Reload callbacks (load from disk and apply ROI)
    function reload_layer_cb(varargin)
        delete(currROI);
        currROI = [];
        loadPBGlomeruliLimits(superobj.Layers(currLayer))
        updateAll
    end
    function reload_MIP_cb(varargin)
        delete(MIPROI);
        MIPROI = [];
        loadPBGlomeruliLimits(superobj.MIP)
        updateAll
    end
    function apply_MIP2layer_cb(varargin)
        pushMIP2Obj
        superobj.Layers(currLayer).layerMask = superobj.MIP.layerMask;
        delete(currROI)
        currROI = [];
        updateAll
    end
    function apply_layer2MIP_cb(varargin)
        pushCurr2Obj
        superobj.MIP.layerMask = superobj.Layers(currLayer).layerMask;
        delete(MIPROI)
        MIPROI = [];
        updateAll
    end
    function apply_MIP2ALL_cb(varargin)
        pushMIP2Obj
        for Lidx = 1:superobj.numZPlanes
            superobj.Layers(Lidx).layerMask = superobj.MIP.layerMask;
        end
        delete(currROI)
        currROI = [];
        updateAll
    end
    function clear_ALL_cb(varargin)
        for Lidx = 1:superobj.numZPlanes
            superobj.Layers(Lidx).layerMask = [];
        end
        delete(currROI)
        currROI=[];
        savePBGlomeruliLimits([superobj.Layers])
        updateAll
        disp('Saved.')
        
    end
    function copy_Layer2Curr_cb(varargin)
        superobj.Layers(currLayer).layerMask = superobj.Layers(duplicate_ctrl.Value).layerMask;
        delete(currROI)
        currROI = [];
        updateAll
    end
    function layer_ctrl_cb(varargin)
        newLayer = layer_ctrl.Value;
        pushCurr2Obj
        updateAll
    end
    function up_ctrl_cb(varargin)
        if currLayer > 1
            newLayer = currLayer-1;
            pushCurr2Obj
            updateAll
        end
    end
    function down_ctrl_cb(varargin)
        if currLayer < length(superobj.Layers)
            newLayer = currLayer+1;
            pushCurr2Obj
            updateAll
        end
    end
    function activitiy_MIP_cb(varargin)
        selectivity_MIP_ctrl.State = 'off';
        pushMIP2Obj
        updateMIPFrame
        updateAll
        if strcmp(activity_MIP_ctrl.State,'on') && ~isempty(MIPROI)
            MIPROI.setColor([0 0 0]);
        elseif ~isempty(MIPROI)
            MIPROI.setColor([1 1 1]);
        end
    end
    function selectivity_MIP_cb(varargin)
        activity_MIP_ctrl.State = 'off';
        pushMIP2Obj
        updateMIPFrame
        updateAll
        if strcmp(selectivity_MIP_ctrl.State,'on') && ~isempty(MIPROI)
            MIPROI.setColor([0 0 1]);
        elseif ~isempty(MIPROI)
            MIPROI.setColor([1 1 1]);
        end
    end
    function activitiy_image_cb(varargin)
        selectivity_ctrl.State = 'off';
        pushCurr2Obj
        updateLayerFrame
        updateAll
        if strcmp(activity_ctrl.State,'on') && ~isempty(currROI)
            currROI.setColor([0 0 0]);
        elseif ~isempty(currROI)
            currROI.setColor([1 1 1]);
        end
    end
    function selectivity_image_cb(varargin)
        activity_ctrl.State = 'off';
        pushCurr2Obj
        updateLayerFrame
        updateAll
        if strcmp(selectivity_ctrl.State,'on') && ~isempty(currROI)
            currROI.setColor([0 0 1]);
        elseif ~isempty(currROI)
            currROI.setColor([1 1 1]);
        end
    end
    function updateAll(varargin)
        % Update layer info
        updateLayer
        % Update mask info
        %         updateMask
        
        plotCurrentPoints
        pushLims2Obj
        % Update UI info
        updateUI
    end
    function updateLayer(varargin)
        if newLayer ~= currLayer
            %             activity_ctrl.State = 'off';
            %             selectivity_ctrl.State = 'off';
            currLayer = newLayer;
            updateLayerFrame
        end
        
    end
    function updateLayerFrame(varargin)
        if strcmp(activity_ctrl.State,'off') && strcmp(selectivity_ctrl.State,'off')
            
            if isempty(superobj.Layers(currLayer).ActivityFrame)
                getFrames(superobj.Layers(currLayer))
                %                             getActivityFrameMax(superobj.Layers(currLayer));
                superobj.Layers(currLayer).ActivityFrame = std(superobj.Layers(currLayer).Frames,[],3);
                layerFrame = superobj.Layers(currLayer).ActivityFrame;
                getPolMaps(superobj.Layers(currLayer))
                superobj.Layers(currLayer).Frames = [];
                superobj.Layers(currLayer).Daq = [];
            else
                layerFrame = superobj.Layers(currLayer).ActivityFrame;
            end
            
            layer_im = image(layerFrame, 'Parent', layer_ax);
            layer_im.CDataMapping = 'scaled';
            layer_ax.CLim(2) =superFunc(superobj,'mean','ActivityFrame') + 4*superFunc(superobj,'std','ActivityFrame');
            colormap(layer_ax,viridis);
            delete(currROI)
            currROI =[];
            
        elseif strcmp(activity_ctrl.State,'on') && strcmp(selectivity_ctrl.State,'off')
            
            
            layer_im = plotTuningMap(superobj.Layers(currLayer),layer_ax);
            
            delete(currROI)
            currROI =[];
        elseif strcmp(activity_ctrl.State,'off') && strcmp(selectivity_ctrl.State,'on')
            
            
            layer_im = plotPolSelMap(superobj.Layers(currLayer),layer_ax);
            
            delete(currROI)
            currROI =[];
            
        end
        drawMaskLayer
        layer_im.UIContextMenu = layermenu;
        layer_ax.Color = [0 0 0];
        layer_ax.XTick = [];
        layer_ax.YTick = [];
        
    end
    function updateMIPFrame(varargin)
        if strcmp(activity_MIP_ctrl.State,'off') && strcmp(selectivity_MIP_ctrl.State,'off')
            if isempty(superobj.MIP.ActivityFrame)
                getFrames(superobj.MIP)
                %                 getActivityFrameMax(superobj.MIP)
                superobj.MIP.ActivityFrame = max(superobj.MIP.Frames,[],3);
                getPolMaps(superobj.MIP)
                superobj.MIP.Frames = [];
                superobj.MIP.Daq = [];
            end
            
            MIPFrame = superobj.MIP.ActivityFrame;
            
            %         getFrames(superobj.MIP)
            %         MIPFrame = max(superobj.MIP.Frames,[],3);
            MIP_im = image(MIPFrame, 'Parent', MIP_ax);
            colormap(MIP_ax,viridis);
            MIP_im.CDataMapping = 'scaled';
            MIP_ax.CLim(2) = mean(MIPFrame(:)) + 4*std(MIPFrame(:));
            delete(MIPROI)
            MIPROI =[];
        elseif strcmp(activity_MIP_ctrl.State,'on') && strcmp(selectivity_MIP_ctrl.State,'off')
            
            
            MIP_im = plotTuningMap(superobj.MIP,MIP_ax);
            
            delete(MIPROI)
            MIPROI =[];
        elseif strcmp(activity_MIP_ctrl.State,'off') && strcmp(selectivity_MIP_ctrl.State,'on')
            
            
            MIP_im = plotPolSelMap(superobj.MIP,MIP_ax);
            
            delete(MIPROI)
            MIPROI =[];
            
        end
        MIP_im.UIContextMenu = MIPmenu;
        MIP_ax.Color = [0 0 0];
        MIP_ax.XTick = [];
        MIP_ax.YTick = [];
        drawMaskMIP
    end
    function drawMaskMIP(varargin)
        hold(MIP_ax,'on')
        plot(MIP_ax,[superobj.MIP.layerMask.position(:,1); superobj.MIP.layerMask.position(1,1)],[superobj.MIP.layerMask.position(:,2); superobj.MIP.layerMask.position(1,2)],'Color', ROIcolor((2)),'LineWidth',2);
    end
 function drawMaskLayer(varargin)
          hold(layer_ax,'on')
        plot(layer_ax,[superobj.MIP.layerMask.position(:,1); superobj.MIP.layerMask.position(1,1)],[superobj.MIP.layerMask.position(:,2); superobj.MIP.layerMask.position(1,2)],'Color', ROIcolor((2)),'LineWidth',2);
    end
    function updateMask(varargin)
        if ~isempty(currROI)
            if ~isempty(superobj.Layers(currLayer).layerMask)
                if isequal ( currROI.getPosition, superobj.Layers(currLayer).layerMask.position)
                    currROI.setColor([1 1 1]);
                else
                    currROI.setColor([1 0 0]);
                end
            else
                currROI.setColor([1 0 0]);
            end
        elseif ~isempty(superobj.Layers(currLayer).layerMask)
            delete(currROI);
            currROI = impoly(layer_ax,superobj.Layers(currLayer).layerMask.position);
            currROI.setColor([1 1 1]);
            fcn = makeConstrainToRectFcn('impoly', ...
                get(layer_ax,'XLim'), get(layer_ax,'YLim'));
            setPositionConstraintFcn(currROI, fcn);
            addNewPositionCallback(currROI, ...
                @currROI_changed_cb);
        end
        function currROI_changed_cb(varargin)
            if strcmp(activity_ctrl.State,'on') || strcmp(selectivity_ctrl.State,'on')
                currROI.setColor([1 1 1]);
            else
                currROI.setColor([1 0 0]);
            end
        end
        if ~isempty(MIPROI)
            if ~isempty(superobj.MIP.layerMask)
                if isequal ( MIPROI.getPosition, superobj.MIP.layerMask.position)
                    MIPROI.setColor([1 1 1]);
                else
                    MIPROI.setColor([1 0 0]);
                end
            else
                MIPROI.setColor([1 0 0]);
            end
        elseif ~isempty(superobj.MIP.layerMask)
            delete(MIPROI);
            MIPROI = impoint(MIP_ax,superobj.MIP.layerMask.position);
            MIPROI.setColor([1 1 1]);
            fcn = makeConstrainToRectFcn('impoint', ...
                get(MIP_ax,'XLim'), get(MIP_ax,'YLim'));
            setPositionConstraintFcn(MIPROI, fcn);
            addNewPositionCallback(MIPROI, ...
                @MIPROI_changed_cb);
        end
        function MIPROI_changed_cb(varargin)
            MIPROI.setColor([1 0 0]);
        end
        
    end
    function updateUI(varargin)
        layer_ctrl.Value = currLayer;
                    save_ctrl.Enable = 'on';

        if isempty(MIPROI)
            poly_MIP_ctrl.Enable = 'on';
%             save_ctrl.Enable = 'off';
        else
            poly_MIP_ctrl.Enable = 'off';
%             save_ctrl.Enable = 'on';
        end
        if isempty(currROI)
            poly_layer_ctrl.Enable = 'on';
%             save_ctrl.Enable = 'off';
        else
            poly_layer_ctrl.Enable = 'off';
%             save_ctrl.Enable = 'on';
        end
        
        
        switch currLayer
            case 1
                up_ctrl.Enable = 'off';
                down_ctrl.Enable = 'on';
            case {layerVals{2:superobj.numZPlanes-1}}
                up_ctrl.Enable = 'on';
                down_ctrl.Enable = 'on';
            case superobj.numZPlanes
                up_ctrl.Enable = 'on';
                down_ctrl.Enable = 'off';
        end
    end
    function pushMIP2Obj(varargin)
        if ~isempty(MIPROI)
            superobj.MIP.layerMask.position = getPosition(MIPROI);
            superobj.MIP.layerMask.mask = createMask(MIPROI);
        end
        delete(MIPROI)
        MIPROI = [];
        updateAll
    end
    function pushCurr2Obj(varargin)
        if ~isempty(currROI)
            superobj.Layers(currLayer).layerMask.position = getPosition(currROI);
            superobj.Layers(currLayer).layerMask.mask = createMask(currROI);
        end
        delete(currROI)
        currROI = [];
        updateAll
    end

    function scanAllLayers(varargin)
        
        h = waitbar(0,'1','Name','Getting frames...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0);
        
        for sIdx = 1:superobj.numZPlanes
            if getappdata(h,'canceling')
                break
            end
            % Report current status in the waitbar's message field
            waitbar(sIdx/superobj.numZPlanes,h,sprintf(['Layer ' num2str(sIdx) '/' num2str(superobj.numZPlanes)]))
            
            if isempty(superobj.Layers(sIdx).ActivityFrame)
                getFrames(superobj.Layers(sIdx))
                % getActivityFrameMax(superobj.Layers(sIdx))
                superobj.Layers(sIdx).ActivityFrame = std(superobj.Layers(sIdx).Frames,[],3);
                getPolMaps(superobj.Layers(sIdx))
                superobj.Layers(sIdx).Frames = [];
                superobj.Layers(sIdx).Daq = [];
            end
        end
        delete(h)       % DELETE the waitbar; don't try to CLOSE it.
        
    end
%     function addROIToFig(varargin)
%         addNewPositionCallback(roi, @roiChangedCallback);
%     end
%     function roiChangedCallback
%     end

    function figDelete(src,callbackdata)
        if DELETE_MIP_FLAG
            superobj.MIP = [];
        end
        delete(fig)
    end

    function layer_im = plotTuningMap(obj,ax)
        discretize = 1;
        if isempty(obj.avgPolImg)
            getPolMaps(obj)
            obj.Frames =[];
            obj.Daq = [];
        end
        
        avgFrame = 1-sqrt(obj.avgPolImg./max(obj.avgPolImg(:)));
        tuningFrame = obj.polPix;
        polAngImg = obj.polAngImg; % discretized version
        
        tuningFrame = polAngImg;
        tuningFrame(obj.polSelImg<obj.polSelThreshold) = nan;
        %         tuningFrame(obj.fftMagImg<0.5) = nan;
        
        combFrame = zeros(size(avgFrame)    );
        if discretize
            combFrame(~isnan(tuningFrame)) = polAngImg(~isnan(tuningFrame));
        else
            combFrame(~isnan(tuningFrame)) = obj.polTuningImg(~isnan(tuningFrame));
        end
        combFrame(isnan(tuningFrame)) = -180.*avgFrame(isnan(tuningFrame));
        
        %
        % % testing new method:
        %
        %    avgFrame = sqrt(obj.avgPolImg./max(obj.avgPolImg(:)));
        %     combFrame = zeros(size(avgFrame)    );
        % combFrame(stdfilt(polAngImg,true(3))<16) = polAngImg(stdfilt(polAngImg,true(3))<16);
        % combFrame(stdfilt(polAngImg,true(3))>=16) = -180.*avgFrame(stdfilt(polAngImg,true(3))>=16);
        %
        
        layer_im = imagesc(combFrame, 'Parent', ax);
        
        
        % colormap([gray(180);[1 1 1];(colorcet('R2','N',180, 'shift', 0.25))])
        % colormap([gray(180);[1 1 1];equalisecolourmap('RGB',flipud(hsv(180)),'CIE76',[1 1 0],0,1)])
        colormap(ax,[gray(180);[1 1 1];flipud(hsv(180))])
        
        set(ax,'CLim',[-180 180])
        
        
    end


    function layer_im = plotPolSelMap(obj,ax)
        if isempty(obj.avgPolImg)
            getPolMaps(obj)
            obj.Frames =[];
            obj.Daq = [];
        end
        fmi = obj.fftMagImg;
        %         fmi(fmi>1) = 1;
        fmi = fmi;
        psi = obj.polSelImg.*fmi;
        %         psi(obj.fftMagImg<0.1) = nan;
        %         layer_im = imagesc(psi./nanmax(psi(:)), 'Parent', ax);
        layer_im = imagesc(psi, 'Parent', ax);
        
        
        % colormap([gray(180);[1 1 1];(colorcet('R2','N',180, 'shift', 0.25))])
        % colormap([gray(180);[1 1 1];equalisecolourmap('RGB',flipud(hsv(180)),'CIE76',[1 1 0],0,1)])
        colormap(ax,(magma(256)))
        
        set(ax,'CLim',[0 0.5])
        
        
    end

end



