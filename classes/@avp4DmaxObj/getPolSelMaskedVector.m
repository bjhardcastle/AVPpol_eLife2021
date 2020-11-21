function varargout = getPolSelMaskedVector(obj,topBottomAllBoth)
if nargin > 1 && ~isempty(topBottomAllBoth)
    
    if ( ischar(topBottomAllBoth) && strcmp(topBottomAllBoth,'top') )  || ...
            ( isnumeric(topBottomAllBoth) && topBottomAllBoth == 1 )
        topBottomAllBoth = 1;
        
    elseif ( ischar(topBottomAllBoth) && strcmp(topBottomAllBoth,'bottom') )  || ...
            ( isnumeric(topBottomAllBoth) && topBottomAllBoth == 2 )
        topBottomAllBoth = 2;
        
    elseif ( ischar(topBottomAllBoth) && strcmp(topBottomAllBoth,'all') )  || ...
            ( isnumeric(topBottomAllBoth) && topBottomAllBoth == 3 )
        topBottomAllBoth = 3;
        
    elseif ( ischar(topBottomAllBoth) && strcmp(topBottomAllBoth,'both') )  || ...
            ( isnumeric(topBottomAllBoth) && topBottomAllBoth == 4 )
        topBottomAllBoth = 4;
        disp('both areas selected: returning 2 arrays [top bottom]')
        
    elseif ( ischar(topBottomAllBoth) && strcmp(topBottomAllBoth,'cell') )  || ...
            ( isnumeric(topBottomAllBoth) && topBottomAllBoth == 5 )
        topBottomAllBoth = 5;

    elseif ( ischar(topBottomAllBoth) && strcmp(topBottomAllBoth,'bkg') )  || ...
            ( isnumeric(topBottomAllBoth) && topBottomAllBoth == 6 )
        topBottomAllBoth = 6;
        
    elseif ( ischar(topBottomAllBoth) && strcmp(topBottomAllBoth,'threshold') )  || ...
            ( isnumeric(topBottomAllBoth) && topBottomAllBoth == 7 )
        topBottomAllBoth = 7;

    else
        topBottomAllBoth = 3;
        disp('no area selected: returning all polsel indices in layermask')
    end
    
else
    topBottomAllBoth = 3;
    disp('no area selected: returning all polsel indices in layermask')
end

loadLayerMasks(obj)

if ~isempty(obj.layerMask)
    if isempty(obj.ActivityFrame)
        if obj.UseMSP
            obj.UseMSP = 0;
            undoFlag = 1;
        else
            undoFlag = 0;
        end
        obj.Unattended = 1;
        getFrames(obj)
        obj.Frames = [];
        if undoFlag
            obj.UseMSP = 1;
        end
        obj.Frames = [];
    end
    
    maskedMaxFrameVector =obj.ActivityFrame(obj.layerMask.mask);
%         maskedMaxFrameVector =obj.polSelImg(obj.layerMask.mask);

    [~,frameIdx] = sort(maskedMaxFrameVector,'descend');
    
    topFraction = frameIdx(1:floor(length(frameIdx)*obj.polSelHistFraction));
    bottomFraction = frameIdx(ceil(length(frameIdx)-length(frameIdx).*obj.polSelHistFraction)+1:end);
    polSelFrameVector = obj.polSelImg(obj.layerMask.mask);
    
    switch topBottomAllBoth
        case 1
            varargout{1} = polSelFrameVector(topFraction);
        case 2
            varargout{1} = polSelFrameVector(bottomFraction);            
        case 3
            varargout{1} = polSelFrameVector;            
        case 4
            varargout{1} = polSelFrameVector(topFraction);
            varargout{2} = polSelFrameVector(bottomFraction);            
        case 5 %cell
            varargout{1} = obj.polSelImg(find(obj.layerMask.mask.*obj.cellMask));
        case 6 %bkg
            varargout{1} = obj.polSelImg(obj.bkgMask);
        case 7 %threshold
            varargout{1} = polSelFrameVector(polSelFrameVector>=obj.polSelThreshold);
    end
end