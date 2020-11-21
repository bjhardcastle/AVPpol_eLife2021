function varargout = getExp5ActivityFrame(obj)
%GETACTIVITYFRAME Gets the average activity frame during trials
% [image] = getActivityFrame(obj,fields,detrend)
%
% This function finds the mean Tiff frame recorded during trials, subtracts
% from it the mean of all extra-trial frames, to enhance activity
% correlated with stimulus presentation, and saves the image to the object.
% If fields is provided as an input, specific trials will be searched for,
% and only their frames will be included ( see findTrials ).
%
% Writes to object (if 'fields' is not specified):
%   obj.ActivityFrame   - can be toggled on/off in 'play' to aid drawing ROIs
%
% Additional images can be also be assigned manually:
%   obj.ActivityFrame1
%   obj.ActivityFrame2
%   etc.
% and can be displayed in the 'play' GUI.
%
% This function accepts the following inputs:
%
%   fields      - (optional) structure with trial parameters and value.
%                 eg. fields.TrialPatNum = 2
%                 If provided, the activity frame will only include the
%                 trials specified in fields.
%
%   detrend     - (optional)  Set to 1 to improve signal to noise by
%                 subtracting constant linear trend or the mean from all
%                 frames. Similar to background subtraction. Takes time.
%                 Default is 0.
%
% This function returns the following outputs:
%
%   image       - (optional) 2D image data of the activity frame for the
%                 trials requested in 'fields'
%
% See also findTrials.

if isempty(obj.TrialStartSample)
    getTrialtimes(obj)
end
if isempty(obj.DaqFile)
    return
end
if isempty(obj.TrialPatNum)
    getParameters(obj)
end
if isempty(obj.Frames)
    if ~obj.Unattended
        undoflag = 1;
    else
        undoflag = 0;
    end
    obj.Unattended = 1;
    getFrames(obj);
    if undoflag
        obj.Unattended = 0;
    end
end

angles = unique(obj.TrialPatNum(obj.TrialSeqNum == 5));
    
fields=[];
fields.TrialSeqNum = 5;
allTrialIdx = findTrials(obj,fields)';

avgTrialLength = median(diff([obj.TrialStartFrame(obj.TrialSeqNum == 5); obj.TrialEndFrame(obj.TrialSeqNum == 5)]));
for aIdx = 1:length(angles)
    
    fields=[];
    fields.TrialSeqNum = 5;
    fields.TrialPatNum = angles(aIdx);
    
    PrefTrialIdx = findTrials(obj,fields)';
    OtherTrialIdx = setdiff(allTrialIdx,PrefTrialIdx);
    
    %     fields.TrialSeqNum =4;
    %     fields.TrialPatNum = mod(angle+90,360);
    %     if fields.TrialPatNum == 0
    %          fields.TrialPatNum = 360;
    %     end
    %     AntiPrefTrialIdx = findTrials(obj,fields)';
    %     fields.TrialPatNum = mod(angle+270,360);
    %     if fields.TrialPatNum == 0
    %          fields.TrialPatNum = 360;
    %     end
    %     AntiPrefTrialIdx = sort([AntiPrefTrialIdx findTrials(obj,fields)']);
    
    if nargin < 3
        if isprop(obj,'Detrend') && obj.Detrend
            detrend = 1;
        else
            detrend = 0;
        end
    end
    
    % Make a 1-D array indicating frames:
    %   0 : from non-selected trials - will be skipped
    %   1 : from pref angle trial frames - will form 'in-trial activity image'
    %  -1 : from inter-trial frames - will form 'post-trial activity image'
     if detrend
        Frames = detrendFrames(obj);
    else
        Frames = obj.Frames;
    end
    t_on = zeros(1,size(obj.Frames,3));
    p_on = zeros(1,size(obj.Frames,3));
    max_on = [];
    mIdx = 0;
    for m = PrefTrialIdx
        t_on( obj.TrialStartFrame(m):obj.TrialEndFrame(m) ) = 1;
        p_on( obj.TrialStartFrame(m):obj.TrialEndFrame(m) ) = -1;
        
        % Add post-trial sequence
        if PrefTrialIdx < m
            t_on( obj.TrialEndFrame(m)+1 : obj.TrialStartFrame(m+1)-1 ) = 1;
            p_on( obj.TrialEndFrame(m)+1 : obj.TrialStartFrame(m+1)-1 ) = 1;
            stopFrame =  obj.TrialStartFrame(m+1)-1 ;
        else
            t_on( obj.TrialEndFrame(m)+1 : obj.TrialEndFrame(m) +avgTrialLength ) = 1;
            p_on( obj.TrialEndFrame(m)+1 : obj.TrialEndFrame(m) +avgTrialLength ) = 1;
                        stopFrame =  obj.TrialEndFrame(m) +avgTrialLength;

        end
        mIdx= mIdx+1;
        max_on(:,:,mIdx) = max(Frames(:,:,obj.TrialStartFrame(m):stopFrame),[],3);
    end
    
    for n = OtherTrialIdx
                t_on( obj.TrialStartFrame(n):obj.TrialEndFrame(n) ) = -1;
                 if OtherTrialIdx < n
            t_on( obj.TrialEndFrame(n)+1 : obj.TrialStartFrame(n+1)-1 ) = -1;
        else
            t_on( obj.TrialEndFrame(n)+1 : obj.TrialEndFrame(n) +avgTrialLength ) = -1;           
        end
    end
        
   
    
    % Get mean frame where this array is 1 (within selected trials)
    maxTrialONframe = max(Frames(:,:,t_on==1),[],3);
    % Get mean frame where the array is -1 (before LED turns on, before first pol rotate trial)
    maxTrialOFFframe = max(Frames(:,:,t_on==0),[],3);
    
    % Get mean frame where this array is 1 (within selected trials)
    avgPostONframe = mean(Frames(:,:,p_on==1),3);
    % Get mean frame where the array is -1 (before LED turns on, before first pol rotate trial)
    avgPostOFFframe = mean(Frames(:,:,p_on==-1),3);
    
    avgPostOFFframe2 = mean(Frames(:,:,p_on==1),3);
    avgPostONframe2 = mean(Frames(:,:,t_on==1),3);
    
    avgOtherTrialframe = mean(Frames(:,:,t_on==-1),3);
    maxOtherTrialframe  = max(Frames(:,:,t_on==-1),[],3);
    
    avgPrefMax = mean(max_on,3);
    % avgPostOFFframe2 = avgPostOFFframe2.*(avgPostOFFframe2>(mean(avgPostOFFframe2(:)) + std(avgPostOFFframe2(:)) ));
    % avgPostONframe2 = avgPostONframe2.*(avgPostONframe2>(mean(avgPostONframe2(:)) + std(avgPostONframe2(:))));
    
    
    %{
% Get mean frame where this array is 1 (within selected trials)
avgTrialONframe = (mean(Frames(:,:,t_on==1),3));
% Get mean frame where the array is -1 (before LED turns on, before first pol rotate trial)
avgTrialOFFframe = (mean(Frames(:,:,t_on==-1),3));


    %}
    %%% When no particular trial activity is sought, it's probably better to
    %%% take the background image as the spontaneous activity before any
    %%% stimuli run, as there can be activity between trials which distorts the
    %%% subtraction
    % if  nargin < 2 || isempty(fields)
    % avgTrialOFFframe = (mean(Frames(:,:,5:10),3));
    % end
    
    % Subtract one from the other
%     aFrame{aIdx} = (maxTrialONframe - maxTrialOFFframe);
 aFrame{aIdx} =  avgPrefMax.*(avgPrefMax>( mean(Frames(:,:,t_on==-1),3) + 4*std(Frames(:,:,t_on==-1),[],3) ));
    pFrame{aIdx} = (avgPostONframe - avgPostOFFframe);
    
    bFrame{aIdx} = (avgPostONframe2 - avgPostOFFframe2);
   
    maxFrame(:,:,aIdx) = avgPrefMax;
    
    se = strel('square',2);
    binaryFrame{aIdx} = imfill(imclose(aFrame{aIdx}>0,se),'holes')  - imfill(imclose(pFrame{aIdx}>3*(std(std([pFrame{:}]))),se),'holes');

end

% And keep only positive values (negative value means less active
% during trials than before/after)
% aFrame(aFrame<0) = 0;

if nargout >0
    varargout{1} = binaryFrame;%
    varargout{2} = angles;
    varargout{3} = maxFrame;
else
    % Push to object
    for aIdx = 1:length(angles)
        obj.(['ActivityFrame' num2str(aIdx)]) = aFrame{aIdx};
    end
    %%
    %and plot
    figure('color','w');
    for aIdx = 1:length(angles)
        subplot(1,length(angles),aIdx)
        imagesc(bFrame{aIdx})
        polarmap
        
        title([num2str(angles(aIdx))])
        axis image
        grid on
        set(gca,'XTickLabel',{});
        set(gca,'YTickLabel',{});
        ax(aIdx) = gca;
        
    end
    % link cmap lims
    for aIdx = 1:length(ax)
        ax(aIdx).CLim = [ min([ax.CLim]) max([ax.CLim])];
        
    end
    
end
%{
    figure('color','w');
for aIdx = 1:length(angles)
    subplot(3,length(angles),aIdx)
%     imagesc(aFrame{aIdx}>0 + -(pFrame{aIdx}>0))
se = strel('square',2);
     imagesc(imclose((aFrame{aIdx}>0).*(median(obj.Frames,3)>10),se)  )
    polarmap
%     colormap(flipud(hot))
    title(['trial ' num2str(angles(aIdx)) ])
    axis image
    grid on
    set(gca,'XTickLabel',{});
    set(gca,'YTickLabel',{});
    ax(aIdx) = gca;
    
        subplot(3,length(angles),length(angles)+aIdx)
se = strel('square',2);
     imagesc(-imclose(pFrame{aIdx}>( 3*(std(std([pFrame{:}]))) ),se))
    polarmap
%     colormap(flipud(hot))
    title(['post ' num2str(angles(aIdx))])
    axis image
    grid on
    set(gca,'XTickLabel',{});
    set(gca,'YTickLabel',{});
    
        subplot(3,length(angles),2*length(angles)+aIdx)
se = strel('square',2);
 imagesc(imclose(aFrame{aIdx}>0,se)  - 2*imclose(pFrame{aIdx}>3*(std(std([pFrame{:}]))),se))
    polarmap
%     colormap(flipud(hot))
    title([num2str(angles(aIdx))])
    axis image
    grid on
    set(gca,'XTickLabel',{});
    set(gca,'YTickLabel',{});

end

%}