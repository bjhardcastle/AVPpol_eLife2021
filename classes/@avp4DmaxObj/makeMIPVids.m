function makeMIPVids(obj,exps)
% input MIP object
assert(contains(obj.File,'_Zmax_'),'Function runs on max intensity projection tiffs.')
% run getparameters if seqnum empty
% TrialSeqNum = exp
% TrialPatNum = trial info
if isempty(obj.TrialSeqNum)
    getParameters(obj)
end
if isempty(obj.Daq)
    getDaqData(obj)
end
if isempty(obj.Frames)
    getFrames(obj)
end
    runBackSub(obj)

if nargin < 2 || isempty(exps)
    exps = unique(obj.TrialSeqNum);
end

% Make crop box by thresholding activity across all frames

activityVal = nanmean(obj.Frames(:)) ;
% activityVal = nanmean(obj.AverageFrame(:)) + 0.5*nanstd(obj.AverageFrame(:)) ;
activityBinaryImg = obj.AverageFrame > activityVal;
% figure,imagesc(activityBinaryImg)

% Dectect areas of activity
bw = regionprops(activityBinaryImg);
% Sort by size
[~,sortIdx] = sort([bw.Area],'descend');
% Find the bounding box of the largest detected region
cropBox = bw(sortIdx(1)).BoundingBox;

for expNum = exps
    if any(ismember(obj.TrialSeqNum,expNum))
        switch expNum
            
            case {1,3}   % ON/OFF UV
                loops = 1;
                [frames,stim] = getEveryFrame(obj,expNum);
                framerate = 8;
                
            case {2,4,8} % Pol angles
                loops = 2;
                [frames,stim] = getTrialFrame(obj,expNum);
                framerate = 8;
                
            case 5       % Bars
                loops = 1;
                [frames,stim] = getEveryFrame(obj,expNum);
                framerate = 8;
                
            case 6       % OF expansion
                loops = 1;
                [frames,stim] = getEveryFrame(obj,expNum);
                framerate = 8;
                
            case {7,9}       % Blue / UV ON/OFF
                loops = 1;
                [frames,stim] = getEveryFrame(obj,expNum);
                framerate = 8;
        end
%         
        makeVideo(obj,frames,stim,expNum,loops,framerate,cropBox,0,0)
        makeVideo(obj,frames,stim,expNum,loops,framerate,cropBox,1,0)
        makeVideo(obj,frames,stim,expNum,loops,framerate,[],0,0)
        makeVideo(obj,frames,stim,expNum,loops,framerate,[],1,0)
        makeVideo(obj,frames,stim,expNum,loops,framerate,cropBox,0,1)
        makeVideo(obj,frames,stim,expNum,loops,framerate,[],0,1)
    end
end

end

function [frames,stim] = getEveryFrame(obj,expNum)
fields = [];
fields.TrialSeqNum = expNum;
trialIdx = findTrials(obj,fields);

startFrame = obj.TrialStartFrame(trialIdx(1)) - 2;
endFrame = obj.TrialEndFrame(trialIdx(end)) + 2;

if endFrame > size(obj.Frames,3)
    endFrame = size(obj.Frames,3);
end

frameSamples = obj.Frametimes([ startFrame : endFrame ]);

stim = obj.Daq(frameSamples,2);


frames = obj.Frames(:,:,startFrame:endFrame);
end

function [frames,stim] = getTrialFrame(obj,expNum)

fields = [];
fields.TrialSeqNum = expNum;
trialIdx = findTrials(obj,fields);

for tidx = 1:length(trialIdx)
    
    trialVec = trialIdx(tidx);
    
    avgFrame = [];
    for vidx = 1:length(trialVec)
        startFrame = obj.TrialStartFrame(trialVec(vidx));
        endFrame = obj.TrialEndFrame(trialVec(vidx));
        if endFrame > size(obj.Frames,3)
            endFrame = size(obj.Frames,3);
        end

        avgFrame(:,:,vidx) = mean(obj.Frames(:,:,startFrame:endFrame),3);
    end
    
    frames(:,:,tidx) = mean(avgFrame,3);
    
    stim(tidx) = obj.TrialPatNum(trialVec);
    
end


end

function [frames,stim] = getAveragedFrame(obj,expNum)
fields = [];
fields.TrialSeqNum = expNum;
trialIdx = findTrials(obj,fields);

uniqueTrials = unique(obj.TrialPatNum(trialIdx));

if expNum == 2 || expNum == 4 || expNum == 8
    for tidx = 1:0.5*length(uniqueTrials)
        trialVec = [];
        
        fields = [];
        fields.TrialSeqNum = expNum;
        fields.TrialPatNum = uniqueTrials(tidx);
        t1 = findTrials(obj,fields);
        
        fields = [];
        fields.TrialSeqNum = expNum;
        fields.TrialPatNum = uniqueTrials(tidx+0.5*length(uniqueTrials));
        t2 = findTrials(obj,fields);
        
        trialVec = sort([t1;t2]);
        
        avgFrame = [];
        for vidx = 1:length(trialVec)
            startFrame = obj.TrialStartFrame(trialVec(vidx));
            endFrame = obj.TrialEndFrame(trialVec(vidx));
            if endFrame > size(obj.Frames,3)
                endFrame = size(obj.Frames,3);
            end
            avgFrame(:,:,vidx) = mean(obj.Frames(:,:,startFrame:endFrame),3);
        end
        
        frames(:,:,tidx) = mean(avgFrame,3);
        
        stim(tidx) = uniqueTrials(tidx);
        
    end
end

end
% make two videos:
%     1) realtime
%     2) average frame for each trial
% Store angle values in vector
% Add stim marker

function makeVideo(obj,frames,stim,expNum,loops,framerate,cropBox,markerONOFF,scaled)
% Get path for saving video
if cropBox
    if markerONOFF
        saveFile = [num2str(expNum),'_markerON_cropped'];
    else
        saveFile = [num2str(expNum),'_markerOFF_cropped'];
    end
else
    if markerONOFF
        saveFile = [num2str(expNum),'_markerON'];
    else
        saveFile = [num2str(expNum),'_markerOFF'];
    end
end
if scaled 
   saveFile = [saveFile '_scaled']; 
end

savePath = obj.Folder;
% Make video file and open
v = VideoWriter(fullfile(savePath,saveFile),'Motion JPEG AVI');
% v.FrameRate = round(obj.AIrate/obj.IFI);
v.FrameRate = framerate;
v.Quality = 100;
open(v)


% Open figure for displaying/capturing frames
VidFig = figure('Color',[0 0 0]);
%;
% set(VidFig,'units','normalized' ...
%     , 'outerposition',[0 0 0.5 0.5]  );
ax = gca;

for n = 1:loops
    for frameidx = 1:length(stim)
        cla(ax)
        
        % Normalize image to brightest value across all frames
        maxVal = mean(frames(:)) + 6*std(frames(:));
        normFrame = frames(:,:,frameidx)./maxVal;
        
        if ~isempty(cropBox)
            normFrame = imcrop(normFrame,cropBox);
        end
        
        % Rotate frames by 90deg cw
        fFrame = imrotate(normFrame,90);
%         fFrame = imresize(fFrameOrig,20*(obj.micronsPerPixel));
imshow(fFrame,'Parent',ax)
axis equal
axis off

set(ax,'position',[0 0 1 1],'units','normalized')
if scaled
    VidFig.Position = [0 0 5*size(fFrame,2)*obj.micronsPerPixel, 5*size(fFrame,2)*obj.micronsPerPixel ];

else
VidFig.Position = [0 0 2*floor((1080*size(fFrame,2)/size(fFrame,1))/2),1080 ];
end
%         caxis(ax,[0 0.3])
        
        if markerONOFF
            
            w = size(fFrame,2);
            wD =w/8;
            wSpace = wD/10;
            wRad = wSpace*8;
            
            switch expNum
                
                case {1,3}  % UV ON/OFF
                    
                    % Add circle for light on off indicator
                    circDia = wRad;
                    xyBottomLeft = [wD + wSpace,size(fFrame,1)-wD+wSpace];
                    % Remember, y-axis is inverted, so bottom left is actually
                    % top left, and the +ve height will extend downwards from there
                    pos = [xyBottomLeft(1) xyBottomLeft(2) circDia circDia];
                    R = rectangle('Position',pos,'Curvature',[1 1]);
                    
                    % Change color depending on stim marker
                    markerPos = (round(stim(frameidx)*5)+1)/2;
                    if markerPos >= 1
                        R.FaceColor = [0.75,0.1,1];
                        R.EdgeColor = 'none';
                    else
                        R.FaceColor = 'none';
                        R.EdgeColor = 'none';
                        
                    end
                    
                    
                    if expNum == 3
                        % Also draw line to indicate pol angle
                        if frameidx <= 0.5*length(stim)
                            aop = 90;
                        else
                            aop = 180;
                        end
                        lineRadius = wRad/2;
                        xyCenter = [wD/2,size(fFrame,1)-wD/2];
                        lineX = lineRadius*cosd(aop);
                        lineY = -lineRadius*sind(aop);
                        xyEnd1 = xyCenter - [lineX, lineY];
                        xyEnd2 = xyCenter + [lineX, lineY];
                        
                        %             L = line([xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                        L = line([xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                        
                        L.Color = 'w';
                        L.LineWidth = 2*wSpace;
                        
                    end
                    
                    
                    
                    
                case {2,4,8} % Pol
                    aop = stim(frameidx);
                    % Draw AoP indicator
                    lineRadius = wRad/2;
                    xyCenter = [wD/2,size(fFrame,1)-wD/2];
                    lineX = lineRadius*cosd(aop);
                    lineY = -lineRadius*sind(aop);
                    xyEnd1 = xyCenter - [lineX, lineY];
                    xyEnd2 = xyCenter + [lineX, lineY];
                    
                    %             L = line([xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                    L = line([xyEnd1(1) xyCenter(1)], [xyEnd1(2) xyCenter(2)]);
                    
                    L.Color = 'w';
                    L.LineWidth = 2*wSpace;
                    
                case 5 % Bars
                    
                    wBox = 2*wD + 2*wSpace;
                    xyBox =  [wBox,wBox/3];
                    lineX1 = linspace(2*wSpace,2*wSpace + wBox,5);
                    lineX2 = lineX1;
                    lineY1 = size(fFrame,1) - repmat(wSpace,1,5);
                    lineY2 = lineY1 - xyBox(2);
                    
                    % Add grey bars for stimulus positions
                    for n = 1:5
                        L(n) = line(ax,[lineX1(n) lineX2(n)], [lineY1(n) lineY2(n)]);
                        L(n).Color = [0.3 0.3 0.3];
                        L(n).LineWidth = round(0.75*mean(diff(lineX1)));
                    end
                    % Add a single blue bar to mark current position
                    markerPos = (round(stim(frameidx)*5)+1)/2;
                    if markerPos >= 1
                        L(markerPos).Color = [0.15,0.1,1];
                    end
                    
                    
                case 6 % OF
                    wBox = 2*wD + 2*wSpace;
                    xyBox =  [wBox,wBox/3];
                    lineX1 = linspace(2*wSpace,2*wSpace + wBox,5);
                    lineX1([1,5]) = [];
                    lineX2 = lineX1;
                    lineY1 = size(fFrame,1) - repmat(wSpace,1,3);
                    lineY2 = lineY1 - xyBox(2);
                    
                    % Add grey bars for stimulus positions
                    for n = 1:3
                        L(n) = line(ax,[lineX1(n) lineX2(n)], [lineY1(n) lineY2(n)]);
                        L(n).Color = [0.3 0.3 0.3];
                        L(n).LineWidth = 2*ceil(mean(diff(lineX1)));
                    end
                    % Add a single blue bar to mark current position
                    markerPos = round(stim(frameidx)*5);
                    switch markerPos
                        case {1,2,3}
                            L(markerPos) = line(ax,[lineX1(markerPos) lineX2(markerPos)], [lineY1(markerPos) lineY2(markerPos)]);
                            L(markerPos).LineWidth = 2*floor(mean(diff(lineX1)));
                            L(markerPos).Color = [0.15,0.1,1];
                        case 4
                            L(1).Color = [0.15,0.1,1];
                            L(2).Color = [0.15,0.1,1];
                            L(3).Color = [0.15,0.1,1];
                    end
                    
                case {7,9}  % UV ON/OFF
                    
                    % Two versions of exp7 exist:
                    %  1) early versions with 3 repetitions
                    %  2) later versions with 3 different trials, 3 reps
                    %  each.
                    % switch code based on number of trials
                    
                    switch length(find([obj.TrialSeqNum]== expNum))
                        
                        case 3
                            % Add circle for blue light on off indicator
                            circDia = wRad;
                            xyBottomLeft = [2*wD + wSpace,size(fFrame,1)-wD+wSpace];
                            
                            pos = [xyBottomLeft(1) xyBottomLeft(2) circDia circDia];
                            R = rectangle('Position',pos,'Curvature',[1 1]);
                            
                            % Change color depending on stim marker
                            markerPos = round(stim(frameidx)*5);
                            if markerPos >= 1
                                R.FaceColor = [0.15,0.1,1];
                                R.EdgeColor = 'none';
                            else
                                R.FaceColor = 'none';
                                R.EdgeColor = 'none';
                            end
                            
                            
                            
                        case 9
                            
                            % Add circle for blue light on off indicator
                            circDia = wRad;
                            xyBottomLeft = [2*wD + wSpace,size(fFrame,1)-wD+wSpace];
                            
                            pos = [xyBottomLeft(1) xyBottomLeft(2) circDia circDia];
                            Rb = rectangle('Position',pos,'Curvature',[1 1]);
                            
                            % Change color depending on stim marker
                            markerPos = round(stim(frameidx)*5);
                            if markerPos >= 1 && frameidx <=(floor(2*length(stim)/3))
                                Rb.FaceColor = [0.15,0.1,1];
                                Rb.EdgeColor = 'none';
                            else
                                Rb.FaceColor = 'none';
                                Rb.EdgeColor = 'none';
                            end
                            
                            % Add circle for UV light on off indicator
                            % Add circle for light on off indicator
                            circDia = wRad;
                            xyBottomLeft = [wD + wSpace,size(fFrame,1)-wD+wSpace];
                            % Remember, y-axis is inverted, so bottom left is actually
                            % top left, and the +ve height will extend downwards from there
                            pos = [xyBottomLeft(1) xyBottomLeft(2) circDia circDia];
                            Ruv = rectangle('Position',pos,'Curvature',[1 1]);
                            
                            % Change color depending on stim marker
                            markerPos = round(stim(frameidx)*5);
                            if markerPos >= 1 && frameidx >=(floor(length(stim)/3))
                                Ruv.FaceColor = [0.75,0.1,1];
                                Ruv.EdgeColor = 'none';
                            else
                                Ruv.FaceColor = 'none';
                                Ruv.EdgeColor = 'none';
                            end
                            
                            
                            if expNum == 9 && frameidx >=(floor(length(stim)/3))
                                
                                
                                % Also draw line to indicate pol angle
                                aop = 0;
                                lineRadius = wRad/2;
                                xyCenter = [wD/2,size(fFrame,1)-wD/2];
                                lineX = lineRadius*cosd(aop);
                                lineY = -lineRadius*sind(aop);
                                xyEnd1 = xyCenter - [lineX, lineY];
                                xyEnd2 = xyCenter + [lineX, lineY];
                                
                                %             L = line([xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                                L = line([xyEnd1(1) xyEnd2(1)], [xyEnd1(2) xyEnd2(2)]);
                                
                                L.Color = 'w';
                                L.LineWidth = 2*wSpace;
                                
                            end
                            
                    end
            end
            
        end
        
        % Grab plot frame and write to video object
        frame = getframe(gcf);
        writeVideo(v,frame);
        
    end
end
close(VidFig)
close(v)


end