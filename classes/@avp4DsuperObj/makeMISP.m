function makeMISP(obj)
% After selecitvity maps are available for all obj.Layers we can construct
% something similar to a max. intensity projection, but using the layer
% with the highest activity for each pixel, and using its
% information for each timepoint
% (Currently used)

superUseMSP(obj,0)
if length(obj.Layers) == 1
    obj.MIP = obj.Layers(1);
    disp('Only one layer exists: MSP is from a single plane')
    return
end
% if any(obj.Exps==2|obj.Exps==4|obj.Exps==8|obj.Exps==10)
     %{
     % removedApril 28
    if isempty(obj.MIP.Frametimes) || isempty(obj.MIP.TrialPatNum)
        getParameters(obj.MIP)
    end

   
    if length(obj.MIP.Frametimes) > size(obj.Layers(1).Frames,3)
            getLogData(obj.MIP)
            numFrames = obj.MIP.logTimePoints;
        else
            numFrames = length(obj.MIP.Frametimes);
        end
    %}
if isempty(obj.Layers(1).Frames)
    getFrames(obj.Layers(1));
end
s1 = size(obj.Layers(1).AverageFrame,1);
s2 = size(obj.Layers(1).AverageFrame,2);

mspSelectStack = nan(s1,s2,length(obj.Layers));
    for Lidx = 1:length(obj.Layers)
        obj.Layers(Lidx).Unattended = 1;
        if Lidx > 1
            getFrames(obj.Layers(Lidx));
        end
        % force refresh some key parameters
        getDaqData(obj.Layers(Lidx))
        getFrametimes(obj.Layers(Lidx))
        getTrialtimes(obj.Layers(Lidx))
        getParameters(obj.Layers(Lidx))

        getInactivityFrame(obj.Layers(Lidx));
        
        % Collect inactivity maps from each layer
        mspSelectStack(:,:,Lidx) = abs(obj.Layers(Lidx).InactivityFrame);
        
        obj.Layers(Lidx).Daq = [];
    end
    num = @(c) size(c,3);
    frameDims = cellfun(num,{obj.Layers.Frames},'UniformOutput',0);
    numFrames = min([frameDims{:}]);
    
    mspFrames = nan(size(obj.Layers(1).AverageFrame,1)*size(obj.Layers(1).AverageFrame,2),numFrames);
    bkgInt = [];

    % Find the position of the max activity across layers
    [~,mspLayerMap] = max(mspSelectStack,[],3);
    
    % Now go through Layer objects and extract time-series across Frames from
    % all pixels included in MSP. Reshape into 2D pixel coords into 1D to
    % make indexing easier
    tempFrame = nan(s1*s2,numFrames);
    for Lidx = 1:length(obj.Layers)
        % reshape each layer's frames
        tempFrame = reshape(obj.Layers(Lidx).Frames(:,:,1:numFrames),size(tempFrame,1),size(tempFrame,2));
        % Extract pixel time-series where they are the maximum avg int
        % across the stack
        mspFrames(mspLayerMap==Lidx,:) = tempFrame(mspLayerMap==Lidx,:);
        obj.Layers(Lidx).Frames = [];
    end
    % Reshape resulting new projection
    mspFrames = reshape(mspFrames,s1,s2,[]);
    obj.MIP.Frames = mspFrames;
    getInactivityFrame(obj.MIP);
    superUseMSP(obj,1)   

end
% end