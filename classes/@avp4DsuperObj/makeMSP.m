function makeMSP(obj)
% After selecitvity maps are available for all obj.Layers we can construct
% something similar to a max. intensity projection, but using the layer
% with the highest pol selectivity for each pixel, and using its
% information for each timepoint
% (Not currently used)

if length(obj.Layers) == 1
    obj.MIP = obj.Layers(1);
    disp('Only one layer exists: MSP is from a single plane')
    return
end
if any(obj.Exps==2|obj.Exps==4|obj.Exps==8|obj.Exps==10)
    if isempty(obj.MIP.Frametimes) || isempty(obj.MIP.TrialPatNum)
        getParameters(obj.MIP)
    end
    if isempty(obj.Layers(1).Frames)
        getFrames(obj.Layers(1));
    end
    if length(obj.MIP.Frametimes) > size(obj.Layers(1).Frames,3)
            getLogData(obj.MIP)  
            numFrames = obj.MIP.logTimePoints;
        else
            numFrames = length(obj.MIP.Frametimes);
        end
    mspSelectStack = nan(size(obj.Layers(1).AverageFrame,1),size(obj.Layers(1).AverageFrame,2),length(obj.Layers));
    mspLayerMap = nan(size(obj.Layers(1).AverageFrame,1),size(obj.Layers(1).AverageFrame,2));
    mspFrames = nan(size(obj.Layers(1).AverageFrame,1)*size(obj.Layers(1).AverageFrame,2),numFrames);
    for Lidx = 1:length(obj.Layers)
        if isempty(obj.Layers(Lidx).Frames)
            getFrames(obj.Layers(Lidx));
        end
        
        % Collect selectivity maps from each layer
        getPolMaps(obj.Layers(Lidx))
        
        if ~isfield(obj.pSet(obj.Layers(1).polExp),'trialRandomizeOrder') || obj.pSet(obj.Layers(1).polExp).trialRandomizeOrder == 0
            % regular exp
            mspSelectStack(:,:,Lidx) = abs(obj.Layers(Lidx).polSelImg.*obj.Layers(Lidx).fftMagImg);
        else
            % randomized angle exp: can't use fft signal
            mspSelectStack(:,:,Lidx) = abs(obj.Layers(Lidx).polSelImg);
        end
    end
    
    % Find the position of the max selectivity across layers
    [~,mspLayerMap] = max(mspSelectStack,[],3);
    
    % Now go through Layer objects and extract time-series across Frames from
    % all pixels included in MSP
    s1 = size(obj.Layers(1).AverageFrame,1);
    s2 = size(obj.Layers(1).AverageFrame,2);
    tempFrame = nan(s1*s2,numFrames);
    for Lidx = 1:length(obj.Layers)
      
      
        tempFrame = reshape(obj.Layers(Lidx).Frames(:,:,1:numFrames),size(tempFrame,1),size(tempFrame,2));
        mspFrames(mspLayerMap==Lidx,:) = tempFrame(mspLayerMap==Lidx,:);
        
        obj.Layers(Lidx).Frames = [];
    end
    mspFrames = reshape(mspFrames,s1,s2,[]);
    obj.MIP.Frames = mspFrames;
    superUseMSP(obj,1);

end
end