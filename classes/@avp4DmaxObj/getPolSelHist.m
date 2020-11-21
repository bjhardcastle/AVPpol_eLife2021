function counts = getPolSelHist(obj)
counts = nan(1,length(obj.polSelHistBins)-1);

    
    maskedPolSelFrameVector = obj.polSelTopVector(obj);
    
    counts = histcounts(maskedPolSelFrameVector,obj.polSelHistBins, 'Normalization', 'probability');

end