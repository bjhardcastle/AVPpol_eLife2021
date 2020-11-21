function tempPlot(obj)
if isempty(obj.polAngRespPooled)
    getPolMaps(obj);    
end

figure
cols = flipud(magma);    

resp0map = getF0map(obj);
selmap = obj.polSelImg;
magmap = obj.fftMagImg;

mag = magmap./max(magmap(:));
w = mag(:).*selmap(:);

pst = obj.polSelThreshold;
obj.polSelThreshold = 0.2;

polAngShift = obj.polAngRespPooled;
for n = 1:6
       
end

for n = 1:12
    hold on 
    
    
    resp = obj.polAngResp(obj.polPixDiscrete==n*30,:);
    
    
    if ~isempty(resp)
        bkgimg = isnan(obj.polPix).*(obj.layerMask.mask);
        
        bkg = obj.polAngResp(logical(bkgimg(:)),:);
        
        
        [~,sortIdx] = sort(max(bkg,[],2),'ascend');
        resp0 = mean(bkg(sortIdx,:),1);
%         respNorm = ((resp - resp0)./mean(resp0,2)) - 1;
        respNorm = ((resp )./mean(resp0,2)) - 1;

%         plot([30:30:360],respNorm,'c')
        
   
        weights=w(obj.polPixDiscrete==n*30);
  wMean  = mean(weights);
  
  
respNorm = circshift(respNorm,-n-2,2);

                plot([30:30:360],mean(respNorm,1),'color',cols(round(wMean*255)+1,:),'LineWidth',2);

%     resp = obj.polAngRespPooled(obj.polPixDiscrete==n*30,:);
%     resp0 = resp0map(obj.polPixDiscrete==n*30);
%     weights=w(obj.polPixDiscrete==n*30);
%      [~,sortIdx] = sort(weights,'descend');
%     if length(sortIdx) >= 50
%         topIdx = sortIdx(1:50);
%     else
%         topIdx = sortIdx;
%     end     
%     topIdx = topIdx(~isnan(w(topIdx)));
%     respN = (resp(topIdx,:)./resp0(topIdx))-1;
%     respMax = max(respN,[],2);
%     respMin = min(respN,[],2);
%     respNorm = respN - respMin;
%     respNorm = respNorm./max(respNorm,[],2);
%      wMean = mean(weights(topIdx));
% respNorm = circshift(respNorm,-n-2,2);
%     if ~isempty(respNorm)
%         lineprops.col = {cols(round(wMean*255)+1,:)};
%         
%         mseb([30:30:360],mean(respNorm,1),std(respNorm,[],1),lineprops,1);
%     end


 end
    
%     plot([30:30:360],respNorm,'c')
%     [~,sortIdx] = sort(max(resp,[],2),'ascend');
%      resp0 = mean(resp(sortIdx(1:floor(length(sortIdx)/10)),:),1);
%     respNorm = resp./resp0 - 1;

end

obj.polSelThreshold = pst;