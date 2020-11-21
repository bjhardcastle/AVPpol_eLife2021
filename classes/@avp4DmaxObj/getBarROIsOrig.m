function varargout = getBarROIsOrig(obj)
numROIs = 1;
f = mean(obj.Frames,3);
[polarImage, patNum] = getExp5ActivityFrame(obj);
polColMap = flipud(hsv(numROIs*2*length(patNum)));

tempROI = struct('barpos',[],'polarity',[],'mask',[],'color',[],'position',[],'response',[]);
 %%   
ridx = 0;

for pIdx = 1:length(polarImage)
    
     u = [-1 1]; % values within k
               

    for l = 1:length(u)
  k = polarImage{pIdx};
        ridx = ridx+1;
        
        col = polColMap(ridx,:);
        
        rk = k==u(l);
%         maskcol = squeeze(cat(3,col(1)*rk,col(2)*rk,col(3)*rk));
        
          fInv =  (1-(f./max(f(:))));
    % Create ~mask, setting pixels that will become colored, polSel
    % regions to zero (since we will add values for color)
    fInvMask = fInv.*(k~=u(l));
    fullMasks = cat(3,fInvMask,fInvMask,fInvMask);
    examplemasks = cat(3,fInvMask,fInvMask,fInvMask);

%         fullMasks = fullMasks + maskcol;
        
            % Get the pixels included at this angle
            rk = k==u(l);
            rp = regionprops(rk,'PixelIdxList','Area');
            [~,rpSrtIdx] = sort([rp.Area],'descend');
            % Sort by area size and for the top X areas, get the idx of
            % all pixels contained within:
            rp = rp(rpSrtIdx);
            % Group the top numROIs areas by size:
            %                 %{
            if isempty(rp)
                ridx = ridx -1;
                continue
            elseif length(rp) < numROIs
                pixForColor = vertcat(rp(1:end).PixelIdxList);
            else
                pixForColor = vertcat(rp(1:numROIs).PixelIdxList);
            end
            %                 %}
            % -- changed to take the top area, only:
            %                 pixForColor = vertcat(rp(1).PixelIdxList);
            
            % Get the idx of all pixels currently set to zero in
            % 'masks', awaiting color (because not all will be colored
            % now)
            pixAtZero = find(k==u(l));
            % Find the pixels which won't be colored and should be
            % reset to original average image (inverted)
            pixToReset = setdiff(pixAtZero,pixForColor);
            
            % Update polAng array to discard all areas but the top X
            k = zeros(size(k));
            k(pixForColor) = u(l);
%             k(pixToReset) = 0;
            
            % Update masks array to discard all areas but the top X
            tempMask = zeros( size(f,1),size(f,2));
            tempMask(pixToReset) = 1;
            toAdd = tempMask.*fInv;
            toAdd3 = cat(3,toAdd,toAdd,toAdd);
            examplemasks = examplemasks + toAdd3;
        
        rk = k==u(l);
        maskcol = squeeze(cat(3,col(1)*rk,col(2)*rk,col(3)*rk));
        examplemasks = examplemasks + maskcol;
          
%         figure
%     imagesc(examplemasks)
    
        maskResponse = scanROI(obj,rk);
                    
        tempROI(ridx).mask = rk;
        tempROI(ridx).color = col;
        tempROI(ridx).position = [1 1 1 1];
        tempROI(ridx).response = maskResponse;
        tempROI(ridx).barpos = pIdx;
        tempROI(ridx).polarity  = u(l);
        
        maskAreas(ridx) = length(find(rk));
    end
  

end
% [~,sortIdx] = sort(maskAreas,'descend');
% tempROI = tempROI(sortIdx);
% % autoROI = tempROI([ find([tempROI.polarity] == 1,2,'first') find([tempROI.polarity] == -1,2,'first') ]); 
autoROI =  tempROI([ find([tempROI.polarity] == 1) find([tempROI.polarity] == -1) ]);
if nargout > 0
   varargout{1} = autoROI;
end
    
%%
end


