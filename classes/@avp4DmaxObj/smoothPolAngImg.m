function connMask = smoothPolAngImg(obj,polAngImg)

polAngImgCopy = polAngImg;
connMask = zeros(size(polAngImg));

minArea = 26;

uniqueAng = unique(polAngImg(polAngImg>0));
for n = 1:length(uniqueAng)
           I = polAngImgCopy == uniqueAng(n);
           connMask = connMask +  uniqueAng(n).*bwareaopen(I,minArea);
end
% connMaskFill = imfill(connMask);
% 
% invAng = flipud(uniqueAng);
% invFill = zeros(size(polAngImg));
% probPix = zeros(size(polAngImg));
% for n = 1:length(uniqueAng)
%         invFill(connMaskFill == uniqueAng(n)) = invAng(n);
%         probPix(connMaskFill == uniqueAng(n)) = probPix(connMaskFill == uniqueAng(n)) +1;
% end
%    revFill = imfill(invFill);
%    polAngImgFill= zeros(size(polAngImg));
% for n = 1:length(uniqueAng)             
%         polAngImgFill(revFill == invAng(n)) = uniqueAng(n);
%          probPix(revFill == invAng(n)) = probPix(revFill == invAng(n)) +1;
% end

figure,subplot(1,3,1), imshow(polAngImgCopy./6);
subplot(1,3,2), imshow(connMask./6);
subplot(1,3,3),imshow(polAngImgCopy.*(polAngImgCopy & ~connMask)./6)

smallAreas = polAngImgCopy.*(polAngImgCopy & ~connMask);
plusminus1Ang(:,1) = circshift(uniqueAng,1);
plusminus1Ang(:,2) = circshift(uniqueAng,-1);
for n = 1:length(uniqueAng)
           I = smallAreas == uniqueAng(n);
           rp = regionprops(I,'PixelIdxList');
                          
           % Examine each small area individually
           for rIdx = 1:length(rp)
               % Get the surrounding perimeter of the area 
               areaPix = zeros(size(polAngImg));
               areaPix(rp(rIdx).PixelIdxList) = 1;
               perimPix = find(imdilate(areaPix, true(3)) - areaPix);
               perimPixVals = polAngImgCopy(perimPix);
               
               % If all surrounding pixels are 0, it's a small area of
               % noise, too small to be a synapse: set to zero.
               if all(perimPixVals == 0)                  
                   connMask(rp(rIdx).PixelIdxList) = 0;
                   
                   % If non-zero surrounding pixels are within one
                   % angle-step, they are small noisy areas likely to
                   % belong to the surrounding area or vice versa:                  
               elseif any( mode(perimPixVals(perimPixVals>0)) == plusminus1Ang(n,:))
                   % Find area which is larger 
                   innerSize = length(rp(rIdx).PixelIdxList);                   
                   [arow,acol] = ind2sub(size(polAngImg),perimPix( perimPixVals==mode(perimPixVals(perimPixVals>0)) ));
                   outerPix = find(grayconnected(polAngImgCopy,acol(1),arow(1)));
                   outerSize = length(outerPix);
                   
                   if outerSize > innerSize
                       connMask(rp(rIdx).PixelIdxList) = mode(perimPixVals(perimPixVals>0));
                   else
                       connMask(outerPix) = uniqueAng(n);
                   end
                   % If non-zero pixels are 
               else
                   % If non-zero surrounding pixels are NOT within one
                   % angle-step, they are small noisy areas which we will
                   % not include in a mask:
                   connMask(rp(rIdx).PixelIdxList) = 0;
               end
           end
end

% 
% maskCounter = 0;
% for n = 1:length(uniqueAng)
%     continueFlag = 1;
%     tempCounter = maskCounter;
%     while continueFlag
%         
%         I = polAngImgCopy == uniqueAng(n);
%         rp = regionprops(imfill(connMask),'Area','Solidity','PixelIdxList','EulerNumber','FilledImage','Image');
% 
%         % Make a new regionprops struct containing only areas larger than
%         % threshold
%         r = rp([rp.Area] >= minArea);
%         
%         %%
%         % Areas less than threshold are set to zero in new polAngImg
%         % (later, if these areas are 'holes' they will be consumed by other
%         % areas, regardless of Ang value)
%         polAngImgCopy(vertcat(rp([rp.Area] < minArea).PixelIdxList)) = 0;
% 
%         [~,sortIdx] = sort([r.Area],'descend');
%         r = r(sortIdx);
%         
%         for m = 1:length(r)
%             if r(m).EulerNumber ==  1
%                 % area is inserted into polAngImg2
%                 tempCounter = tempCounter + 1;
%                 connMask(r(m).PixelIdxList) = tempCounter;
%             else
%                 polAngImgCopy = resolveHoles(polAngImgCopy,r(m));
%             end
%         end
%         
%         continueFlag = ~all([r.EulerNumber]==1);
%         
%     end
%     maskCounter = tempCounter;
% end
% end
% 
% function polAngImgCopy = resolveHoles(polAngImgCopy,rm)
%     
% end