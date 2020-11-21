function superImageSummary(objarray)
discretize=0;
for oidx = 1:length(objarray)
    
    if isempty(objarray(oidx).MIP.polExp)
        continue
    else
        obj = objarray(oidx).MIP;
        titleStr = ['obj(' num2str(oidx) ') ' objarray(oidx).Name];
        figure('Name',titleStr);
        ax(1) = subplot(1,2,1);
        imagesc(ax(1),obj.polPix)
        colormap(ax(1),flipud(hsv(180)))
        
        avgFrame = 1-sqrt(obj.avgPolImg./max(obj.avgPolImg(:)));
        tuningFrame = obj.polPix;
        polAngImg = obj.polTuningImg;
        
        tuningFrame = polAngImg;
        tuningFrame(obj.polSelImg<obj.polSelThreshold) = nan;
        combFrame = zeros(size(avgFrame)    );
        if discretize
            combFrame(~isnan(tuningFrame)) = polAngImg(~isnan(tuningFrame));
        else
            combFrame(~isnan(tuningFrame)) = obj.polTuningImg(~isnan(tuningFrame));
        end
        combFrame(isnan(tuningFrame)) = -180.*avgFrame(isnan(tuningFrame));
        
        layer_im = imagesc(imrotate(combFrame,90), 'Parent', ax(1));
%         colorbar
        
        colormap(ax(1),[flipud(gray(180));[1 1 1];flipud(hsv(180))])
        
        set(ax(1),'CLim',[-180 180])
%         axis image
        axis off
        if ~isempty(objarray(oidx).MIP.layerMask) && ~isempty(objarray(oidx).MIP.layerMask.position)
            hold on
            maskPosX = [objarray(oidx).MIP.layerMask.position(:,2) ; objarray(oidx).MIP.layerMask.position(1,2)];
            maskPosY = 1+size(combFrame,2)-[objarray(oidx).MIP.layerMask.position(:,1); objarray(oidx).MIP.layerMask.position(1,1)];

		plot(maskPosX,maskPosY,'Color','k','LineWidth',2)
        end

        ax(2) = subplot(1,2,2);
        imagesc(ax(2),imrotate(objarray(oidx).MIP.polSelImg,90))
        colormap(ax(2),flipud(magma))
%         axis image
        axis off
                set(ax(2),'CLim',[0 1])

%         colorbar
        
        
        
        tightfig(gcf);
        
        if ~isempty(objarray(oidx).MIP.layerMask.position)
            hold on
		plot(maskPosX,maskPosY,'Color','w','LineWidth',2)
        end
    end
    
end
end