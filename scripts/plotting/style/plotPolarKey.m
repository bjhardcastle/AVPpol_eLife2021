figure('color','w','Name','Polar Key','NumberTitle','off','position',[ 215,575,265,243])
pAx = polaraxes;
hold(pAx,'on');
    set(pAx,'ThetaDir','clockwise','ThetaZeroLocation','right')
    pAx.RAxis.Visible = 'off';
    pAx.RLim = [0 1];
    pAx.RTick = [0 1];
    
    
    
    numAngs = 6;
    cols = flipud(hsv(numAngs));
for rIdx = 1:numAngs
            % add line to polar axes
            theta = deg2rad(180/numAngs*(rIdx));
            if isnan(theta)
                theta = 0;
            end
            rho = 1;
            polarplot(pAx,[theta-pi theta],[rho rho],'-','color',cols(rIdx,:),'LineWidth',3)
            
            
end

addExportFigToolbar(gcf)