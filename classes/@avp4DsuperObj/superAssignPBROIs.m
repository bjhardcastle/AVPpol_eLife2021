function superAssignPBROIs(objarray,refreshDistMask,manualLimits)
if nargin<2 || isempty(refreshDistMask)
    refreshDistMask = 1;
end
if nargin<3 || isempty(manualLimits)
manualLimits = 1;
else
    % There are 16 e-pg glomeruli (8 per side), but we'll give G1 (medial pair)
    % half the width of all others
    gLimits = linspace(0,1,16); % 16 boundaries to give 15 even widths.
    gLimsL = [ gLimits(1:8) 0.5]; % G1 fills remaining space up to centre
    gLimsR = 1-fliplr(gLimsL); % right side is mirror symmetry
end

for oidx = 1:length(objarray)
    if manualLimits
        gLimits = getPBglomeruliLimits(objarray(oidx).MIP);
        gLimsL = gLimits(1:9);
        gLimsR = gLimits(10:18);
    end

    objarray(oidx).MIP.ROI = [];
    %         figure
    
    % Get mask of PB, where each pixel value is it's normalized position
    % along the PB, from left to right [0:1]
    if refreshDistMask || isempty(objarray(oidx).MIP.PBcurvilinearDistMask)
        getCurvilinearDist(objarray(oidx).MIP,[],[],1)
    end
    dist = objarray(oidx).MIP.PBcurvilinearDistMask;
    
    for n = 1:8 % each glomerulus-pair
        
        % Get boundaries of left glomerulus
        thisLmin = gLimsL(n);
        thisLmax = gLimsL(n+1);
        
        % Get boundaries of right paired glomerulus (shifted one to
        % right: G9L pairs with G2R, not G1R, see Wolff et al 2015)
        
        thisRmin = gLimsR(mod(n,length(gLimsR)-1)+1);
        thisRmax= gLimsR(mod(n+1,length(gLimsR)-1)+1);
        
        
        if thisRmax < thisRmin
            thisRmax = 1;
        end
        
        maskL = zeros(size(dist));
        maskL((dist>=thisLmin & dist<=thisLmax)) = 1;
        
        maskR = zeros(size(dist));
        maskR((dist>=thisRmin & dist<=thisRmax)) = 1;
        
        maskC= zeros(size(dist)); % L + R
        maskC((dist>=thisLmin & dist<=thisLmax) | (dist>=thisRmin & dist<=thisRmax)) = 1;
        
        % make ROI sets: combined(1:8), individual left(9:16),individual right(17:24)
        % note: the index refers to groupings from Fig 14 in
        % Wolff et al 2015, so cells in the two most dorsal EB wedges are
        % grouped, which is shifted-by one from the typical locust grouping
        % of L9/R1, L8,R2, etc.
        %
        % ROIindices: [9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
        % glomeruli: [L9,8,7,6,5,4,3,2,1,R2,3,4,5,6,7,8,9,1]
        
        objarray(oidx).MIP.ROI(n).mask = logical(maskC);
        objarray(oidx).MIP.ROI(n+8).mask = logical(maskL);
        objarray(oidx).MIP.ROI(n+16).mask = logical(maskR);
        
        for m = 1:3
            rIdx = (m-1)*8 + n;
            
            
            objarray(oidx).MIP.ROI(rIdx).color = ROIcolor(n);
            
            % Make an outline just so it displays in 'play' or showROIs...
            % shouldn't be used since all but one mask consists of a pair of
            % separated regions
            % Trace outline of filled mask areas
            B = bwboundaries(objarray(oidx).MIP.ROI(rIdx).mask);
            if ~isempty(B)
                roix = [];
                roiy = [];
                for b = 1:length(B)
                    roix = [roix B{b}(:,2)'];
                    roiy = [roiy B{b}(:,1)'];
                end
                % Downsample the outline by skipping every nth point
                skipPts = 3;
                roiPos = [ roix([1:skipPts:end])' roiy([1:skipPts:end])'];
            end
            
            objarray(oidx).MIP.ROI(rIdx).position = roiPos;
            objarray(oidx).MIP.ROI(rIdx).response = scanROI(objarray(oidx).MIP,objarray(oidx).MIP.ROI(rIdx).mask);
            objarray(oidx).MIP.UseFixedResp = 1;
            
            % Also get the tuning preference of the ROI
            polExp = objarray(oidx).MIP.polExp;
            if ~isempty(polExp)
                ROIresp = mean(objarray(oidx).MIP.polAngResp(objarray(oidx).MIP.ROI(rIdx).mask,:),1);
                polAngsOrig = [objarray(oidx).MIP.pSet(polExp).polAngleStep:objarray(oidx).MIP.pSet(polExp).polAngleStep:360];
                [~,axtheta] = circ_rangle(circ_axial(deg2rad(repmat(polAngsOrig(1:length(polAngsOrig)),size(ROIresp,1),1)),2), ROIresp, deg2rad(mode(diff(polAngsOrig))), 2);
                theta = axtheta/2;
                Ang =0.5.*(wrapTo360((2.*(rad2deg(theta)))));
                objarray(oidx).MIP.ROI(rIdx).angle = Ang;
            end
        end
        
        % subplot(1,8,n)
        % imshowpair(objarray(oidx).MIP.ROI(n).mask,objarray(oidx).MIP.InactivityFrame,'falsecolor')
        
    end
    saveROIs(objarray(oidx).MIP)
end
