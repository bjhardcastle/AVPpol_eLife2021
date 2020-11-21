% Get ROI tuning/selectivity data for ROIs in:
%   TuBu_a   (anterior bulb)
%   R4m      (anterior bulb)
%   E-PG     (protocerebral bridge)
% for scatter plots, polarotopy stats, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathsAVP

lineStr = {...
    'R34H10_Bu'; ...
    'R34D03_Bu'; ...
%     'SS00096_PB'; ...
    };
% Note:  This E-PG struct has been superceded by another scatter plot
% Here, the PSI values used were calculated from first 2 cycles of each pol
% mapping experiment. We later understood that these weren't reliable
% (indicated by the low PSI) and sought single cycles where all glomeruli
% contained reliable responses, and then plot their tuning in a scatter
% plot. The results turned out to be qualitatively similar (no correlation
% between position and tuning) but given the inconsistencies in E-PG
% responses, the later version is the more correct version (see
% 'make_various_EPG_panels.m')

for Lidx = 1:length(lineStr)

    % load object array
    eval(['load' lineStr{Lidx}])           
 
    lastRec = 0; % if zero, first matching recording is found for each fly
    inclOnly = 0; % if one, only use x(incl) (as specified in obj load function)
    
    loadROIs([x.MIP])
    flies = unique([x.Fly]);
    ROIstruct = struct();
    cidx = 0; %counter
    for fidx = 1:length(flies)
        if lastRec
            % For each fly (going in rev order)
            thisFly = fliplr(find([x.Fly]==flies(fidx)));
        else
            % For each fly (in chronological order)
            thisFly = (find([x.Fly]==flies(fidx)));
        end
        clear tidxL tidxR tidxPB
        tidxL = [];
        tidxR = [];
        tidxPB = [];
        for tidx = 1:length(thisFly)
            if isempty(tidxL) && contains(x(thisFly(tidx)).Name,' L ') && ~isempty(x(thisFly(tidx)).MIP.polExp) && x(thisFly(tidx)).MIP.polExp >2
                tidxL = thisFly(tidx);
            end
            if isempty(tidxR) && ~contains(x(thisFly(tidx)).Name,' L ') && ~isempty(x(thisFly(tidx)).MIP.polExp) && x(thisFly(tidx)).MIP.polExp>2
                tidxR = thisFly(tidx);
            end
            if isempty(tidxPB) && ~isempty(x(thisFly(tidx)).MIP.polExp) && x(thisFly(tidx)).MIP.polExp>2
                if ~inclOnly || ismember(thisFly(tidx),incl)
                    tidxPB = thisFly(tidx);
                end
            end
        end
        
        switch x(tidx).Line
            case 'SS00096-gal4'
                if ~isempty(tidxPB)
                    cidx = cidx +1;
                else
                    continue
                end
            otherwise
                if ~isempty([tidxL tidxR])
                    cidx = cidx+1;
                else
                    continue
                end
        end
        
        for RLidx = 1:2
            switch x(1).Line
                case 'SS00096-gal4'
                    if isempty(tidxPB)
                        warning('something wrong with tidx')
                        continue
                    else
                        tidx = tidxPB;
                    end
                    clear thisROIcenter % only necessary when running as a script, after looking at non-EPG lines with 2 center coords
                    
                    % find xy center of layermask (will be the same for L and R)
                    m = regionprops(x(tidx).MIP.layerMask.mask,'centroid');
                    thisMaskCenter = m.Centroid;
                    
                    % For each roi
                    for ridx = [1:8]+8*RLidx % look at individual ROIs on left and right side and plot separately
                        
                        % find xy center)
                        c = regionprops(x(tidx).MIP.ROI(ridx).mask,'centroid');
                        if isempty(x(tidx).MIP.PBcurvilinearDistMask)
                            getCurvilinearDist(x(tidx).MIP,[],[],1);
                        end
                        thisROIcenter = x(tidx).MIP.PBcurvilinearDistMask(round(c.Centroid(2)),round(c.Centroid(1))); % 1d coords
                        
                        % get tuning, resultant vector length
                        [thisTuning,thisRlength] = tuningROI(x(tidx).MIP,x(tidx).MIP.ROI(ridx).mask,1);
                        
                        % get psi
                        thisPSI = mean(x(tidx).MIP.polSelImg(x(tidx).MIP.ROI(ridx).mask));
                        
                        % We'll make three plots:
                        % left ROIs [0:0.5]
                        % right ROIs [0.5:1]
                        % joint l/r ROIs [0:1]
                        
                        % Add all ROIs to plot
                        switch RLidx
                            case 1 % left
                                %plot(ax(1),thisROIcenter,thisTuning,'.','color',cols(round(thisPSI*255+1),:));
                                %
                                ROIstruct(cidx).Left.PSI(ridx-8) = thisPSI;
                                ROIstruct(cidx).Left.tuning(ridx-8) = thisTuning;
                                ROIstruct(cidx).Left.posPB(ridx-8) = thisROIcenter; % position [0:1]
                                ROIstruct(cidx).Left.posPBnorm(ridx-8) = (thisROIcenter-(RLidx-1)*0.5)/0.5; % position normalized mod([0:2],1)
                                ROIstruct(cidx).Left.Fly =  x(tidx).Fly;
                                ROIstruct(cidx).Left.DateTime = [x(tidx).DateStr x(tidx).TimeStr];
                            case 2 % right
                                %plot(ax(2),thisROIcenter,thisTuning,'.','color',cols(round(thisPSI*255+1),:));
                                
                                ROIstruct(cidx).Right.PSI(ridx-8) = thisPSI;
                                ROIstruct(cidx).Right.tuning(ridx-8) = thisTuning;
                                ROIstruct(cidx).Right.posPB(ridx-8) = thisROIcenter; % position [0:1]
                                ROIstruct(cidx).Right.posPBnorm(ridx-8) = (thisROIcenter-(RLidx-1)*0.5)/0.5; % position normalized mod([0:2],1)
                                ROIstruct(cidx).Right.Fly =  x(tidx).Fly;
                                ROIstruct(cidx).Right.DateTime = [x(tidx).DateStr x(tidx).TimeStr];
                                
                                % also add joint rois, which are in positions
                                % 1:8 (same phase as ridx)
                                glomPos = mod(ridx-1,8)+1;
                                jointTuning = tuningROI(x(tidx).MIP,x(tidx).MIP.ROI(glomPos).mask,1);
                                %plot(ax(3),glomPos,jointTuning,'.','color',cols(round(thisPSI*255+1),:));
                                
                                ROIstruct(cidx).Joint.PSI(glomPos) = mean(x(tidx).MIP.polSelImg(x(tidx).MIP.ROI(glomPos).mask));
                                ROIstruct(cidx).Joint.tuning(glomPos) = jointTuning;
                                ROIstruct(cidx).Joint.Fly =  x(tidx).Fly;
                                ROIstruct(cidx).Joint.DateTime = [x(tidx).DateStr x(tidx).TimeStr];
                                %                             if isnan( autocorrthresh( tidx, glomPos) ) || isnan( autocorrthresh( tidx, glomPos+8))
                                %                                 ROIstruct(cidx).Joint.PSI(glomPos) = 0;
                                %                             end
                                %                             if ROIstruct(cidx).Joint.PSI(glomPos) <mean(x(tidx).MIP.polSelImg(x(tidx).MIP.bkgMask))
                                %                                 ROIstruct(cidx).Joint.PSI(glomPos) = 0;
                                %                             end
                        end
                        
                    end
                otherwise
                    % For left / right separately
                    if RLidx == 1 && ~isempty(tidxL)
                        tidx = tidxL;
                    elseif RLidx == 2 && ~isempty(tidxR)
                        tidx = tidxR;
                    else
                        continue
                    end
                    % find xy center of layermask (will be the same for L and R)
                    m = regionprops(x(tidx).MIP.layerMask.mask,'centroid');
                    thisMaskCenter = m.Centroid;
                    
                    % For each roi
                    for ridx = 1:length(x(tidx).MIP.ROI)
                        % find xy center
                        c = regionprops(x(tidx).MIP.ROI(ridx).mask,'centroid');
                        thisROIcenter = c.Centroid;
                        
                        % get tuning, resultant vector length
                        [thisTuning,thisRlength] = tuningROI(x(tidx).MIP,x(tidx).MIP.ROI(ridx).mask,1);
                        
                        % get psi
                        thisPSI = mean(x(tidx).MIP.polSelImg(x(tidx).MIP.ROI(ridx).mask));
                        
                        % find angle of each roi
                        thisAng = atan2d(thisROIcenter(1)-thisMaskCenter(1), thisROIcenter(2)-thisMaskCenter(2));
                        
                        % find vert/horiz limits of mask for normalizing position of x y roi centers
                        [i,j] = ind2sub(size(x(tidx).MIP.layerMask.mask),find(x(tidx).MIP.layerMask.mask));
                        horizmin = min(i);
                        horizwidth = max(i)-min(i); % horiz in brain, not image
                        vertmin = min(j);
                        vertwidth = max(j)-min(j);
                        
                        % find x/y center of mass of ROI, in terms of normalized mask
                        % width/height [0:1]. Instead of cutting centers outside the
                        % mask, clip them to the limits (should be close anyway, and
                        % this way outliers won't skew the position of all others)
                        thisYpos = (thisROIcenter(1) - vertmin)/vertwidth;
                        thisXpos = (thisROIcenter(2) - horizmin)/horizwidth;
                        thisYpos(thisYpos > 1) = 1;
                        thisYpos(thisYpos < 0) = 0;
                        thisXpos(thisXpos > 1) = 1;
                        thisXpos(thisXpos < 0) = 0;
                        
                        % Add all ROIs to plot
                        switch RLidx
                            case 1 % left
                                % circ
                                %plot(ax(1),thisAng,thisTuning,'.','color',cols(round(thisPSI*255+1),:));
                                % vert
                                %plot(ax(3),thisYpos,thisTuning,'.','color',cols(round(thisPSI*255+1),:));
                                % horiz
                                %plot(ax(5),thisXpos,thisTuning,'.','color',cols(round(thisPSI*255+1),:));
                                
                                ROIstruct.Left(cidx).PSI(ridx) = thisPSI;
                                ROIstruct.Left(cidx).tuning(ridx) = thisTuning;
                                ROIstruct.Left(cidx).posX(ridx) = thisXpos;
                                ROIstruct.Left(cidx).posY(ridx) = thisYpos;
                                ROIstruct.Left(cidx).posAng(ridx) = thisAng;
                                ROIstruct.Left(cidx).Fly = x(tidx).Fly;
                                ROIstruct.Left(cidx).DateTime = [x(tidx).DateStr x(tidx).TimeStr];
                                
                            case 2 % right
                                
                                
                                ROIstruct.Right(cidx).PSI(ridx) = thisPSI;
                                ROIstruct.Right(cidx).tuning(ridx) = thisTuning;
                                ROIstruct.Right(cidx).posX(ridx) = thisXpos;
                                ROIstruct.Right(cidx).posY(ridx) = thisYpos;
                                ROIstruct.Right(cidx).posAng(ridx) = thisAng;
                                ROIstruct.Right(cidx).Fly = x(tidx).Fly;
                                ROIstruct.Right(cidx).DateTime = [x(tidx).DateStr x(tidx).TimeStr];
                                
                        end
                    end
                    
            end
        end
    end
    
    switch x(1).Line
        case 'SS00096-gal4'
            %ax(1) corr: circ pos, circ tuning
            %ax(2) corr: circ pos, circ tuning
            % numFlies are always 16 per fly
            %ax(3) corr: circ pos, circ tuning (if at all)
            % numFlies are always 8 per fly
            ROIstructLeft = [ROIstruct.Left];
            ROIstructRight = [ROIstruct.Right];
            ROIstructPooled = [ROIstruct.Joint];
            
            % number of flies is the same for all
            ROIstats = struct;
            ROIstats.All.numFlies = length(ROIstructLeft(cellfun(@length,{ROIstructLeft.PSI})~=0));
              
            save(...
                ROIstruct_path(lineStr{Lidx}), ...
                'ROIstructLeft',...
                'ROIstructRight',...
                'ROIstructJoint',...
                'ROIstats'...
                );
            
        otherwise
            %ax(1:2) corr: circ pos, circ tuning
            %ax(3:4) corr: linear pos, circ tuning
            %ax(5:6) corr: linear pos, circ tuning
            ROIstructLeft = [ROIstruct.Left];
            ROIstructRight = [ROIstruct.Right];
            ROIstats = struct;

            % store some info on ROIs in each struct
            ROIstats.Left.numFlies = length(ROIstructLeft(cellfun(@length,{ROIstructLeft.PSI})~=0));
            ROIstats.Right.numFlies = length(ROIstructRight(cellfun(@length,{ROIstructRight.PSI})~=0));
            ROIstats.Left.numROIsPerFly = cellfun(@length,{ROIstructLeft.PSI});
            ROIstats.Right.numROIsPerFly = cellfun(@length,{ROIstructRight.PSI});
            ROIstats.Left.numROIsNonZero = ROIstats.Left.numROIsPerFly(ROIstats.Left.numROIsPerFly>0);
            ROIstats.Right.numROIsNonZero = ROIstats.Right.numROIsPerFly(ROIstats.Right.numROIsPerFly>0);
            ROIstats.Left.numROIsPerFlyAvg = mean(ROIstats.Left.numROIsNonZero);
            ROIstats.Right.numROIsPerFlyAvg = mean(ROIstats.Right.numROIsNonZero);
            ROIstats.Left.numROIsPerFlyStd = std(ROIstats.Left.numROIsNonZero);
            ROIstats.Right.numROIsPerFlyStd = std(ROIstats.Right.numROIsNonZero);
            ROIstats.Left.numROIsTotal = sum(ROIstats.Left.numROIsPerFly);
            ROIstats.Right.numROIsTotal = sum(ROIstats.Right.numROIsPerFly);

            % assemble pooled structure from right data and then add mirrored left data 
            ROIstructPooled = [ROIstruct.Right];
            for jidx = 1:length(ROIstructLeft)                 
                if ~isempty(ROIstructLeft(jidx).PSI)
                    leftData = ROIstructLeft(jidx);
                    leftData.posX = 1-leftData.posX;
                    leftData.posY = 1-leftData.posY;
                    leftData.posAng = wrapTo360(leftData.posAng)-180; % range [-180:180]
                
                    % if data from this fly already exists, pool both sides
                [~,flyIdx] = ismember(ROIstructLeft(jidx).Fly,[ROIstructPooled.Fly]);
                if flyIdx > 0
                    ROIstructPooled(flyIdx).PSI = [ROIstructPooled(flyIdx).PSI leftData.PSI];
                    ROIstructPooled(flyIdx).tuning = [ROIstructPooled(flyIdx).tuning leftData.tuning];
                    ROIstructPooled(flyIdx).posX = [ROIstructPooled(flyIdx).posX leftData.posX];
                    ROIstructPooled(flyIdx).posY = [ROIstructPooled(flyIdx).posY leftData.posY];
                    ROIstructPooled(flyIdx).posAng = [ROIstructPooled(flyIdx).posAng leftData.posAng];

                else % otherwise append data
                    ROIstructPooled(end+1) = leftData;
                end
                end
            end
            ROIstats.Pooled.numFlies = length(ROIstructPooled(cellfun(@length,{ROIstructPooled.PSI})~=0));
            ROIstats.Pooled.numROIsPerFly = cellfun(@length,{ROIstructPooled.PSI});
            ROIstats.Pooled.numROIsNonZero = ROIstats.Pooled.numROIsPerFly(ROIstats.Pooled.numROIsPerFly>0);
            ROIstats.Pooled.numROIsPerFlyAvg = mean(ROIstats.Pooled.numROIsNonZero);
            ROIstats.Pooled.numROIsPerFlyStd = std(ROIstats.Pooled.numROIsNonZero);
            ROIstats.Pooled.numROIsTotal = sum(ROIstats.Pooled.numROIsNonZero);
            
            save(...
                ROIstruct_path(lineStr{Lidx}), ...
                'ROIstructLeft',...
                'ROIstructRight',...
                'ROIstructPooled',...
                'ROIstats'...
                );
            
    end
    
end