function plotVectormap(obj)
colormapPlot = 0;
vectorPlot = 1;

% Ensure correct type of object:
assert(strcmp( class(obj), 'avp4DmaxObj'), 'Object must be of type ''avpSubsetObj''')

if isempty(obj.Frames)
    obj.getFrames;
end
if isempty(obj.TrialPatNum)
   getParameters(obj) 
end
% Make a NaN vector the length of the number of frames for this subset
            polVec = nan(1,size(obj.Frames,3));
            
             % Populate vector with polarizer angle at each frame
             AngleTrials = find(obj.TrialSeqNum == 4);
            for tidx = 1:length(AngleTrials)
                polVec( obj.TrialStartFrame(AngleTrials(tidx)) :  obj.TrialEndFrame(AngleTrials(tidx)) ) = obj.TrialPatNum(AngleTrials(tidx));
            end
            
            % Make additional vectors for each pref+antifpref angle
            % combination within the experiment. This will aid 
            % finding max response quickly:
            angs = unique(obj.TrialPatNum(obj.TrialSeqNum == 4));
            respVec = nan(0.5*length(angs), size(obj.Frames,3));
            for aidx = 1:0.5*length(angs)
                respVec( aidx, (polVec == angs(aidx)) ) = 1;
                respVec( aidx, (polVec == mod(angs(aidx) + 180, 360) ) ) = 1;
            end
            
            % Get 2D frame indices:
            [row2D,col2D] = find(obj.Frames(:,:,1)>-9999);
            PS1D = nan(1, length(row2D));
            PD1D = nan(1, length(row2D));
            
            for pixIdx = 1:length(row2D)
                
                % Extract pixel resp vector:
                respSeries = squeeze(obj.Frames(row2D(pixIdx),col2D(pixIdx),:));
                
                respMean = [];
                for rIdx = 1:size(respVec,1)
                    respMean(rIdx) = nanmean(respSeries.*respVec(rIdx,:)');
                end
                if all(respMean == 0 ) || nanmax(respMean) == 0
                    PS1D(pixIdx) = 0;
                    PD1D(pixIdx) = NaN;
                else
                    [prefResp,rMaxIdx] = max(respMean);
                    
                    % Get the pol direction preference:
                    prefAng = angs(rMaxIdx);
                    
                    % Extract pixel's mean value during:
                    
                    % frames whose angle matches the angle that occured at the max resp
                    prefResp1 = mean( respSeries( polVec == prefAng ) );
                    
                    % frames whose angle matches the angle that occured
                    % at the max resp +180
                    prefResp2 = mean( respSeries( polVec == mod(prefAng+180, 360 ) ));
                    
                    % frames whose angle matches the angle that occured at the max resp
                    % + 90 deg
                    nullResp1 = mean( respSeries( polVec == mod(prefAng+90, 360 ) ) );
                    % frames whose angle matches the angle that occured at the max resp
                    % - 90 deg
                    nullResp2 = mean( respSeries( polVec == mod(prefAng-90, 360 ) ) );
                    
                    
                    % Calculate polarization selectivity:
                    PS1D(pixIdx) = (  max([prefResp1 prefResp2]) -  mean([nullResp1 nullResp2]) );% / max([prefResp1 prefResp2]) ;
                    %PS1D(pixIdx) = 1 - ( nanmean([nullResp1 nullResp2]) / nanmax([prefResp1 prefResp2]) );
                    
                    % And store the pol direction preference:
                    prefAngs = [prefAng mod(prefAng+180, 360)];
                    [~,prIdx] = max([ prefResp1 prefResp2 ]);
                    PD1D(pixIdx) = prefAngs(prIdx);
                    
                end
                
            end
            
            % Create pol selectivity image:
            PS1Dpos = PS1D;
            PS1Dpos(PS1D<0) = 0;
            PS1Dpos(isnan(PS1D))  = 0;
            PS2D = reshape( PS1Dpos , size(obj.Frames,1) , size(obj.Frames,2) );
            PSimg = PS2D;
            
            %figure,imshow(PSimg,[])
            
            % Create pol Direction tuning image, with selectivity image as a filter:
            PD1DnanReplaced = PD1D;
            PD1DnanReplaced(isnan(PD1D)) = 360+mean(diff(angs)); % Pixels with NaN direction preference (outside trial) are set to grey
            %                 PD1DnanReplaced( (PS1Dpos) <= 0 ) = 360+mean(diff(angs)); % Pixels with low or negative PSelectivity are also set to grey
            PD1DnanReplaced( (PS1Dpos) <= nanmedian(PS1Dpos(:)) ) = 360+mean(diff(angs)); % Pixels with low or negative PSelectivity are also set to grey
            PD2D = reshape( PD1DnanReplaced , size(obj.Frames,1) , size(obj.Frames,2) );
            
%             obj.PDFrame = PD2D;
%             obj.PSFrame = PS2D;
%             

        
        if colormapPlot
            
            % Create colormap for 0:180deg
            c = hsv(360/mean(diff(angs)));
            % Cat with itself to give color for 0 = color for 180,
            % 90=270 etc
            c(end+1:end+length(c),:) = c;
            % And add grey, for excluded pixels
            c(end+1,:) = [0.2 0.2 0.2];
            
            % Apply colormap (colors are in discrete 10deg steps in pol angle: divide angles by 10)
            PDcol =  ind2rgb((PD2D./mean(diff(angs)) + 360/mean(diff(angs))),c);
            PDimg = PDcol;
            PDimg = [];
            % Remove noise from image on RGB channels seperately:
            filtSize = [3 3];
            %hh = fspecial('average', filtSize);
            %PDimg(:,:,1) = imfilter(PDcol(:,:,1) ,hh);
            %PDimg(:,:,2) = imfilter(PDcol(:,:,2), hh);
            %PDimg(:,:,3) = imfilter(PDcol(:,:,3), hh);
            PDimg(:,:,1) = medfilt2(PDcol(:,:,1) ,filtSize);
            PDimg(:,:,2) = medfilt2(PDcol(:,:,2), filtSize);
            PDimg(:,:,3) = medfilt2(PDcol(:,:,3), filtSize);
            
            
            % Any pixel RGB vals which fall between the discrete colors
            % in the colormap must have a mix of direction tuning
            % values: to smooth the image and discard further noise,
            % only keep pixels who have commonly tuned neighbours
            for pixIdx = 1:length(row2D)
                if ~(ismember( squeeze(PDimg(row2D(pixIdx),col2D(pixIdx),:))', c(1:end-1,:), 'rows'))
                    PDimg(row2D(pixIdx),col2D(pixIdx),:) = c(end,:);
                end
            end
            
            % overlay = imfuse(obj.AverageFrame, PDimg ,'method','blend','scaling','joint');
            
            figure('color',[0.3 0.3 0.3],'Position',get(0, 'Screensize'));
            subplot(1,4,3)
            imshow( PDimg ,[])
            title('Pol Angle Preference')
            subplot(1,4,2)
            imshow(PSimg,[])
            title('Pol Selectivity')
            subplot(1,4,1)
            imshow(obj.AverageFrame,[])
            title(obj.File,'interpreter','none');
            
            r = linspace(0,1,10);
            theta = linspace(2*pi, 0, 100);
            [rg, thg] = meshgrid(r,theta);
            [x,y] = pol2cart(thg,rg);
            subplot(1,4,4)
            % Reverse direction, so +ve angle from 0 is cw
            pcolor(-x,y,thg);
            % Shift colorwheel 90deg ccw to put 0deg at top
            shiftedcmap = circshift(c(1:end-1,:),(720/mean(diff(angs)))/4);
            colormap( shiftedcmap );
            shading flat;
            axis equal tight
            axis off
            ax = gca;
            ax.Position(3) = ax.Position(3).*0.2;
            title('0')
            
%             export_fig(fullfile(newstorePath, [num2str(objIdx),'cmap','Exp',num2str(subIdx),'.png']))
            
%             close gcf
        end
        
        if vectorPlot
            % PDimg to PEllipse
% PD2D = obj.PDFrame;
% PS2D = obj.PSFrame;
pd = PD2D;
ps = PS2D;

if isempty(obj.AverageFrame)
    obj.AverageFrame = mean(obj.Frames,3);
end


downRate = floor(0.0075*size(pd,2)); % Sampling will be a square of side 2*downRate + 1

% % % Downsample, then upsample, to give a coarser representation of data
% % downRate = 5;
psDown=imresize(ps,1/(downRate*2+1),'method','nearest');
psUp=imresize(psDown,(downRate*2+1),'method','nearest');
% % pdDown=imresize(pd,1/downRate,'method','nearest');
% % pdUp=imresize(pdDown,downRate,'method','nearest');


% Upsample the 'base' intensity image onto which vectors will be drawn:
upRate = 4;

avgF = imresize(obj.AverageFrame,upRate); % Expand avg frame 4x
gridp = 0.5*(downRate*2+1)*upRate; % 1 pix (data point) in pdDown  =  downRate pix in original (pd)  =   upRate*downRate pix in upsampled img (avgF)

xp = [gridp:2*gridp:size(avgF,2)];
yp = [gridp:2*gridp:size(avgF,1)];

xd = [downRate+1:2*downRate+1:size(pd,2)-(downRate+1)];
yd = [downRate+1:2*downRate+1:size(pd,1)-(downRate+1)];

% [xg, yg] = meshgrid(xp,yp);
figure, imshow( avgF,[0 1.5*max(obj.AverageFrame(:))]);
colormap(gca,copper)
hold on
vecScale = gridp*1.2/nanmax(psDown(:));
lineCol = [0.8 0.8 1];
for yn = 1:length(yd)
    for xn = 1:length(xd)
        
        % Get the nxn square of pol direction data points
        pdSample = pd( [yd(yn)-downRate : yd(yn)+downRate] , [xd(xn) - downRate : xd(xn)+downRate] );
        pdSample = pdSample(:);
        pdSample(pdSample==360 + mean(diff(angs))) = NaN;
        
        % theta = mod(rad2deg(circ_median(deg2rad(pdSample))),360);
        
        % FInd the pixel with maximum polarization selectivity
        psSample = ps( [yd(yn)-downRate : yd(yn)+downRate] , [xd(xn) - downRate : xd(xn)+downRate] );
        psSample = psSample(:);
        [maxPS, maxPSIdx] = nanmax(psSample);
        
        % Use that pixel's direction for the vector angle:
        if maxPS ~= 0 && ~isnan(maxPS)
            theta =  pdSample( maxPSIdx );
        else
            theta = nan;
        end
        % And find the mean selectivity of the sample square for the
        % magnitude. Vector length is then a mesaure of selectivity for
        % the area, and direction is the most representative angle. Mean
        % angle would be meaningless when there's lots of variation,
        % particularly at low angle-step sampling
        
        vecLength = psDown(yn,xn)*vecScale;
        
        % Only draw the line if it has non-zero length:
        if ~isnan(theta) && ~isnan(vecLength) && vecLength ~= 0
            % Create one half of line, from origin in preferred direction:
            x1 = xp(xn);
            y1 = yp(yn);
            
            x2 = x1 + vecLength*cosd(theta);
            y2 = y1 + vecLength*sind(theta);
            
            line([x1 x2], [y1 y2],'LineWidth',1,'Color',lineCol)
            
            % Create other half, from origin in anti-preferred direction:
            x1 = xp(xn);
            y1 = yp(yn);
            
            x2 = x1 + vecLength*cosd(mod(theta+180,360));
            y2 = y1 + vecLength*sind(mod(theta+180,360));
            
            line([x1 x2], [y1 y2],'LineWidth',1,'Color',lineCol)
            
        end
        
        
    end
end

%             export_fig(fullfile(newstorePath, [num2str(objIdx),'vecmap','Exp',num2str(subIdx),'MAX.png']),'-m2')
%             close gcf
        end
        
        
    end
    
   