% Generates plot with normalized tuning curves (for angle of polarization) 
% for TuBu_a and R4m ROIs in the anterior bulb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear x
lineStr = { ...
    'R34H10_Bu';...
    'R34D03_Bu';...
    };
numFlies=[];
for sidx = 1:2
    clear r 
    inclObj = [];
    ROIct = 0;
    eval(['load' lineStr{sidx} ''])
    
    loadLayerMasks([x.MIP])
    loadROIs([x.MIP])

    for oidx = 1:length(x)
        % only examine regular pol mapping exps with 12 angles
        if size(x(oidx).MIP.polAngResp,2)==12 && x(oidx).MIP.polExp > 2
            
            loadROIs(x(oidx).MIP)
            % For each ROI
            for ridx = 1:length(x(oidx).MIP.ROI)
                ROIct = ROIct+1;
                % Get preferred angle 
                ROImask = x(oidx).MIP.ROI(ridx).mask;
                [tun,res] = tuningROI(x(oidx).MIP,x(oidx).MIP.ROI(ridx).mask);
                ROIresp = x(oidx).MIP.polAngResp(ROImask,:);
                
                % Get normalized response profile 
                meanresp = mean(ROIresp,1);
                meanrespnorm = (meanresp - min(meanresp(:)))./(max(meanresp(:)) - min(meanresp(:)));
                % interpolate in fourier domain so we can shift tuning angles to be
                % centered on the tuning peak
                m = interpft(meanrespnorm,360);
                M = circshift(m,30); % first angle in polangresp is 30deg (regardless of motor direction)
                r(ROIct,:) = circshift(M,-round(tun));
            end
            inclObj(end+1) = oidx;
        end
    end
    rs{sidx} = r;
    numFlies(sidx) = length(unique([x(inclObj).Fly]));
end

%% Plot all normalized responses and their mean
%{
printpath = fig6path;
savename = 'TuBu_R4m';

prefix = 'normalized_resp';

figure('color','w')
hold on
for sidx = 1:2
    r = rs{sidx};
    linecolor = ROIcolor(sidx);
    brightcolor = linecolor+(1-linecolor)*0.55;
plot(r','color',brightcolor,'LineWidth',0.5)
pm(sidx) = plot(mean(r,1),'color',linecolor);

end
uistack(pm,'top')
addExportFigToolbar
ylim([-0.2 1.2])
xlim([0 360])
xticks([0:90:360])
pbaspect([1,1,1])
xlabel('angle');
ylabel('norm resp');
setAVPaxes(gca,2)
tightfig(gcf)
suffix = 'traces';
printAVP
%}

%% Plot mean normalized response with a shaded error bar
printpath = fig6path;
savename = 'TuBu_R4m';

prefix = 'normalized_resp';

figure('color','w')
hold on
for sidx = 1:2
    lineprops.col = {[ROIcolor(sidx)]};
    lineprops.width = 1;
    mseb([],mean(rs{sidx},1),std(rs{sidx},[],1),lineprops,1)
end
addExportFigToolbar
ylim([-0.2 1.2])
xlim([0 360])
xticks([0:90:360])
pbaspect([1,1,1])
xlabel('angle');
ylabel('norm resp');
setAVPaxes(gca,2)
tightfig(gcf)
suffix = 'errorbar';
printAVP

%% print info
disp(['[TuBua (blue)]: ' num2str(size(rs{1},1)) ' ROIs, N=' num2str(numFlies(1))])
disp(['[R4m (red)]: ' num2str(size(rs{2},1)) ' ROIs, N=' num2str(numFlies(2))])


