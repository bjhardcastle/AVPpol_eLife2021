function superPolCtrl(objarray)
% Plot the distribution of polarization selectivity for pixels in control
% experiments with no polarizer and experiments with polarizer, for
% bkground areas and cells. Then run stats to see if there's an
% effect of having the polarizer

intensityFraction = 0.10;  % Proportion of pixels in mask that will be used to form bkgrnd area (lowest) and cells (highest)
% backgroundOutsideMask = 1; % default is background pixels within mask area
binwidth = 1/50;
pooledPolSelValues(1:10) = struct('background',[],'cells',[],'flies',[]);
meanPolSelValues(1:10) = struct('background',nan,'cells',nan,'flies',nan);

% Make sure there are some recordings with and without polarizer
if ~any([objarray.Exps]==2) 
    disp('No control experiments (exp 2) exist')
    
elseif ~any([objarray.Exps]==2 | [objarray.Exps]==4 | [objarray.Exps]==10) 
    disp('No pol experiments (exp 2 or 4) exist')
    return
end

for oidx = 1:length(objarray)
    % Find which experiment to use for this recording: exp4 and exp2 are
    % mutually exclusive
    if any([objarray(oidx).Exps]==2)
        expNum = 2;
    elseif any([objarray(oidx).Exps]==4)
        expNum = 4;
    elseif any([objarray(oidx).Exps]==10)
        expNum = 10;
    else
        continue
    end
    
    if length(objarray(oidx).Layers) == 1 && isempty(objarray(oidx).MIP)
        objarray(oidx).MIP = objarray(oidx).Layers;
    end
%     getPolMaps(objarray(oidx).MIP)

    if isempty(objarray(oidx).MIP.layerMask)
        loadLayerMasks(objarray(oidx).MIP)
        if isempty(objarray(oidx).MIP.layerMask)
            maskLayers(objarray(oidx))
            [~]=input('Waiting for layer mask - press any key to continue');
        end
    end
    
    % Get original MIP frames for looking at average brightness values
    if isempty(objarray(oidx).MIP.ActivityFrame)
        if objarray(oidx).MIP.UseMSP
            objarray(oidx).MIP.UseMSP = 0;
            undoFlag = 1;
        else
            undoFlag = 0;
        end
        objarray(oidx).MIP.Unattended = 1;
        getFrames(objarray(oidx).MIP)
        objarray(oidx).MIP.Frames = [];
        if undoFlag
            objarray(oidx).MIP.UseMSP = 1;
        end
        objarray(oidx).MIP.Frames = [];
    end

    %{
    maskedMaxFrameVector = objarray(oidx).MIP.ActivityFrame(objarray(oidx).MIP.layerMask.mask);        
%     maskedMaxFrameVector = objarray(oidx).MIP.polSelImg(objarray(oidx).MIP.layerMask.mask);        

    [~,frameIdx] = sort(maskedMaxFrameVector,'descend');

    topQuart = frameIdx(1:floor(length(frameIdx)*intensityFraction));
    bottomQuart = frameIdx(ceil(length(frameIdx)-length(frameIdx).*intensityFraction)+1:end);

    maskedPolSelFrameVector = objarray(oidx).MIP.polSelImg(objarray(oidx).MIP.layerMask.mask);

    pooledPolSelValues(expNum).background = [pooledPolSelValues(expNum).background;maskedPolSelFrameVector(bottomQuart)];
    pooledPolSelValues(expNum).cells = [pooledPolSelValues(expNum).cells;maskedPolSelFrameVector(topQuart)];
    pooledPolSelValues(expNum).flies = [pooledPolSelValues(expNum).flies;objarray(oidx).Fly];
        
    meanPolSelValues(expNum).background(end+1) = median(maskedPolSelFrameVector(bottomQuart));
    meanPolSelValues(expNum).cells(end+1) = median(maskedPolSelFrameVector(topQuart));
    meanPolSelValues(expNum).flies(end+1) = objarray(oidx).Fly;
    
    %}
 %{ 
      maskedPolSelFrameVectorCells = getPolSelMaskedVector( objarray(oidx).MIP, 'top');
    maskedPolSelFrameVectorBkg = getPolSelMaskedVector( objarray(oidx).MIP, 'bottom');

    pooledPolSelValues(expNum).background = [pooledPolSelValues(expNum).background;maskedPolSelFrameVectorBkg];
    pooledPolSelValues(expNum).cells = [pooledPolSelValues(expNum).cells;maskedPolSelFrameVectorCells];
    pooledPolSelValues(expNum).flies = [pooledPolSelValues(expNum).flies;objarray(oidx).Fly];
        
    meanPolSelValues(expNum).background(end+1) = median(maskedPolSelFrameVectorBkg);
    meanPolSelValues(expNum).cells(end+1) = median(maskedPolSelFrameVectorCells);
    meanPolSelValues(expNum).flies(end+1) = objarray(oidx).Fly;
   %}
 
    maskedPolSelFrameVectorCells = objarray(oidx).MIP.polSelImg(objarray(oidx).MIP.cellMask);
    maskedPolSelFrameVectorBkg = objarray(oidx).MIP.polSelImg(objarray(oidx).MIP.bkgMask);

    pooledPolSelValues(expNum).background = [pooledPolSelValues(expNum).background;maskedPolSelFrameVectorBkg];
    pooledPolSelValues(expNum).cells = [pooledPolSelValues(expNum).cells;maskedPolSelFrameVectorCells];
    pooledPolSelValues(expNum).flies = [pooledPolSelValues(expNum).flies;objarray(oidx).Fly];
        
    meanPolSelValues(expNum).background(end+1) = median(maskedPolSelFrameVectorBkg);
    meanPolSelValues(expNum).cells(end+1) = median(maskedPolSelFrameVectorCells);
    meanPolSelValues(expNum).flies(end+1) = objarray(oidx).Fly;
   
end

%%
f = figure('color','w');

titleStr{1} = 'no polarizer';
titleStr{2} = 'with polarizer';
titleStr{3} = 'background';
titleStr{4} = 'cells';

cols(1,:) = ROIcolor(16);
cols(2,:) = [0.7 0.7 0.7];
cols(3,:) = ROIcolor(1);
cols(4,:) = ROIcolor(2);


for a = 1:4
    ax(a) = subplot(2,2,a);
    hold(ax(a),'on')
    
    if a < 3 %  exp2 (bkgrnd vs cell)  |  exp4 (bkgrnd vs cell)
        aExp = a*2;
        
        vals1 = pooledPolSelValues(aExp).background;
        vals2 = pooledPolSelValues(aExp).cells;
        if  isempty(vals1) || isempty(vals2)
            cla(ax(a))
             axis(ax(a),'off')
            continue
        else
        histogram(ax(a),vals1,'binwidth',binwidth,'normalization','probability','FaceColor',cols((a-1)*2+1,:))
        histogram(ax(a),vals2,'binwidth',binwidth,'normalization','probability','FaceColor',cols((a-1)*2 + 2,:))
        
        % equal group sizes since we're testing group a vs b in the same
        % recording
%         k = kruskalwallis([vals1,vals2],titleStr(3:4),'off');
          k = ranksum(vals1,vals2);
        p(a) = k;
        
        t = title(ax(a),sprintf('%s, mean=%1.3f N=%d vs. mean=%1.3f N=%d, p=%1.6f',titleStr{a},mean(vals1),length(unique(pooledPolSelValues(2).flies)),mean(vals2),length(unique(pooledPolSelValues(4).flies)),p(a) ));
        L=legend(titleStr{3:4});
        end
        
    else   %  bkgrnd (exp2 vs exp4)  |  cell (exp2 vs exp4)        
            
        vals1 = pooledPolSelValues(2).(titleStr{a});
        vals2 = pooledPolSelValues(4).(titleStr{a});
         if isempty(vals1) || isempty(vals2)
               cla(ax(a))
               axis(ax(a),'off')
            continue
        else
        histogram(ax(a),vals1,'binwidth',binwidth,'normalization','probability','FaceColor',cols(a-2,:))
        histogram(ax(a),vals2,'binwidth',binwidth,'normalization','probability','FaceColor',cols(a,:))
        
        % unequal group sizes now we're comparing across different
        % recordings (to add/remove polarizer between exp2 and exp4 req. stopping recording)
        group1 = repmat(titleStr(1),length(vals1),1);
        group2 = repmat(titleStr(2),length(vals2),1);
        
%         k  = kruskalwallis([vals1;vals2],[group1;group2],'off');
        k = ranksum(vals1,vals2);
        p(a) = k;
        
        t = title(ax(a),sprintf('%s, mean=%1.3f N=%d vs. mean=%1.3f N=%d, p=%1.6f',titleStr{a},mean(vals1),length(unique(pooledPolSelValues(2).flies)),mean(vals2),length(unique(pooledPolSelValues(4).flies)),p(a) ));
        L=legend(titleStr{1:2});
        
         end
    end
    
    set(ax(a),'xlim',[0 1])
    xlabel(ax(a),'pol selectivity index')
    ylabel(ax(a),'probability')
    set(ax(a),'TickDir','out')
    pbaspect(ax(a),[1 1 1])
    t.FontSize = 10;
    ylims = ylim;
    ylim2(a) = ylims(2);

    L.Box = 'off';
    L.Location = 'southoutside';
    L.Orientation = 'horizontal';
    
    set(ax(a),'ylim',[0 max(ylim2)])

end
bentitle(sprintf('%s',objarray(1).Layers(1).Line));

% linkaxes(ax(1:4))

%%
f2 = figure('color','w');
inclFlies = intersect(meanPolSelValues(2).flies,meanPolSelValues(4).flies);
[~,exp2Idx] = ismember(inclFlies,meanPolSelValues(2).flies);
[~,exp4Idx] = ismember(inclFlies,meanPolSelValues(4).flies);
for a = 1:2
    ax(a) = subplot(1,2,a);
    hold(ax(a),'on')
    exp2data = [meanPolSelValues(2).(titleStr{a+2})(exp2Idx)];
    exp4data = [meanPolSelValues(4).(titleStr{a+2})(exp4Idx)];
    plot(ax(a),[zeros(1,length(exp2data))',ones(1,length(exp4data))']',[exp2data',exp4data']','Color',cols(2,:))
    plot(ax(a),zeros(1,length(exp2data)),exp2data,'o','MarkerSize',4,'MarkerFaceColor',cols(a+2,:),'MarkerEdgeColor',cols(a+2,:))
    plot(ax(a),ones(1,length(exp4data)),exp4data,'o','MarkerSize',4,'MarkerFaceColor',cols(a+2,:),'MarkerEdgeColor',cols(a+2,:))
    
%     cohens_d = nanmean( exp4data - exp2data ) / nanstd( exp4data - exp2data )
    xlim([0 1])
    ylim([0 1])
    [~,t2p] = ttest2(exp2data,exp4data);
%     t =  title(sprintf('%s, p=%1.5f',titleStr{a+2},t2p));
    t = title(ax(a),sprintf('%s, mean=%1.3f N=%d vs. mean=%1.3f N=%d, p=%1.6f',titleStr{a+2},mean(exp2data),length(exp2data),mean(exp4data),length(exp4data), t2p ));

    t.FontSize = 10;

    ax(a).XTick = [0 1];
    ax(a).XTickLabel = {titleStr{1};titleStr{2}};
    ax(a).XTickLabelRotation = 30;
    
    ax(a).YTick = [0:0.25:1];
    
    if a ==1
        ax(a).YAxis.Label.String = 'polarization selectivity index';
    end
    pbaspect([1 1.5 1])
    
end
bentitle(sprintf('%s, N=%d',objarray(1).Layers(1).Line,length(exp4Idx)));
% offsetAxesXoff(ax(1))
% offsetAxesXoff(ax(2))
%%
%{
% 
% a = pooledPolSelValues(2).background;
% [~,maxidx] = max(a);
% a(maxidx) = median(a);
% k1 = kruskalwallis([pooledPolSelValues(2).background,a],{'A';'A'''});


histogram(ax(2),pooledPolSelValues(4).background,'binwidth',1/50,'normalization','probability')
histogram(ax(2),pooledPolSelValues(4).cells,'binwidth',1/50,'normalization','probability')
set(ax(2),'xlim',[0 1])
xlabel(ax(2),'pol selectivity index')
ylabel(ax(2),'probability')
set(ax(2),'TickDir','out')

k2 = kruskalwallis([pooledPolSelValues(4).background,pooledPolSelValues(4).cells],{'background';'cells'});

title(ax(2),sprintf('with polarizer, N=%d, p=%1.3f',length(unique(pooledPolSelValues(4).flies)),k2 ))

pbaspect(ax(2),[1 1 1])
%}

%{
d1 = [maskedPolSelFrameVector(topQuart)];
[~,t1] = ttest([maskedPolSelFrameVector(bottomQuart)],[maskedPolSelFrameVector(topQuart)])




maskedAvgFrameVector = x(objidx_test).MIP.AverageFrame(x(objidx_test).MIP.layerMask.mask);
[~,frameIdx] = sort(maskedAvgFrameVector,'descend');
topQuart = frameIdx(1:floor(length(frameIdx)*intensityFraction));
bottomQuart = frameIdx(ceil(length(frameIdx)-length(frameIdx).*intensityFraction)+1:end);

maskedPolSelFrameVector = x(objidx_test).MIP.polSelImg(x(objidx_test).MIP.layerMask.mask);

p2 = ranksum([maskedPolSelFrameVector(bottomQuart)],[maskedPolSelFrameVector(topQuart)],'tail','left')
k2 = kruskalwallis([maskedPolSelFrameVector(bottomQuart),maskedPolSelFrameVector(topQuart)])
d2 = [maskedPolSelFrameVector(topQuart)];
[~,t2] = ttest([maskedPolSelFrameVector(bottomQuart)],[maskedPolSelFrameVector(topQuart)])

subplot(1,2,2)
hold on
histogram(maskedPolSelFrameVector(bottomQuart),'binwidth',1/50,'normalization','probability')
histogram(maskedPolSelFrameVector(topQuart),'binwidth',1/50,'normalization','probability')
xlim([0 1])
xlabel('pol selectivity index')
ylabel('probability')
% L=legend('background ROI','cell ROI');
% L.Box = 'off';
% L.Location = 'northeast';
set(gca,'TickDir','out')
title('with polarizer')
pbaspect([1 1 1])


p3 = ranksum(d1,d2,'tail','left')
%}

