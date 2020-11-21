% Generate plots for Fig8 (E-PG responses in protocerebral bridge)
% Each cell can be run independently

pathsAVP
if exist(fullfile(fig8path),'dir')
    try rmdir(fullfile(fig8path),'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tuning and selectivity maps in PB
loadSS00096_PB
printpath = fig8path;

superUseMSP(x,1)
superPolThreshold(x,-1); 

useTuningCols = 1;
useLayerMask = 1;


for sIdx = [7,1,3,4]
    suffix = num2str(sIdx);
    
    prefix = 'polSel';
    plotPolSelImg( x(selectObj(sIdx)).MIP ,1,[8+1:8+16])
    printAVP
    
    prefix = 'polTuning';
    plotCombPolImgManual( x(selectObj(sIdx)).MIP ,useTuningCols,[8+1:8+16],~useLayerMask,1)
    printAVP
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSI values in E-PG (and R4m EB for comparison)
pathsAVP
printpath = fig8path;
prefix = 'psi_';

savename = 'EPG_R4m';
objStr = {'SS00096_PB';'R34D03_EB'};

% probability density
pdfplotPSI(objStr,'add','none');
% legend off
suffix = 'pdf';
printAVP

%abs values
b = boxplotPSI(objStr);
suffix = 'box';
printAVP

% mean-of-controls-subtracted values
b = boxplotPSI(objStr,'pol','compare');
suffix = 'diff';
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example glomerulus-pair ROI traces 
loadSS00096_PB
printpath = fig8path;

exampleRec = selectObj(3);
exampleROIs = [2,6];

plotPolMapTrialManual(x(selectObj(3)),[exampleROIs+8 exampleROIs+16])

prefix = 'polMapping_';
suffix = ['ROI_glom[' num2str(exampleROIs(1)) ',' num2str(exampleROIs(2))  ']'];
printAVP


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example cycle tuning vector grid with PSI colormap 
loadSS00096_PB
printpath = fig8path;

exampleRec = selectObj(3);
exampleROIs = [2,6];

load(snapshotStruct_HalfCycle_path('SS00096_PB'))

for n = length(snapStruct):-1:1

    if ~isempty(snapStruct(n).individual) && strcmp(num2str(snapStruct(n).individual.ID),[x(exampleRec).DateStr x(exampleRec).TimeStr] )        
        tuning = flipud(snapStruct(n).individual.ROImeanTuningAng);
        psi = flipud(snapStruct(n).individual.ROIpsi);
        rowStr = cellstr(num2str([size(psi,1):-1:1]'));
        ax = pbgridplot(tuning,psi,rowStr,0);
%         ax.YAxis.TickLabel = {strcat( ['cycle'], num2str([length(rowStr):-1:1]') )};
        ax.YAxis.TickLabel = {};
        ax.YAxis.TickLabel{1} = ['cycle ' num2str(length(rowStr))];
        ax.YAxis.TickLabel{length(rowStr)} = 'cycle 1';
        r1 = rectangle(ax,'Position',[exampleROIs(1)-1,0,1,length(rowStr)],'FaceColor','none','EdgeColor',ROIcolor(exampleROIs(1)));
        r2 = rectangle(ax,'Position',[exampleROIs(2)-1,0,1,length(rowStr)],'FaceColor','none','EdgeColor',ROIcolor(exampleROIs(2)));
    end
end
% getAVPplotParams
% setAVPaxes(ax,defaultAxisHeight_cm,defaultAxisHeight_cm)
tightfig(gcf)
prefix = 'regular';
savename = 'PBgrid';
suffix = ['PSI_ROI_glom[' num2str(exampleROIs(1)) ',' num2str(exampleROIs(2))  ']'];
printAVP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% All panels using single cycle responses 
plot_PB_single_cycle