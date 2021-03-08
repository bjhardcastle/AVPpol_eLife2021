% Generate plots for Appendix1Fig1 (overall pathway diagram, tuning curves, psi vals )
% Each cell can be run independently

pathsAVP
thisfigpath = figX1path;
if exist(thisfigpath,'dir')
    try rmdir(thisfigpath,'s'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSI vals (All cells)
% printpath = fig10s1path;
pathsAVP
printpath = figX1path;

prefix = 'psi_';

savename = 'pathway_';
objStr = {...
    'SS00096_PB';...
    '';'R34D03_EB';...
    '';'R49E09_Bu';'R88A06_Bu';'R88A06_Bu_ant';'R34H10_Bu';'R34D03_Bu';'R19C08_Bu';...
    '';'R49E09_AOTU';'R88A06_AOTU';'R34H10_AOTU';'R17F12_AOTU';'R73C04_AOTU';'R56F07_AOTU';...
    '';'DmDRA';...
    };

% probability density
% pdfplotPSI(objStr,'add','none');
% % legend off
% suffix = 'pdf';
% printAVP

%abs values
% cells only (brightest pix)
b = boxplotPSI(objStr);
suffix = 'box';
printAVP
%{
% layer mask (all pix within)
b = boxplotPSI(objStr,'dist','mask');
suffix = 'layerMask';
printAVP
%}

% mean-of-controls-subtracted values (where controls exist)
% cells only (brightest pix)
b = boxplotPSI(objStr,'pol','compare');
suffix = 'diff';
printAVP
%{
% layer mask (all pix within)
b = boxplotPSI(objStr,'pol','compare','dist','mask');
suffix = 'diff_layerMask';
printAVP
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tuning curves for all cells 
pathsAVP
printpath = figX1path;
getAVPplotParams
alpha = 0.3;
for nidx = 1:size(objnames,1)
    eval(['load' objnames{nidx}])
    superUseMSP(x,-1)
    
    %{
    useCellMask = 0;
    output = superTuningCurve(x,useCellMask);
    for lidx = 1:length(output.Lines)
        output.Lines(lidx).Color = [objCols.(objnames{nidx}) alpha];
    end
    try 
       output.meanLine.Color = [objCols.(objnames{nidx}) 0.9];
    catch
    end
    prefix = 'tuning_curve';
    savename = objnames{nidx};
    suffix = '_layerMask';
    printAVP
    %}
    
    useCellMask = 1;
    output = superTuningCurve(x,useCellMask);
    for lidx = 1:length(output.Lines)
        output.Lines(lidx).Color = [objCols.(objnames{nidx}) alpha];
    end
    try
        output.meanLine.Color = [objCols.(objnames{nidx}) 0.9];
    catch
    end
    prefix = 'tuning_curve';
    savename = objnames{nidx};
    suffix = '_cellMask';
    printAVP

end


