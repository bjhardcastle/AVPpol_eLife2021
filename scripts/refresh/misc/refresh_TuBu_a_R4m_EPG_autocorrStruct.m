% Get first peak of autocorrelation function of each ROI's response during
% pol mapping exp, for TuBu_a and R4m in anterior bulb and E-PGs in PB.
% Store in a format that can be fed into 'notboxplot' to display the shift
% of the peak relative to the period of the stimulus for each celltype
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathsAVP

autocorrStruct = struct('Peakshift',[],'GrpIdx',[],'StimLocked',[],...
    'Line',[],'Exp',[],'FlyNumPeriodic',[],...
    'FlyNum_AllExamined',[],'GrpIdx_AllExamined',[]);

lineCell{1} = 'R34H10_Bu';  
lineCell{2} = 'R34D03_Bu';
lineCell{3} = 'SS00096_PB';
% this order reflects x-axis positions. exp2 (control) and exp4 will be
% plot for each, side by side

for Lidx = 3:-1:1 % since this is often run while E-PG array is aready loaded
    
    lineStr = lineCell{Lidx};
    eval(['load' lineStr])
    loadROIs([x.MIP])
    
    for cIdx = 1:2 % run ctrls (exp2) and regular mapping (exp4)
        
        pmidx = Lidx*3+(cIdx)-1; % x-axis position (includes some spacing between lines)
        
        ctr = [];
        for oidx = 1:length(x)
            if ~isempty(x(oidx).MIP.polExp) && x(oidx).pSet(2+(cIdx-1)*2).trialTestPauseLength == 4 && ...
                    x(oidx).pSet(2+(cIdx-1)*2).trialRandomizeOrder == 0 && ismember(2+(cIdx-1)*2,x(oidx).Exps)
                loadROIs( x(oidx).MIP)
                
                pexp = x(oidx).MIP.polExp;
                
                % Take 2 cycles worth of timeseries data
                p1 = x(oidx).MIP.expStart(x(oidx).MIP.polExp);
                [~,closestTrialIdx] = min(abs(x(oidx).MIP.TrialStartFrame- p1) );
                p2 = x(oidx).MIP.TrialEndFrame(closestTrialIdx+11+(cIdx-1)*12); % 12 trials per cycle, hardcoded
                
                if strcmp(x(1).Line,'SS00096-gal4')
                    ROIidx = [9:24]; % left and right ROIs separately
                else
                    ROIidx = 1:length(x(oidx).MIP.ROI);
                end
                
                for n = ROIidx
                    
                    % for each glomerulus position, we get the ROI time-series,
                    resp = x(oidx).MIP.ROI(n).response;
                    resp = resp(p1:p2);
                    resp = detrend(resp);
                    fps = x(oidx).MIP.AIrate/x(oidx).MIP.IFI;
                    
                    halfCycTimeS =(180/x(oidx).MIP.pSet(pexp).polAngleStep)*(x(oidx).MIP.pSet(pexp).trialTestPauseLength + 1.33*x(oidx).MIP.pSet(pexp).trialMotorPauseLength ) ;
                    
                    maxlag = ceil(1.5*halfCycTimeS*fps);
                    
                    [a,lags] = xcorr(resp,maxlag,'coeff');
                    a_S = a(lags>0);
                    lags_S = lags(lags>0)/fps;
                    [pks,locs] = findpeaks(a_S,'MinPeakProminence',0.4);
                    
                    if length(locs)>1
                        % deal with detection of two peaks
                        [~,closestidx] =  min(abs(lags_S(locs)-halfCycTimeS)); % cyctime was prev hardcoded as 31
                        locs = locs(closestidx);
                        pks = pks(closestidx);
                        %{
                           figure
                           hold on
                           plot(lags_S,a_S,'color', x(oidx).MIP.ROI(n).color)
                           plot(lags_S(locs),pks,'r*')
                           xlabel('Lag (s)')
                           ylabel('Autocorrelation')
                        %}
                    end
                    
                    if ~isempty(locs)
                        
                        % does the peak occur at the period of the
                        % stimulus? (within +/- half the static angle
                        % presentation time)
                        pkshift = halfCycTimeS - lags_S(locs);
                        if abs(pkshift) <= 0.5*x(oidx).MIP.pSet(pexp).trialTestPauseLength
                            ctr(end+1) = 0;
                        else
                            ctr(end+1) = 1;
                        end
                        autocorrStruct.Peakshift(end+1) = pkshift; % relative to period of stimulus
                        autocorrStruct.GrpIdx(end+1) = pmidx; % x-index for notboxplot
                        autocorrStruct.StimLocked(end+1) = ctr(end); % logical, peakshift within limits to indicate same periodicity as stimulus
                        autocorrStruct.Line{end+1} = [lineStr]; % array label
                        autocorrStruct.Exp(end+1) = 2+(cIdx-1)*2; % expnum, 2 or 4
                        autocorrStruct.FlyNumPeriodic(end+1) = x(oidx).Fly;
                        
                    end
                    autocorrStruct.FlyNum_AllExamined(end+1) = x(oidx).Fly; % store the number of all flies examined
                    % Controls in particular may have so little
                    % periodicity their autocorr peak does not get
                    % detected
                    autocorrStruct.GrpIdx_AllExamined(end+1) = pmidx;
                    
                    %{
                        hold on
                        plot(lags_S,a_S,'color', x(oidx).MIP.ROI(n).color)
                        plot(lags_S(locs),pks,'r*')
                        xlabel('Lag (s)')
                        ylabel('Autocorrelation')
                    %}
                    
                end
            end
        end
        %         percentageNonCorrelated(pmidx) = 100*sum(ctr)/length(ctr);
    end
end
autocorrStruct.lineCell = lineCell;
save(autocorrStruct_path,'autocorrStruct');
