function [selOut,aopOut] = superPBMIP(objarray)
% Rough code for checking cycle-by-cycle variation in selectivity index and
% tuning angle
% We just keep a selection of trials from the MIP object and run the usual
% anaylses (for pol tuning /selectivity) and quantify the average value
% across pixels or the change from the previous cycle 
fwdIncreasingSel = nan(length(objarray),8);
bkIncreasingSel = fwdIncreasingSel;
eachIndividSel = fwdIncreasingSel;


aopOut = nan(length(objarray),8,3); % for fwd, bk, individ AOP diffs, up to 8 half cycles

for oidx = 1:length(objarray)
    if ismember(4,objarray(oidx).Exps) ...
            && ( ~isfield('polOffBetweenTrials',objarray(oidx).MIP.pSet(4)) || objarray(oidx).MIP.pSet(4).polOffBetweenTrials == 0 ) ...
            && ( ~isfield('trialRandomizeOrder',objarray(oidx).MIP.pSet(4)) || objarray(oidx).MIP.pSet(4).trialRandomizeOrder == 0 )
        
        %  examine every half cycle
        numCycles =  objarray(oidx).MIP.pSet(4).trialReps*2;
        %            getTrialtimes(objarray(oidx).MIP)
        %            getParameters(objarray(oidx).MIP)
        getPolMaps(objarray(oidx).MIP)
        trialsPerCycle = 180/objarray(oidx).MIP.pSet(4).polAngleStep;
        
        objarray(oidx).MIP.polSelHistFraction = 0.15;
        maskedPolSelFrameVector = getPolSelMaskedVector(objarray(oidx).MIP,'top');
        initVal = median(maskedPolSelFrameVector);
        
        numTrialsOrig = length(find(objarray(oidx).MIP.TrialSeqNum==4));
        
        fwdIncreasingAOP = [];
        bkIncreasingAOP = [];
        eachIndividAOP = [];
        
        try
            for n = 1:numCycles
                keepExpTrials(objarray(oidx).MIP,4,1:(n)*trialsPerCycle)
                if n==1
                    getPolMapsHalfCycle(objarray(oidx).MIP)
                else
                    getPolMaps(objarray(oidx).MIP)
                end
                %                    superImageSummary(objarray(oidx))
                %superPolCtrl(objarray(oidx))
                restoreAllTrials(objarray(oidx).MIP)
                
                maskedPolSelFrameVector = getPolSelMaskedVector(objarray(oidx).MIP,'top');
                fwdIncreasingSel(oidx,n) = median(maskedPolSelFrameVector);
                
                fwdIncreasingAOP(:,n) = getPolTuningMaskedVector(objarray(oidx).MIP,'top');
                
            end
            
            for n = numCycles:-1:1
                keepExpTrials(objarray(oidx).MIP,4,numTrialsOrig-(n*trialsPerCycle)+1:numTrialsOrig)
                if n==1
                    getPolMapsHalfCycle(objarray(oidx).MIP)
                else
                    getPolMaps(objarray(oidx).MIP)
                end
                %                    superImageSummary(objarray(oidx))
                %superPolCtrl(objarray(oidx))
                restoreAllTrials(objarray(oidx).MIP)
                
                maskedPolSelFrameVector = getPolSelMaskedVector(objarray(oidx).MIP,'top');
                bkIncreasingSel(oidx,n) = median(maskedPolSelFrameVector);
                
                bkIncreasingAOP(:,n) = getPolTuningMaskedVector(objarray(oidx).MIP,'top');
                
            end
            
            for n = 1:numCycles
                keepExpTrials(objarray(oidx).MIP,4,[1:trialsPerCycle] + (n-1)*trialsPerCycle)
                
                getPolMapsHalfCycle(objarray(oidx).MIP)
                
                
                %                    superImagSummary(objarray(oidx))
                %superPolCtrl(objarray(oidx))
                restoreAllTrials(objarray(oidx).MIP)
                
                maskedPolSelFrameVector = getPolSelMaskedVector(objarray(oidx).MIP,'top');
                eachIndividSel(oidx,n) = median(maskedPolSelFrameVector);
                
                eachIndividAOP(:,n) = getPolTuningMaskedVector(objarray(oidx).MIP,'top');
                
               plotCombPolImg(objarray(oidx).MIP)
                
            end
            objarray(oidx).MIP.Frames=[];
            objarray(oidx).MIP.Daq=[];
            
            
            for cycIdx = 2:n % calculating diff
                aopOut(oidx,cycIdx,1) = abs( rad2deg( nanmedian( wrapToHalfPi( circ_dist( deg2rad(fwdIncreasingAOP(:,cycIdx-1)), deg2rad(fwdIncreasingAOP(:,cycIdx)) )  ) ) ));
                aopOut(oidx,cycIdx,2) = abs( rad2deg( nanmedian( wrapToHalfPi( circ_dist( deg2rad(bkIncreasingAOP(:,cycIdx-1)), deg2rad(bkIncreasingAOP(:,cycIdx)) )  ) ) ));
                aopOut(oidx,cycIdx,3) = abs( rad2deg( nanmedian( wrapToHalfPi( circ_dist( deg2rad(eachIndividAOP(:,cycIdx-1)), deg2rad(eachIndividAOP(:,cycIdx)) )  ) ) ));
            end
            
            
        catch ME
            restoreAllTrials(objarray(oidx).MIP)
            rethrow(ME)
        end
        
    end
end

selOut(:,:,1)=fwdIncreasingSel;
selOut(:,:,2)=bkIncreasingSel;
selOut(:,:,3)=eachIndividSel;

dirInd = zeros(1,length(objarray));
for oidx =  1:length(objarray)
if isfield(objarray(oidx).MIP.pSet(4),'StepDIR') % if not, stepdir is 0
    dirInd(oidx) = objarray(oidx).MIP.pSet(4).StepDIR;
end
end


figure,

ax(1)= subplot(2,3,1);
plot(fwdIncreasingSel','color',ROIcolor(1))
hold on
plot(nanmean(fwdIncreasingSel,1),'b','LineWidth',2)

ax(2)= subplot(2,3,2);
plot(bkIncreasingSel','color',ROIcolor(2))
hold on
plot(nanmean(bkIncreasingSel,1),'r','LineWidth',2)

ax(3)= subplot(2,3,3);
plot(eachIndividSel','color',ROIcolor(8))
hold on
plot(nanmean(eachIndividSel,1),'k','LineWidth',2)
ylabel('pol sel shift individual cycles')

ax(4)= subplot(2,3,4);
plot((aopOut(:,:,1))','color',ROIcolor(1))
hold on
plot(nanmean(aopOut(:,:,1),1),'b','LineWidth',2)
ylabel('tuning shift using with increasing cycles (from start)')

ax(5)= subplot(2,3,5);
plot((aopOut(:,:,2))','color',ROIcolor(2))
hold on
plot(nanmean(aopOut(:,:,2),1),'r','LineWidth',2)
ylabel('tuning shift using with increasing cycles (from last)')

ax(6)= subplot(2,3,6);
hold on
% plot(aopOut((dirInd == 1),:,3)','color',ROIcolor(8))
% plot(aopOut((dirInd == 5),:,3)','color',ROIcolor(2))
% plot(nanmean(aopOut((dirInd == 1),:,3),1),'k','LineWidth',2)
% plot(nanmean(aopOut((dirInd == 5),:,3),1),'r','LineWidth',2)
plot(aopOut(:,:,3)','color',ROIcolor(8))
plot(nanmean(aopOut(:,:,3),1),'k','LineWidth',2)
xlabel('cycle')
ylabel('tuning shift compared to previous')
linkaxes(ax(1:3),'x')

end


function lambda = wrapToHalfPi(lambda)
%wrapToHalfPi Wrap angle in radians to [0 0.5*pi]
%
% modified from built-in matlab function:
%   lambdaWrapped = wrapTo2Pi(LAMBDA) wraps angles in LAMBDA, in radians,
%   to the interval [0 2*pi] such that zero maps to zero and 2*pi maps
%   to 2*pi. (In general, positive multiples of 2*pi map to 2*pi and
%   negative multiples of 2*pi map to zero.)
%
%   See also wrapToPi, wrapTo180, wrapTo360.

% Copyright 2007-2008 The MathWorks, Inc.

% WrapToPi:
q = (lambda < -0.5*pi) | (0.5*pi < lambda);
lambda(q) = wrapToPi(lambda(q) + pi);
end
