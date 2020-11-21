classdef avp4DmaxObj < SlidebookObj
    
    
    properties (Dependent)
        % %     Link    % Link to open Folder in Windows File Explorer - depends on 'Folder'
        DaqFolder % Daq files are stored separately to .tiffs until an object is made - they are then copied to the folder with .tiffs
        ZPlane
        
        
        Line % gal4 line
        Cell % celltype recorded
        Area % Area of brain recorded
        fly  % alternative to 'Fly'. auto-generated from filename, not including date
        
    end
    
    properties
        
        % From log file:
        numZPlanes
        Channels
        micronsPerPixel
        micronsPerZStep
        logTimePoints
        
        %From getparameters
        LEDvolt
        
        ROI4D
        
        polROI
        fullROI  % ROIs for all above threshold pixels for a particular angle
        % Regular ROIs are for a subset of these pixels (the largest connected
        % area)
        barROI
        fullAngMaskImg
        exampleAngMaskImg
        
        k
        
        polAngImg
        polSelImg
        
        fftAngImg
        fftMagImg
        
        polMagThreshold = 25;
        polSelThreshold = 0.35;
        
        avgPolImg % median
        stdPolImg
        polActivityImg % no longer used
        polTuningImg % mean response tuning for each pixl
        varImg % mean x std across all frames
        
        InactivityFrame% mean outside of experiments 
        
        polAngResp  % each pixel's mean resp at each angle presented
        polAngRespPooled % as above, pooled for angles+180
        polVecLength
        colorblindMode = 0;
        
        layerMask
        
        
        barMaxImg   % average of max activity in each rep of bar position (3Darray)
        barAvgImg
        barStdImg
        barTuningImg
        barSelMap
        
        thrustMaxImg
        thrustAvgImg
        thrustStdImg
        
        UseMSP % Take frames from the max selectivity projection (stored in .mat file)
        % note: MSP only exists for multilayer recs which incl a pol exp. getFrames will fail if set to useMSP and frames don't exist
        % use superUseMSP() on super obj arrays to be safe.
                
        polSelHistBins = linspace(0,1,51); % histogram bins
        polSelHistFraction = 0.1; % fraction of total mask pixels counted in top/bottom portions
        
        storedTrialVariables  % temporary storage of all trial data % see restoreExpTrials 
        % %     Fly     % Fly number. A single fly can be linked to more than one object
        % %     File    % Tiff file name, without extension
        % %     Folder  % Tiff file folder path
        % %     DateStr % Tiff capture date: 'yymmdd'
        % %     TimeStr % Tiff capture time: 'hhmm'
        % %     DaqFile % DAQ data filename: aq_EXP_MASTER_20{DateStr}_{TimeStr}.daq
        
     snapStruct = struct(); % structure with snapshots of pol sel/tuning maps using different combinations of pol exp cycles, for investigating their stability over time

        altMaskStr  = '';    % adding a string here causes an alternative 
% layerMask to be loaded/saved, without affecting anything else.
                             % Used for analysing ant/sup bulb responses
                             % separately in R88A06. Note: new alternative
                             % polROIs can be saved/loaded, but regular
                             % ROIs cannot (superclass function)

     limitCycles = 1; % enabled by default, this limits the extent of pol experiments analysed with getPolMaps, in order to give a fair comparison of PSI values across experiments
     PBcurvilinearDistMask
     PBglomLims
    end
    
    properties (Hidden = true, Dependent = true)
        polPix % Main array of pol tunings by pixel, after applying amplitude/polSelecitivity thresholds and masks to exclude cell bodies etc
        polPixDiscrete % as above, discretized
        bkgMask
        cellMask
        
        barPix
        barSelImg
        thrustPix
        thrustSelImg
        polExp  % (2,4,8,10) reference to pol experiment protocol
       
        SpecifiedDaqTime % Some recordings made in a loop generated DAQ files with times which made them difficult to match
        % - the solution was usually to use the daq file prior to the one
        % matched in 'getDaqFile'.. this option enforces a specified time
        % string at the end of the daq file name (getting the date wrong
        % wouldn't be possible)

    end
    
    properties (Transient)
        % Only stored temporarily. When the object is saved to disk, these
        % properties are not saved.
        % %
        % %     Daq             % Timeseries data extracted from DaqFile
        % %     Frames          % Tiff frames. Clear to save memory
        % %     Unattended = 0  % Set to 1 to run functions without requesting user input
    end
    
    properties (Transient = true, SetAccess = protected)
        % As above, but can only be modified by a class method
        % %
        % %     BackgroundSubtracted = 0  % Read-only marker that indicates background has been subtracted from stored Frames, if set to 1.
        % %
    end
    
    properties (Hidden)
        pSet  % from parameters file
    end
    
    methods
        %%%% Constructor %%%%
        
        function obj = avp4DmaxObj(pathIN)
            
            if nargin == 0
                pathIN = [];
            end
            
            % Call superclass constructor:
            obj@SlidebookObj(pathIN);
            
            % For some recordings Slidebook crashed but frame markers
            % continued to be recorded, here we just cut short the amount
            % of the DAQ file we look at for recorded stuff
            if strcmp([obj.DateStr obj.TimeStr], '1810091409')
                obj.TrialSettings.setlimits = [1 2130000];
            end
            if strcmp([obj.DateStr obj.TimeStr], '1902191025')
                obj.TrialSettings.setlimits = [1 2152000];
            end
            if strcmp([obj.DateStr obj.TimeStr], '1810091416')
                obj.TrialSettings.setlimits = [1 805600];
            end
            if strcmp([obj.DateStr obj.TimeStr], '1810041300')
                obj.TrialSettings.setlimits = [1 4734659];
            end
            
            % Run some additional code on creation:
            try
                getParameterFile(obj)
                
                % Override default values:
                obj.TrialSettings.chan = 3;
                loadLayerMasks(obj)
            catch ME
                msgText = getReport(ME);
                disp('Constructor error:')
                disp(msgText)
            end
            
        end
        
    end
    
    
    methods
        %%% Functions for dependent variables:
        
        function value = get.ZPlane(obj)
            % Get imaging plane of this file:
            plane = str2double( regexp(obj.File,'(?<=.*_Z)\d*(?=_C)','match'));
            if ~isempty(plane) && ~isnan(plane)
                % Range from Slidebook export is 0:numZPlanes-1
                % Here we convert to 1:numZPlanes
                value = plane + 1;
            else
                value = [];
            end
        end
        
        function value = get.DaqFolder(obj)
            % .daq files are stored along with .sld file.
            % When an object is made, getDaqFile uses the DaqFolder
            % parameter to locate the .daq file, then copies it to the same
            % location as the .tiff file. After that, the .daq file can be
            % located as usual
            value = char(strcat(regexpi(obj.Folder,'.*18_data\','match'), 'dump\', obj.DateStr, '\'));
        end        
        function value = get.Cell(obj)
            parts = strsplit(obj.Folder,filesep);
            value = parts{end-5};
        end
        function value = get.Area(obj)
            parts = strsplit(obj.Folder,filesep);
            areaStr = parts{end-6};
            value = areaStr(3:end);
        end
        function value = get.Line(obj)
            parts = strsplit(obj.Folder,filesep);
            value = parts{end-4};
        end
        function value = expStart(obj,expNum)
            if isempty(obj.TrialSeqNum)
                getParameters(obj)
            end
            value = obj.TrialStartFrame(find(obj.TrialSeqNum == expNum,1,'first'));
        end
        % Make unique fly number from 6-digit date (start of folder name)
        % and fly number from that day ("fly1" "fly2" etc)
        function fly = get.fly(obj)
            FlyStr = regexpi(obj.File,'(?<=(fly|Capture)\D*)(\d*)','match');
            fly = str2double(FlyStr);
        end
        
        function value = expStop(obj,expNum)
            if isempty(obj.TrialSeqNum)
                getParameters(obj)
            end
            value = obj.TrialEndFrame(find(obj.TrialSeqNum == expNum,1,'last'));
        end
        function value = expPolAngs(obj,expNum)
            if isempty(obj.TrialSeqNum)
                getParameters(obj)
            end
            value = unique(obj.TrialPatNum(obj.TrialSeqNum == expNum));
        end
        function value = polColMap(obj)
            
            switch obj.colorblindMode
                case 0
                    if isempty(obj.TrialSeqNum)
                        getParameters(obj)
                    end
                    if any(obj.TrialSeqNum == 4)
                        value = flipud(hsv(0.5*length(obj.expPolAngs(4))));
                    elseif  any(obj.TrialSeqNum == 2)
                        value = flipud(hsv(0.5*length(obj.expPolAngs(4))));
                    else
                        value = flipud(hsv(6));
                    end
                    value(size(value,1)+1,:) = ROIcolor(8); % Gray for control ROI in last place
                    
                case 1 % Assumes 6 angles + control
                    cbfMap(1,:) = [66,126,148]./255;
                    cbfMap(2,:) = [196,100,46]./255;
                    cbfMap(3,:) = [134,193,157]./255;
                    cbfMap(4,:) = [255,219,56]./255;
                    cbfMap(5,:) = [149,139,95]./255;
                    cbfMap(6,:) = [0,0,0];
                    cbfMap(7,:) = [0.4 0.4 0.4];
                    value = cbfMap;
            end
        end
        function value = get.barSelImg(obj)
            if isempty(obj.barMaxImg)
                getBarMaps(obj);
            end
            % Selectivity: (avgmax in pref bar pos - mean avgmax in all other bar pos)/pref
            [maxframe,~] = max(obj.barMaxImg,[],3);
            s = sort(obj.barMaxImg,3,'descend');
            value = (maxframe - mean(s(:,:,2:end),3))./maxframe;
            % Alternative is to use min(obj.barMaxImg,[],3)
        end
        function value = get.barPix(obj)
            if isempty(obj.barMaxImg)
                getBarMaps(obj);
            end
            [maxframe,idxframe] = max(obj.barMaxImg,[],3);
            % convert pref pos idx into approx bar azimuth position on screen
            pos = linspace(-90,90,5);
            posArray = pos(idxframe);
            posArray(obj.barSelImg<obj.polSelThreshold) = nan;
            if ~isempty(obj.layerMask)
                posArray(obj.layerMask.mask==0) = nan;
            end
            value = reshape( posArray , size(idxframe) );
        end
        function value = get.thrustSelImg(obj)
            if isempty(obj.thrustMaxImg)
                getThrustMaps(obj);
            end
            % Selectivity: (avgmax in pref thrust pos - mean avgmax in all other thrust pos)/pref
            [maxframe,~] = max(obj.thrustMaxImg,[],3);
            s = sort(obj.thrustMaxImg,3,'descend');
            %             value = (maxframe - mean(s(:,:,2:end),3))./maxframe;
            value = (maxframe -  mean(obj.thrustMaxImg(:,:,1:3),3))./maxframe;
            
            % Alternative is to use min(obj.thrustMaxImg,[],3)
        end
        function value = get.thrustPix(obj)
            if isempty(obj.thrustMaxImg)
                getThrustMaps(obj);
            end
            [maxframe,idxframe] = max(obj.thrustMaxImg,[],3);
            % convert pref pos idx into approx bar azimuth position on screen
            %             pos = linspace(-90,90,5);
            posArray = (idxframe);
            posArray(obj.thrustSelImg<obj.polSelThreshold) = nan;
            if ~isempty(obj.layerMask)
                posArray(obj.layerMask.mask==0) = nan;
            end
            value = reshape( posArray , size(idxframe) );
        end
        function value = get.polPix(obj)
            value = nan(size(obj.AverageFrame));
            if ~isempty(obj.polExp)
                if isempty(obj.polTuningImg)
                    getPolMaps(obj);
                end
                
                polAngImg = obj.polTuningImg;
                
                % set pixels below pol selectivity threshold to nans
                polAngImg(obj.polSelImg<obj.polSelThreshold) = nan;
                
                % if a layer mask has been drawn, set pixels outside it to nans
                if isempty(obj.layerMask)
                    loadLayerMasks(obj);
                end
                if ~isempty(obj.layerMask)
                    polAngImg(~obj.layerMask.mask) = nan;
                end
                
                value = polAngImg;
            end
        end
%         function value = get.polSelHist(obj)
%             value = getPolSelHist(obj);            
%         end

        function value = get.polPixDiscrete(obj)
            value = nan(size(obj.AverageFrame));
            if ~isempty(obj.polExp) && ~isempty(obj.polPix)
            % This section allows for discretization of the continuous
            % tuning map, for auto finding ROIs etc with common tuning.
            % Descretize func can't wrap 0/180, and since we want bin centers to
            % fall on presented angles (of which 0/180 is one), we instead have to
            % shift angle data in the range [0:0.5*binWidth] to [180:180 + 0.5*binWidth]
            binWidth = 30;
            binCenters = [0:binWidth:180-binWidth];  % length n
            shiftedBinEdges = [0.5*binWidth:binWidth:180+0.5*binWidth];  % length n+1, first/last are equivalent
            % all angle data must be in the range [shiftedBinEdges(1) : shiftedBinEdges(end)]
            Ang =0.5.*(wrapTo360((2.*(obj.polPix))));
            % shift angles from [0:180] to  [15:195] so we can look at bins centered on
            % each presented angle [0:30:180], with 165 through 15 as a single bin
            shiftAng =0.5.*(wrapTo360((2.*(Ang - shiftedBinEdges(1)))))  + shiftedBinEdges(1);
            discreteGrpVal = discretize(shiftAng,shiftedBinEdges);
            value = discreteGrpVal.*binWidth; % convert back to angle values
            end
        end
        function value = get.bkgMask(obj)
            bkgintensityfraction = 0.10;%obj.polSelHistFraction; 
            % ie 10% dimmest pixels (over whole exp) constitute background

            if ~isempty(obj.stdPolImg)
                            avgimg = obj.avgPolImg;
                        elseif ~isempty(obj.polExp)
                            mspState = obj.UseMSP;
                            obj.UseMSP = 1;
                            getFrames(obj)
                            getPolMaps(obj)
                            obj.Frames = [];
                            obj.UseMSP = mspState;
                            avgimg = obj.avgPolImg;
                        else
                            if isempty(obj.ActivityFrame1)
                                obj.Unattended = 1;
            
                                getFrames(obj);
                                obj.Frames = [];
                            end
                            avgimg = obj.ActivityFrame1;
                        end
                        
            % "crop" sides of image tight to layerMask to try to deal with registration
            % blanks which cause erroneous pol Sel values (~0.5)
            if isempty(obj.layerMask)
                cropFraction = 0.15;
                
                yCrop = ceil(size(avgimg,1).*cropFraction);
                xCrop = ceil(size(avgimg,2).*cropFraction);
                avgimg(1:yCrop,:) = 2^16;
                avgimg(end-yCrop+1:end,:) = 2^16;
                avgimg(:,1:xCrop) = 2^16;
                avgimg(:,end-xCrop+1:end) = 2^16;
            else
                m = obj.layerMask.mask;
                lCrop = find(sum(m,1),1,'first') -1; % left
                rCrop = find(sum(m,1),1,'last') + 1;  % right
                tCrop = find(sum(m,2),1,'first') -1; % top
                bCrop = find(sum(m,2),1,'last') +1;  % bottom
                avgimg(1:tCrop,:) = 2^16;
                avgimg(bCrop:end,:) = 2^16;
                avgimg(:,1:lCrop) = 2^16;
                avgimg(:,rCrop:end) = 2^16;
            end
            
            if ~isempty(obj.layerMask)
                avgimg(obj.layerMask.mask) = 2^16;
            end
            avgimg(avgimg<=10) = 2^16; %try to deal with edges 
            
            remainingPix = sum(avgimg(:)~=2^16);
            [bkg,sortIdx] = sort(avgimg(:),'ascend');
            idx = sortIdx(1:floor(remainingPix*bkgintensityfraction));
            mask = zeros(size(avgimg));
            mask(idx) = 1;
            value = logical(mask);  
            
        end

        function value = get.cellMask(obj)
        
            if ~isempty(obj.stdPolImg)
                avgimg = obj.stdPolImg;
                
            elseif ~isempty(obj.polExp)
                mspState = obj.UseMSP;
                obj.UseMSP = 1;
                getFrames(obj)
                getPolMaps(obj)
                obj.Frames = [];
                obj.UseMSP = mspState;
                avgimg = obj.stdPolImg;
                
            else
                if isempty(obj.ActivityFrame1)
                    obj.Unattended = 1;
                    
                    getFrames(obj);
                     obj.Frames = [];
                end
                avgimg = obj.ActivityFrame1;
            end
                    
            if ~isempty(obj.layerMask)
                avgimg(~obj.layerMask.mask) = 0;
            end
            
            [~,sortIdx] = sort(avgimg(:),'descend');
            idx = sortIdx(1:sum(obj.bkgMask(:))); % get the same number of pixels as the bkgMask
                        
            mask = zeros(size(avgimg));
            mask(idx) = 1;
            value = logical(mask);            
        end
        function value = get.SpecifiedDaqTime(obj)
           value = [];
           value = getStoredDaqTimes(obj); % stored 1524, 200428
           % fetches 'storedDaqTimeStr' and 'datetimeID'
           
        end
        function value = get.polExp(obj)
            value = [];
            if isempty(obj.TrialSeqNum)
                getParameters(obj)
            end
            if ismember(4,obj.TrialSeqNum)
                value = 4;  % regular pol mapping (polarizer attached)
            elseif ismember(2,obj.TrialSeqNum)
                value = 2; % regular control pol mapping (without polarizer)
            elseif ismember(8,obj.TrialSeqNum)
                value = 8; % hi-res, longer motor-pause pol mapping (polarizer attached)
            elseif ismember(10,obj.TrialSeqNum)
                value = 10; % reverse-direction pol mapping (polarizer attached)
            end
        end
        
        %%%% Sublclass functions with their own .m file %%%%
        getParameterFile(obj)
        assignStackROIs(objarray)
        [xyRectangle,zLim] = getCropBox(objarray)
        makeMIPVids(obj,exps)
        makeStackVid(objarray)
        makeColorStackVid(objarray)
        plotColorPol(obj)
        
        updateROIs(objarray) % modified to make ROI mask if missing - original SlidebookObj version could be replaced with this
        runBackSub(obj)
        
        makePolBarComparison(obj)
        varargout = getActivityFrameMax(obj,fields,detrend)
        
        getFrametimes(obj) % modified to keep every tenth frame (corresponding to z-plane 0)
        getDaqFile(obj) % Allow search for DAQ file with earlier timestamp (up to 5mins before rec) and look in ..\18_data\dump folders for .daq files, NOT folder with the tiff file
        getLogFile(obj) % Looser search for matching log file. Extract additional Zplane info
        getParameters(obj) % Mods for pol angle
        % In 4D stacks, some planes may not have been captured during
        % movement of the pattern (a longer presentation / fewer planes
        % should have been used)
        [ExpStimOnTime, ExpStimOffTime, TrialStimOnFrame, TrialStimOffFrame] = detectPanelsMovement(obj,gainStr)
        
        runTifRegStatic(obj, poolSiz, walkThresh, usfac) % Copies 'regDispXY.txt' along with registered tiff..
        runTifRegApply(obj,regDispXYfullpath) % Applies 'regDispXY.txt' specified by regDispXYfullpath to a new object tiff
        
        runAcrossPlanes(objarray,funcstring) % Share .Daq and .Frames from object 1:end and execute method
        
        assignExpActivityFrames(objarray)
        
        getInactivityFrame(objarray)

        plotVectormap(obj)
        
        F0 = findF0(obj, mask, trialData,  trialIdx, timeVector, normFPS)
        
        plotExp3ColPol(obj) % legacy
        
        plotStacks(objarray)
        varargout = getPolAngleActivityFrame(obj,angle)
        [Ftrace] = scanROI4D(objarray,maskIdx,scan_extent)
        [Ftrace] = scanROI(obj, mask, scan_extent) % Modified to allow traces full of zeros to be output
        
        varargout = plotTrials(objarray, ROImaskidx, fields, errorbar, plotcolor, tfig)
        
        varargout = plotTrialsExp3(objarray, ROImaskidx, fields, errorbar, plotcolor, tfig)
        [responseArray, timeVector, F0Array, ObjIdx] = findRespArrayExp3(objarray, ROImaskidx, fields)
        F0 = findF0Exp3(obj, mask, trialData,  trialIdx, timeVector, normFPS)
        
        varargout = plotTrialsExp4(objarray, ROImaskidx, fields, errorbar, plotcolor, tfig)
        [responseArray, timeVector, F0Array, ObjIdx,stimTimeArray] = findRespArrayExp4(objarray, ROImaskidx, fields)
        F0 = findF0Exp4(obj, mask, trialData,  trialIdx, timeVector, normFPS)
        
        varargout = plotTrialsExp4Rev(objarray, ROImaskidx, fields, errorbar, plotcolor, tfig)
        [responseArray, timeVector, F0Array, ObjIdx,stimTimeArray] = findRespArrayExp4Rev(objarray, ROImaskidx, fields)

        varargout = plotTrialsExp5(objarray, ROImaskidx, fields, errorbar, plotcolor, tfig)
        [responseArray, timeVector, F0Array, ObjIdx,stimTimeArray] = findRespArrayExp5(objarray, ROImaskidx, fields)
        F0 = findF0Exp5(obj, mask, trialData,  trialIdx, timeVector, normFPS)
       
        varargout = plotTrialsExp6(objarray, ROImaskidx, fields, errorbar, plotcolor, tfig)
        [responseArray, timeVector, F0Array, ObjIdx,stimTimeArray] = findRespArrayExp6(objarray, ROImaskidx, fields)
        F0 = findF0Exp6(obj, mask, trialData,  trialIdx, timeVector, normFPS)
  
        varargout  = findROIs(obj, numROIs, fields, detrend,saveplot) % modified to not use activity frame - copied from kmL2VisObj
        
        plotPolData(obj)
        polSel = findPolSelImg(obj)
         
         
        % findF0 is copied from the @SlidebookObj folder, and defined here
        % exactly as in the superclass. The contents of the .m file can now
        % be modified.
        % Now, whenever findF0 is called to be run on an object of this
        % subclass-type, it is the subclass version of the function that
        % will execute.
        [responseArray, timeVector, F0Array, ObjIdx] = findRespArray(objarray, ROImaskidx, fields)
        
        
        varargout = getExp3ActivityFrame(obj)
        varargout = getExp5ActivityFrame(obj)
        plotColorBars(obj)
        varargout = getExp6ActivityFrame(obj)
        plotColorDots(obj)
        varargout = getExp7ActivityFrame(obj)
        
        varargout = getBarROIs(obj)
        varargout = getPolROIs(obj, numROIs, manualROItoggle)
        varargout = getPolROIsOrig(obj, numROIs, manualROItoggle)
               
        loadPolROIs(objarray)
        savePolROIs(objarray)

        plotNatalieData(obj)
        
        plotXYcorr(obj,Yoff,Xoff,rotateOFF)
        plotRThetacorr(obj,Yoff,Xoff,rotateOFF)
        
        loadLayerMasks(objarray)
        saveLayerMasks(objarray)
        
        getPolMaps(obj) % Preferred tuning maps for exp2/exp4
        getPolMapsHalfCycle(obj) % For looking at cycle-by-cycle tuning
        avpcmap(obj)
        
        F0map = getF0map(obj) % F0 from before first trial, for exp2/exp4
        tempPlot(obj)
        
        getBarMaps(obj) % max activity maps for exp5
        getThrustMaps(obj)% max activity maps for exp6
        
        getPolMapsUnpooled(obj)
        
        varargout = getFrames(obj,scan_extent)  % modified to load MSP frames, if requested
        varargout = getPolSelMaskedVector(obj,topBottomAllBoth)
        varargout = getPolTuningMaskedVector(obj,topBottomAllBoth)

        addToTrialStart(objarray,numSec)  
        addToTrialEnd(objarray,numSec)
        
        keepExpTrials(obj,expNum,cycleNum)
        restoreAllTrials(obj)
        
        plotCombPolImg(obj,showPolCols,ROIangs,noMask)
                  plotCombPolImgManual(obj,showPolCols,ROIidx,noMask,noiseFilter)

        addScalebar(obj,ax,micronLength,invertColor)
        plotPolSelImg(obj,invertCols,ROIangs)
        
        varargout = getCurvilinearDist(obj,I,J,refreshMask) 
gLimits = getPBglomeruliLimits(obj,refreshLimits)
savePBGlomeruliLimits(objarray)
loadPBGlomeruliLimits(objarray)
value = getStoredDaqTimes(obj)
[tuningAng,resultantLength] = tuningROI(obj,ROImask,convertToPubAngs)
    end
    
end