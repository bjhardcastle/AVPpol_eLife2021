classdef avp4DsuperObj < handle_light
    % Parent class for avp4D... objects
    % Objects are made with a folder
    % Within that folder, all MIP tiffs will be found, along with individual
    % layers, and seperate avp4DmaxObj objects made for each. This super
    % object simply organizes the avp4DmaxObj arrays and facilitates
    % analysis across multiple layers.
    % If multiple channels were recorded, the GCaMP acitvity recording will
    % be prioritised, and static tdTomato channel made available for
    % registration, background subtraction, etc
    properties
        Area
        Cell
        Line
    end
    properties (Dependent)
        Link
        Active
    end
    properties
        Name
        ZSelect % defunct
        Layers    % array of avp4DmaxObj objects, one for each layer
        MIP       % single avp4DmaxObj object, MIP of layers
%         MSP       % single avp4DmaxObj object, max selectivity projection of layers
    end
    properties (Dependent)
        Masks
    end
    properties (Hidden)
        Folder
        LogFile
        TimeStr
        DateStr
        Channels  % number of
        numZPlanes
        micronsPerPixel
        micronsPerZStep
        activeChan
        activeChanStr
        staticChan
        staticChanStr
        staticMIP   % object
        staticLayers  % object array
        pSet % exp parameters saved at time of recording in .mat file
        
    end
    properties (Hidden, Dependent)
        Fly
        Exps
        containsPolMapExp
        bkgPSIstruct % across all layers
        cellPSIstruct
        PSIstruct
        checkDaqs
    end
    
    methods
        % Constructor
        function obj = avp4DsuperObj(pathIN)
            
            % If no path is specified:
            if nargin == 0 || isempty(pathIN)
                try
                    [PathName] = uigetdir('Y:\ben\avp_pol\18_data\');
                    pathContents = cellstr( ls( fullfile(PathName)));
                    % Check this has the correct subfolders for an avp4D
                    % recording
                    if any(contains(pathContents,'layers')) ...
                            ||  any(contains(pathContents,'projection_max'))
                        obj.Folder = PathName;
                        [~,obj.Name] = fileparts(PathName);
                    else
                        disp('Folder does not meet specifications for this object: check constructor')
                        disp('No object made')
                        clear obj
                        return
                    end
                catch
                    disp('No object made')
                    clear obj
                    return
                end
            end
            
            % If path is specified:
            if nargin >0 && ~isempty(pathIN)
                if isa(pathIN,'char')
                    PathName = pathIN;
                    if exist(PathName)==7
                        pathContents = cellstr( ls( fullfile(PathName)));
                        % Check this has the correct subfolders for an avp4D
                        % recording
                        if any(contains(pathContents,'layers')) ...
                                ||  any(contains(pathContents,'projection_max'))
                            obj.Folder = PathName;
                            [~,obj.Name] = fileparts(PathName);
                        else
                            disp('Folder does not meet specifications for this object: check constructor')
                            disp('No object made')
                            clear obj
                            return
                        end
                    end
                else
                    disp('Please input path to folder as a string, or leave empty')
                    disp('No object made')
                    clear obj
                    return
                    
                end
            end
            
            % Get linked files and properties:
            getLogData(obj)     % Finds .log file and extracts info about recording
            getChans(obj)       % Work out which channel was active/GCaMP, and which static/tdTom (if present)
            getRecProps(obj)    %
            getSubObjArrays(obj)% Make object arrays of tiffs in 'layers' and 'projection_max' folders
            getParameterFile(obj)
        end
        
        % Method for dependent property giving hyperlink to object path in explorer
        function Link = get.Link(obj)
            Link = ['<a href="matlab:winopen(''' obj.Folder ''')">open folder</a>'];
        end
        
        % Make unique fly number from 6-digit date (start of folder name)
        % and fly number from that day ("fly1" "fly2" etc)
        function Fly = get.Fly(obj)
            assert(isnumeric(str2double(obj.Name(1:6))), 'Check first part of foldername is date string')
            DateNum = str2double(obj.Name(1:6));
            FlyStr = regexpi(obj.Name,'(?<=(fly|Capture)\D*)(\d*)','match');
            FlyNum = str2double(FlyStr);
            Fly = DateNum*10 + FlyNum;
        end
        
        % Forward obj.activeChanStr to property with shorter name for
        % displaying in command window
        function Active = get.Active(obj),Active = obj.activeChanStr;end
        
        % Quick check of whether any layerMasks have been made
        function value = get.Masks(obj)
            value = 0;
            if ~isempty([obj.Layers.layerMask])
                value = 1;
            else
                loadLayerMasks(obj.Layers(1))
                loadLayerMasks(obj.Layers(end))
                if ~isempty([obj.Layers.layerMask])
                    value = 1;
                elseif ~isempty([obj.MIP.layerMask])
                    value = -1;
                end
            end
        end
        
        function value = get.Exps(obj)
            if isempty(obj.Layers(end).TrialSeqNum)
                getParameters(obj.Layers(end))
            end
            value = unique([obj.Layers(end).TrialSeqNum]);
        end
        function value = get.containsPolMapExp(obj)
            value = 0;
            for p = [2,4,8,10]
               if ismember(p,obj.Exps)
                  value = p; 
                  return
               end
            end
        end
        function value = get.bkgPSIstruct(obj)
            value = [];
            

            for lIdx = 1:length(obj.Layers)
                % Temporarily apply MIP layer mask to
                % individual layers
                obj.Layers(lIdx).layerMask = obj.MIP.layerMask;
                vectorOut = getPolSelMaskedVector(obj.Layers(lIdx),'bkg');
                mean(vectorOut)
                if length(vectorOut) > 1
                    value = [value;vectorOut];
                end
                
                % Restore layer mask
                loadLayerMasks(obj.Layers(lIdx))
            end
        end
        function value = get.cellPSIstruct(obj)
        value = [];
            if obj.Masks == 1
                for lIdx = 1:length(obj.Layers)
                    vectorOut = getPolSelMaskedVector(obj.Layers(lIdx),'top');
                    if length(vectorOut) > 1
                        value = [value;vectorOut];
                    end
                    
                end
            elseif obj.Masks == -1
                disp('No layer masks exist')
            else
                disp('Some layer masks missing')           
            end
        end
        
        function value = get.PSIstruct(obj)
        value = struct('cellMask',[],'bkgMask',[],'layerMaskTop',[],'layerMaskBottom',[],'layerMaskAll',[]);
            
                for lIdx = 1:length(obj.Layers)
                    % Temporarily apply MIP layer mask to
                    % individual layers
                    obj.Layers(lIdx).layerMask = obj.MIP.layerMask;
                    
                    vectorOut = getPolSelMaskedVector(obj.Layers(lIdx),'cell');
                    value.cellMask = [value.cellMask;vectorOut];
                    
                    vectorOut = getPolSelMaskedVector(obj.Layers(lIdx),'bkg');
                    value.bkgMask = [value.bkgMask;vectorOut];
                    
                    vectorOut = getPolSelMaskedVector(obj.Layers(lIdx),'top');
                    value.layerMaskTop = [value.layerMaskTop;vectorOut];
                    
                    vectorOut = getPolSelMaskedVector(obj.Layers(lIdx),'bottom');
                    value.layerMaskBottom = [value.layerMaskBottom;vectorOut];
                    
                    vectorOut = getPolSelMaskedVector(obj.Layers(lIdx),'all');
                    value.layerMaskAll = [value.layerMaskAll;vectorOut];

                    % Restore layer mask
                    loadLayerMasks(obj.Layers(lIdx))
                end
            
        end
        
        function value = get.checkDaqs(obj)
            if isempty(obj.MIP)
                value = [];
            else
            value = all(strcmp(obj.MIP.DaqFile, {obj.Layers.DaqFile}));
            end
        end 
    end
    
    methods
        getLogData(obj)
        getChans(obj)
        getSubObjArrays(obj)        
        getParameterFile(obj)

        superChangeRoot(objarray,oldroot,newroot)

        superGetPolMaps(obj)
        maskLayers(obj) % GUI for making/checking masks on each layer

        MIPsummary(obj)
        Layersummary(obj,layerSelect)
        
        
        [timeVector,fittedData] = getAlignedTimeseries(obj,ROI,expNum,inclPre,inclPost)
        F0 = findPolF0(obj,ROI,expNum)
        
        superXY(obj,Yoff,Xoff,rotateOFF,darkMode,applyTuningCols,plotIndividFits,plotPooledFit,UseMSPifLayerMasksEmpty)
        superMIPXY(obj,Yoff,Xoff,rotateOFF, darkMode)
        superPBXY(obj,Yoff,Xoff,rotateOFF, darkMode)
                superRTheta(obj,Yoff,Xoff,rotateOFF, darkMode)
        superBarXY(obj,Yoff,Xoff,rotateOFF, darkMode)
        superPCA(obj)
        superPCA2(obj)

        
        % Tools
                varargout = sort(objarray,direction)
makeMSP(obj)
        superSaveMSP(objarray,overwrite) % max selectivity projection
        superSaveMISP(objarray,overwrite) % max intensity projection (with fixed layers)
        superTimes(objarray)
        superNames(objarray)
        superVal = superFunc(obj,funcString,propString)
        superUseMSP(objarray,value)
        superPolThreshold(objarray,thresh)
        superRefresh(objarray)

        makeMISP(obj) % like MSP, but using intensity, not psi
        superPolCtrl(objarray)
        
        superBarMap(obj,layerSpec,darkMode)
        superBarPolMap(obj,layerSpec,darkMode)
        superPBTMap(obj,layerSpec,noMask)  % pol bar thrust together
        
        
        varargout = plotSinglePolarResp(objarray,roiAngs,applyTuningColorsToManualROIs,addMeanROIresultants,plotCtrl)
        outputData = plotFlashResp(objarray,roiAngs,errorbar)  % using auto pol ang ROIs (obj.MIP.polROI)
        outputData = plotFlashRespManual(objarray,ROIidx,angs,plotCtrl) % using manually made ROIs (obj.MIP.ROI)
        plotPolMapTrial(objarray,roiAngs,errorbar,numCycles)
        plotPolMapTrialManual(objarray,ROIidx,errorbar,numCycles) % using manually made ROIs (obj.MIP.ROI)
        plotPolMapTrialManualRev(objarray,ROIidx,errorbar,numCycles)
        plotPolMapCtrlTrial(objarray,roiAngs,errorbar,numCycles)
        outputData = plotFlashRespExp9(objarray,roiAngs,showErrorBar)
        plotBarMapTrial(objarray,roiAngs,errorbar,numCycles)
        plotThrustTrial(objarray,roiAngs,errorbar,numCycles)
        
        superTuningMap(obj,layerSpec,discretize,darkMode)
        superSelectMap(obj,layerSpec,darkMode)
        superImageSummary(objarray)
        
        varargout =  superPolAngHist(objarray,useTuningColors,usePSIthreshold)
        varargout = superPolAngHistROI(objarray,useTuningColors,usePSIthreshold)
superPolAngResultant(objarray,usePSIweighting,usePSIthreshold)
superPolAngResultantCirc(objarray,usePSIweighting,usePSIthreshold)
superPolAngResultantROI(objarray,usePSIweighting,usePSIthreshold)

        PSI = getPSIstruct(objarray,excludeSingleLayerRecs)
        
        % E-PG protocerebral bridge functions
        superPBTrial(objarray)
        [selOut,aopOut]  = superPBMIP(objarray)
         snapStruct = superCycleSnapshots(objarray)
         snapStruct = superCyclePolarResp(objarray,roiIdx)
         snapStruct = superCyclePolarHalfResp(objarray,roiIdx)
         superAssignPBROIs(objarray,refreshDistMask,manualLimits) % get approximate rois in PB
         setPBlims(superobj)
assignPolPixRGB(objarray)

 varargout = superTuningCurve(objarray,useCellMask)

    end
    
end