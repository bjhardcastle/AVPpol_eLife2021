function superPolThreshold(objarray,thresh)
assert(isnumeric(thresh),'Input a threshold value')
if thresh < 0
    % use avg pol sel value in the background as threshold
    useBkgThresh = 1;
    stdMultiplier = abs(thresh); % added in for testing 1 vs 2 std.
else
    useBkgThresh = 0;
end

if useBkgThresh && strcmp(objarray(1).Line,'SS00096-gal4')
    pbthresh = getPBCtrlPolSelThreshold(stdMultiplier);
end
for oidx = 1:length(objarray)
    if objarray(oidx).containsPolMapExp
        for lidx = 1:length(objarray(oidx).Layers)
            if useBkgThresh
                switch objarray(oidx).Line
                    case 'SS00096-gal4'
                        thresh = pbthresh;
                    otherwise
                        thresh = getBkgPolSelThreshold(objarray(oidx).Layers(lidx) , stdMultiplier);
                end
            end
            objarray(oidx).Layers(lidx).polSelThreshold = thresh;
        end
        if ~isempty(objarray(oidx).MIP)
            if useBkgThresh
                switch objarray(oidx).Line
                    case 'SS00096-gal4'
                        thresh = pbthresh;
                    otherwise
                        thresh = getBkgPolSelThreshold(objarray(oidx).MIP , stdMultiplier);
                end
            end
            objarray(oidx).MIP.polSelThreshold = thresh;
        end
    end
end
end
function thresh = getPBCtrlPolSelThreshold(stdMultiplier)
PSI = loadPSIstruct;
a={PSI.SS00096_PB.ctrlCellData};
p = [];
for n = 1:length(a)
    p = [p;a{n}];
end
thresh = nanmean(p) + stdMultiplier*nanstd(p);
end
function thresh = getBkgPolSelThreshold(obj,stdMultiplier)
%{
% now incorporated into obj.bkgMask
bkgintensityfraction = 0.10; % ie 10% dimmest pixels (over whole exp) constitute background
avgimg = obj.avgPolImg;
% "crop" sides by some % of frame dimension to try to deal with registration
% blanks which cause high pol Sel values (~0.5)
cropFraction = 0.075;
yCrop = ceil(size(avgimg,1).*cropFraction);
xCrop = ceil(size(avgimg,2).*cropFraction);
avgimg(1:yCrop,:) = 2^16;
avgimg(end-yCrop+1:end,:) = 2^16;
avgimg(:,1:xCrop) = 2^16;
avgimg(:,end-xCrop+1:end) = 2^16;
[bkg,sortIdx] = sort(avgimg(:),'ascend');
idx = sortIdx(1:floor(length(bkg)*bkgintensityfraction));
mask = zeros(size(obj.avgPolImg));
mask(idx) = 1;
%}
thresh = mean(obj.polSelImg(obj.bkgMask))+stdMultiplier*std(obj.polSelImg(obj.bkgMask));
if thresh>0.75
    thresh = 0.75;
end
end