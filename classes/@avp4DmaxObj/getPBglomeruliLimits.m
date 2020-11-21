function gLimits = getPBglomeruliLimits(obj,refreshLimits)
if isempty(obj.PBglomLims)
    loadPBGlomeruliLimits(obj)
end
if nargin<2 || isempty(refreshLimits)
    refreshLimits=0;
end
if ~refreshLimits && ~isempty(obj.PBglomLims)
    gLimits = obj.PBglomLims;
else
%     setPBlims(obj)
%     gLimits = getPBglomeruliLimits(obj,0);
end