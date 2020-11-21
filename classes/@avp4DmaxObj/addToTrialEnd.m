function addToTrialEnd(objarray,numSec)
assert(nargin == 2,'A number of seconds must be input as numSec')
assert(isnumeric(numSec),'numSec argument must be a numerical value (time, in seconds)')

	for oidx = 1:length(objarray)
		if isempty(objarray(oidx).TrialEndFrame)
		getParameters(objarray(oidx))
		objarray(oidx).Daq = [];
		end
		
		% Shift trial end samples:
		objarray(oidx).TrialEndSample = objarray(oidx).TrialEndSample + numSec*objarray(oidx).AIrate;
		
		% Find new trial end frames:
		for tidx = 1:length(objarray(oidx).TrialEndSample)
			if objarray(oidx).TrialEndSample(tidx) > objarray(oidx).Frametimes(end)
				error('numSec is too large: non-existent samples requested, past last frame')
% 			elseif tidx <length(objarray(oidx).TrialEndSample) && objarray(oidx).TrialEndSample(tidx)>objarray(oidx).TrialEndSample(tidx+1)
% 				error('numSec is too large: trial end times are overlapping starts of next trials')
			elseif objarray(oidx).TrialEndSample(tidx) <= objarray(oidx).TrialStartSample(tidx)
				error('numSec is too large: trial end time cannot be less than start time')
			else
			tryFrame = prevFrame(objarray(oidx),objarray(oidx).TrialEndSample(tidx));
				if tryFrame > 0
					objarray(oidx).TrialEndFrame(tidx) =  tryFrame;
				end
			end
		end
		
	end
	disp(['Done. ', num2str(numSec) 's added at end of each trial'])
end
