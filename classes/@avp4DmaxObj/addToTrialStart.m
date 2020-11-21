function addToTrialStart(objarray,numSec)
assert(nargin == 2,'A number of seconds must be input as numSec')
assert(isnumeric(numSec),'numSec argument must be a numerical value (time, in seconds)')

	for oidx = 1:length(objarray)
		if isempty(objarray(oidx).TrialStartFrame)
		getParameters(objarray(oidx))
		objarray(oidx).Daq = [];
		end
		
		% Shift trial start samples:
		objarray(oidx).TrialStartSample = objarray(oidx).TrialStartSample - numSec*objarray(oidx).AIrate;
		
		% Find new trial start frames:
		for tidx = 1:length(objarray(oidx).TrialStartSample)
			if objarray(oidx).TrialStartSample(tidx) < 1
				error('numSec is too large: non-existent samples requested, below sample(1)')
% 			elseif tidx > 1 && objarray(oidx).TrialStartSample(tidx) <objarray(oidx).TrialEndSample(tidx-1)
% 				error('numSec is too large: trial start times are overlapping ends of previous trials')
			elseif objarray(oidx).TrialStartSample(tidx) >= objarray(oidx).TrialEndSample(tidx)
				error('numSec is too large: trial start time cannot be greater than end time')
			else
			tryFrame = nextFrame(objarray(oidx),objarray(oidx).TrialStartSample(tidx));
				if tryFrame > 0
					objarray(oidx).TrialStartFrame(tidx) =  tryFrame;
				end
			end
		end
		
	end
	disp(['Done. ', num2str(numSec) 's added at start of each trial'])
end
