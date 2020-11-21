% Get flash responses for three TuBu drivers in AOTU and BU
% We only store the exp9 responses (blue, UVpol+blue, UVpol; 3 reps each)
% but some remnants of code for storing exp3 responses remain (UV pol at
% 0/90deg). Responses are found for the cellmask (X% brightest pixels in
% layer mask - see avp4DmaxObj)

lineStr = { ...
    'R34H10_AOTU';...
    'R88A06_AOTU';...
    'R49E09_AOTU';...
    'R34H10_Bu';...
    'R88A06_Bu';...
    'R49E09_Bu';...
    };

% flash3 = struct;
flash9 = struct;
for sidx = 1:length(lineStr)
    eval(['load' lineStr{sidx} ''])
    
    loadLayerMasks([x.MIP])
    
    for oidx = 1:length(x)
        x(oidx).MIP.ROI = [];
        x(oidx).MIP.ROI(1).mask = x(oidx).MIP.cellMask;
        x(oidx).MIP.ROI(1).color = [1 1 1];        
    end
    
    
    %{
    % % polarized UV flash responses
    flash3.(lineStr{sidx}) = plotFlashResp(x,1);
    flash3.(lineStr{sidx}) = rmfield(flash3.(lineStr{sidx}), {'fig';'ax'});
    close(gcf)
    %}
    
    % blue flash responses
    flash9.(lineStr{sidx}) = plotFlashRespExp9(x,1);
    flash9.(lineStr{sidx}) = rmfield(flash9.(lineStr{sidx}), {'fig';'ax'});
    close(gcf)
    
end
% save(TuBu_flash_data_path,'flash3','flash9')
save(TuBu_flash_data_path,'flash9')
