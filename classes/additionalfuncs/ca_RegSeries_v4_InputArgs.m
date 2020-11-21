function ca_RegSeries_v4_InputArgs(im0, path0, poolSiz, walkThresh, usfac)
%CA_REGSERIES_V4_INPUTARGS Performs image registration on a specified Tiff file.
% ca_RegSeries_v4_InputArgs(im0,path0,poolSiz, walkThresh, usfac)
% 
% output: the registered version of the input file, saved with suffix
% '_reg' before the file extension.
%
% v_4 was created so that this modified version with input arguments could 
% be found by SlidebookObjs (along with sint_SaveStack_v2 below)
% It also makes use of Matlab's recent addition of Tifflib class for faster
% read/write speeds. Unnecessary output files also suppressed.
%
% Tiff stack is registered in two stages:
%
% in the first stage (loop 1), a "poolSiz" number of consecutive frames are averaged
% and a condensed stack of the series, poolStk, is generated from these
% pooled frames. consecutive frames of the poolStk are iteratively registered
% to each other using fft based sub-pixel registration (loop 2).
%
% in the second stage (loop 3), individual frames from the full scale series are
% registered to the pooled average frames that they contributed to, again
% with fft based sub-pixel registration. the xy displacements needed to bring
% the full series into register are the sum of the pooled stack displacements
% and the single frame displacements.
%
% the registered stack is generated at the end (loop 4) with appropriate
% zeros-padding to accommodate movement generated due to the displacements.
%
% there are three built-in parameters that can be varied within the script:
% 1. poolSiz: this can be changed based on empirical experience. for low
% time resolution, long imaging sessions, poolSiz set to 12 minutes of
% imaging works well.
% 2. walkThesh: frame-to-frame displacement magnitude, in pixels, considered a
% 'jerk' that is likely to be an error. if a jerk is detected, the frame is
% filtered to reduce intensity variation and registered again. if the jerk
% persists, the frame is 'skipped'; that is, it is not displaced relative
% to the preceeding frame. reducing this value should not have a negative
% impact on the quality of the registration; the time taken to process the
% data would increase, however.
% 3. usfac: up-sampling factor used in sub-pixel registration. increasing
% this may improve the quality but will also increase run time. 
%
% See also runTifReg.
%%
if nargin < 5
    usfac = 10; % up-sampling factor for sub-pixel registration.
end
if nargin < 4
    walkThresh = 5; % hard threshold to identify big jerks that need extra attention
end
if nargin <3
    poolSiz = 400;
end 

im1 = [im0(1:end-4) '_reg.tif'];

% New method using tifflib
w = warning('off','all');
tifobj = Tiff([path0 im0],'r+'); % Open object in read mode

% First get the number of frames:
setDirectory(tifobj,1)          % Make sure we're at the first frame.
while ~lastDirectory(tifobj)
    nextDirectory(tifobj);
end
objframes = currentDirectory(tifobj);

% If PoolSiz is larger than number of frames in tif, reduce it to match
if poolSiz > objframes
    poolSiz = objframes;
end

% Initialize inStk:
inStk = zeros( getTag(tifobj,'ImageLength'), getTag(tifobj,'ImageWidth'), objframes , 1);

% Now read tiff image data from each directory:
setDirectory(tifobj,1)  % Make sure we're at the first frame (default),
for n = 1:objframes
    inStk(:,:,n) = read(tifobj);
    if n ~= objframes
        nextDirectory(tifobj);
    end
end

close(tifobj); clear tifobj

% Set warnings back to previous status
warning(w);

stkX = size(inStk,2);
stkY = size(inStk,1);
stkZ = size(inStk,3);
    
inStk(end,:,:) = 0; % clear flyback line

% variables for fft-based translation of the images. code lifted from
% dftregistration demo.
nr = stkY; nc = stkX;
Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
[Nc,Nr] = meshgrid(Nc,Nr);
phase = 0;

% calculating the first and last frames (i.e. loopLims) of each pool. the
% last pool may be larger than poolSiz, but smaller than 2*poolSize.
poolNum = floor(stkZ/poolSiz);
poolLims = zeros(poolNum,2);
poolLims(:,1) = 1:poolSiz:poolSiz*poolNum';
poolLims(:,2) = poolLims(:,1)+poolSiz-1;
poolLims(end,2) = stkZ;
   
poolStk = zeros(stkY,stkX,poolNum,'single');

for p0 = 1:poolNum % loop 1: generating the poolStk
    poolStk(:,:,p0) = mean(inStk(:,:,poolLims(p0,1):poolLims(p0,2)),3);
end

y_allP = zeros(poolNum,1); % frame-to-frame y displacement 
x_allP = zeros(poolNum,1); % frame-to-frame x displacement 
walkP_x = 0; % cumulative y displacement; loop specific variable 
walkP_y = 0; % cumulative x displacement; loop specific variable  

poolStkR = poolStk;

disp(['registering pooled series...']);

for p0 = 2:poolNum % loop 2: registering the poolStk

    slc0 = poolStk(:,:,p0-1);
    slc1 = poolStk(:,:,p0);
    output = dftregistration(fft2(slc0),fft2(slc1),usfac);
    y_offset = round(usfac*output(3))/usfac; % limiting the precision to 1/usfac pix
    x_offset = round(usfac*output(4))/usfac;
    
    y_allP(p0) = y_offset;
    x_allP(p0) = x_offset;
    
    walkP_y = walkP_y+y_offset;
    walkP_x = walkP_x+x_offset;
   
    % using fft to apply the translations. code lifted from dftregistration
    % demo.
    deltar = -walkP_y;
    deltac = -walkP_x;
    g = ifft2(fft2(slc1).*exp(1i*2*pi*(deltar*Nr/nr+deltac*Nc/nc))).*exp(-1i*phase);
    g0 = real(g);
    
    % zero-flooding the shifts, which are wrapped around with the
    % fft-based translation
    if walkP_y<0
        g0(end-ceil(abs(walkP_y))+1:end,:) = 0;
    elseif walkP_y>0
        g0(1:ceil(abs(walkP_y)),:) = 0;
    end
    if walkP_x<0
        g0(:,end-ceil(abs(walkP_x))+1:end) = 0;
    elseif walkP_x>0
        g0(:,1:ceil(abs(walkP_x))) = 0;
    end
    
    poolStkR(:,:,p0) = g0;

end

% optional: can save pre- and post-registration pooled Stk. comment out if
% not necessary.
% sint_SaveStack_v2(poolStk,16,path0,'poolStk.tif',tifExt);
% sint_SaveStack_v2(poolStkR,16,path0,'poolStk_reg.tif',tifExt);
%%   
y_all = zeros(stkZ,1);
x_all = zeros(stkZ,1);

y_allP_sum = cumsum(y_allP); % cumulative y displacement from loop 2. 
x_allP_sum = cumsum(x_allP); % cumulative x displacement from loop 2. 

disp(['registering full series to pooled series...']);

for p0 = 1:poolNum % loop 3: registering the full Stk
    slc0 = poolStk(:,:,p0);
    for i0 = poolLims(p0,1):poolLims(p0,2)
        if mod(i0,1000)==0 % progress report
            disp(['...frame ',num2str(i0),', ',num2str(100*(i0/stkZ),'%6.2f'),'% complete,']);
        end
        
        slc1 = inStk(:,:,i0);
        
        output = dftregistration(fft2(slc0),fft2(slc1),usfac);
        
        if norm([output(4) output(3)])>walkThresh % try again with field flattening if there's a big jerk
%             disp(['!: jerk in frame ',num2str(i0)]);
            slc0f = medfilt2(slc0);
            slc0g = imgaussian(slc0f,10);
            mx0 = max(slc0g(:));
            n0 = (slc0g/mx0).^-1;
            
            slc1f = medfilt2(slc1);
            slc1g = imgaussian(slc1f,10);
            mx1 = max(slc1g(:));
            n1 = (slc1g/mx1).^-1;
            
            output = dftregistration(fft2(slc0.*n0),fft2(slc1.*n1),usfac);
            
            if norm([output(4) output(3)])>walkThresh
                output(3) = 0;
                output(4) = 0;
%                 disp('--> frame skipped');
            end
        end
        
        
        y_offset = round(usfac*output(3))/usfac; % limiting the precision to 1/usfac pix
        x_offset = round(usfac*output(4))/usfac;
        
        y_all(i0) = y_offset+y_allP_sum(p0);
        x_all(i0) = x_offset+x_allP_sum(p0);
             
    end

end

padYtop = abs(floor(min(y_all)));
padYbot = abs(ceil(max(y_all)));

padXtop = abs(floor(min(x_all)));
padXbot = abs(ceil(max(x_all)));

inStkP = padarray(inStk,[padYtop padXtop],'pre');
inStkP = padarray(inStkP,[padYbot padXbot],'post');

clear inStk

nrP = size(inStkP,1); ncP = size(inStkP,2);
NrP = ifftshift([-fix(nrP/2):ceil(nrP/2)-1]);
NcP = ifftshift([-fix(ncP/2):ceil(ncP/2)-1]);
[NcP,NrP] = meshgrid(NcP,NrP);
phase = 0;

disp(['applying displacements to full series...']);

for i0 = 2:stkZ % loop 4: applying the registration to full Stk with zero padding
    if mod(i0,1000)==0 % progress report
        disp(['...frame ',num2str(i0),', ',num2str(100*(i0/stkZ),'%6.2f'),'% complete,']);
    end
    
    slc1 = inStkP(:,:,i0);

    walk_y = y_all(i0);
    walk_x = x_all(i0);
        
    % using fft to apply the translations. code lifted from dftregistration
    % demo.
    deltar = -walk_y;
    deltac = -walk_x;
    g = ifft2(fft2(slc1).*exp(1i*2*pi*(deltar*NrP/nrP+deltac*NcP/ncP))).*exp(-1i*phase);
    g0 = real(g);
    
    % zero-flooding the shifts, which are wrapped around with the
    % fft-based translation
    if walk_y<0
        g0(end-ceil(abs(walk_y))+1:end,:) = 0;
    elseif walk_y>0
        g0(1:ceil(abs(walk_y)),:) = 0;
    end
    if walk_x<0
        g0(:,end-ceil(abs(walk_x))+1:end) = 0;
    elseif walk_x>0
        g0(:,1:ceil(abs(walk_x))) = 0;
    end
    
    inStkP(:,:,i0) = g0;

end

% saving the output of second round 
disp(['final loop complete, saving output stack...']);
sint_SaveStack_v2(inStkP,16,path0,im1,im0);

dlmwrite([path0 'regDispXY.txt'],[x_all y_all],'\t');
disp('cumulative displacements saved.');

