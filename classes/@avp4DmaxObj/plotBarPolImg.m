function plotBarPolImg(obj,noMask,noiseFilter)
%plotBarPolImg Plot std of exp4/exp5 in polar(blue/red)cols
% [] = plotBarPolImg(obj,noMask,noiseFilter)

if  ~any(obj.TrialSeqNum == 4)
    disp('No pol exp exists')
    return
elseif  ~any(obj.TrialSeqNum == 5) 
    disp('No bar exp exists')
    return
end

if nargin<2 || isempty(noMask)
    noMask = 0;
end
if nargin<3 || isempty(noiseFilter)
    noiseFilter = 1;
end

titleStr = [obj.File];
f = figure('Name',titleStr,'color','w');
ax(1) = subplot(1,1,1);

loadLayerMasks(obj)
if isempty(obj.AverageFrame)
    mspState = obj.UseMSP;
    unattendedState = obj.Unattended;
    
    obj.UseMSP = 0;
    obj.Unattended = 1;
    
    getFrames(obj)
    
    obj.UseMSP = mspState;
    obj.Unattended = unattendedState;
    
    obj.Frames = [];
end
if isempty(obj.InactivityFrame)
    getInactivityFrame(obj)
end
    im = obj.InactivityFrame;
lim_bright = max(im(obj.cellMask));
lim_dark = mean(im(obj.bkgMask));
% normalize:
avgFrame = im - lim_dark;
avgFrame = avgFrame./(mean(avgFrame(obj.layerMask.mask(:))) + 3*std(avgFrame(obj.layerMask.mask(:))) );
avgFrame = 1-avgFrame;
% avgFrame = 1-( (im - lim_dark) ./ max(im(:) - lim_bright - lim_dark )  );
avgFrame(avgFrame>1) = 1;
avgFrame(avgFrame<0) = 0;

for expNum = [4,5]
    e1 = obj.expStart(expNum);
    e2 = obj.expStop(expNum);
   if isempty(obj.Frames)
       getFrames(obj)
   end
   framerange = obj.Frames(:,:,e1:e2);
       
   img = std(framerange,[],3);
   img = (img - mean(img(:)))/std(img(:));

   expimg{expNum-3} = img;  
end
  
    polimg = expimg{1};
    polimg = polimg.*(polimg>0);
    barimg = expimg{2};
    barimg = barimg.*(barimg>0);
    barimg = barimg./max(barimg(:));
    polimg = polimg./max(polimg(:));
    combimg = -barimg + polimg;

    
    if noiseFilter

        s{2} = [0 1 0; 0 1 0; 0 0 0];
        s{3} = [0 0 0; 0 1 0; 0 1 0];
        s{4} = [0 0 0; 0 1 1; 0 0 0];
        s{5} = [0 0 0; 1 1 0; 0 0 0];
        s{6} = [1 0 0; 0 1 0; 0 0 0];
        s{7} = [0 0 1; 0 1 0; 0 0 0];
        s{8} = [0 0 0; 0 1 0; 1 0 0];
        s{9} = [0 0 0; 0 1 0; 0 0 1];
        s{end+1} = [0 0 0; 0 1 0; 0 0 0];
        
        for n = 2
            for sIdx = 2:length(s)
                
                % Hit or miss removal of isolated single pixels in center of 3x3 sq
                notNan = ~isnan(combimg);
                a = imerode(notNan,s{sIdx});
                b = imerode(~notNan,~s{sIdx});
                
                no_islands = (a&b);
                combimg(no_islands) = nan;
                
           
            end
        end
    end
    
    if ~noMask
       combimg(~obj.layerMask.mask) = nan;
        combFrame = zeros(size(avgFrame)    );
    combFrame(isnan(combimg)) = -180.*avgFrame(isnan(combimg));
        combFrame(~isnan(combimg)) = 90.*combimg(~isnan(combimg)) + 90;

    layer_im = imagesc(imrotate(combFrame,90), 'Parent', ax(1));
    colormap(ax(1),[flipud(gray(181));[0 0 0];polarmap(180)])
    set(ax(1),'CLim',[-180 180])
    
else
    layer_im = imagesc(imrotate(combimg,90), 'Parent', ax(1));

    polarmap
    set(ax(1),'CLim',[-1 1])
end


addScalebar(obj,ax(1),10)
axis(ax(1),'image')
axis(ax(1),'off')
getAVPplotParams
setAVPaxes(ax(1),defaultImageHeight_cm)
tightfig(f)
addExportFigToolbar(f)