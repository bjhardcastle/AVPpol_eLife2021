function plotColorPol(obj)
assert(size(obj,1) == 1,'run on single obj only')

if isempty(obj.Frames)
    getFrames(obj)
end

% Make crop box by thresholding activity across all frames 

activityVal = nanmean(obj.Frames(:)) ;
% activityVal = nanmean(obj.AverageFrame(:)) + 0.5*nanstd(obj.AverageFrame(:)) ;
activityBinaryImg = obj.AverageFrame > activityVal;
% figure,imagesc(activityBinaryImg)

% Dectect areas of activity
bw = regionprops(activityBinaryImg);
% Sort by size
[~,sortIdx] = sort([bw.Area],'descend');
% Find the bounding box of the largest detected region
cropBox = bw(sortIdx(1)).BoundingBox;


%%

x1 = getPolAngleActivityFrame(obj,30);
x15 = getPolAngleActivityFrame(obj,60);
x2 = getPolAngleActivityFrame(obj,90);
x25 = getPolAngleActivityFrame(obj,120);
x3 = getPolAngleActivityFrame(obj,150);
x35 = getPolAngleActivityFrame(obj,180);

%% inverted, white bkgrnd

whiteThresh = 0;
x1scaled = x1.*(x1>0)./max(x1(:));
x1scaled = x1scaled.*(x1scaled>whiteThresh);
X1 = cat(3,zeros(size(x1)),x1scaled,zeros(size(x1))); 
% X1 = imcomplement(X1);

x15scaled = x15.*(x15>0)./max(x15(:));
x15scaled = x15scaled.*(x15scaled>whiteThresh);
X15 = cat(3,zeros(size(x15)),x15scaled,zeros(size(x15))); % Green > Purple
% X15 = imcomplement(X15);
x15scaled = x15scaled.*(x15scaled>0.1);

x2scaled = x2.*(x2>0)./max(x2(:));
x2scaled = x2scaled.*(x2scaled>whiteThresh);
X2 = cat(3,x2scaled,zeros(size(x2)),x2scaled); % Purple > Green
% X2 = imcomplement(X2);

x25scaled = x25.*(x25>0)./max(x25(:));
x25scaled = x25scaled.*(x25scaled>whiteThresh);

X25 = cat(3,x25scaled,zeros(size(x25)),x25scaled); % Purple > Green
% X25 = imcomplement(X25);


x3scaled = x3.*(x3>0)./max(x3(:));
x3scaled = x3scaled.*(x3scaled>whiteThresh);
X3 = cat(3,x3scaled,x3scaled,zeros(size(x3))); % Yellow > blue
% X3 = imcomplement(X3);

x35scaled = x35.*(x35>0)./max(x35(:));
x35scaled = x35scaled.*(x35scaled>whiteThresh);

X35 = cat(3,x35scaled,x35scaled,zeros(size(x35))); % Yellow > BLue
% X35 = imcomplement(X35);

X12 = imfuse(X1,X2,'blend','Scaling','joint');
X123 = imfuse(X12,X3,'blend','Scaling','joint');

X12b = imfuse(X15,X25,'blend','Scaling','joint');
X123b = imfuse(X12b,X35,'blend','Scaling','joint');

X123ab = imfuse(X123,X123b,'blend','Scaling','joint');

invFactor = 4;

cropName{1} = 'nocrop';
cropName{2} = 'cropped';
for n = 1:2


f = figure('color',[1 1 1]);
frame = imcomplement(invFactor*X123);
% frame = (2*invFactor*X123);

try
    if n == 2

frame = imcrop(frame,cropBox);
    end
catch
end
imshow(imrotate(frame,90))
% title('30/90/150','color','k')
axis equal
axis off
ax = gca;
ax.Position = [0 0 1 1];
f.Position(3) = 6*size(frame,1);
f.Position(4) = 6*size(frame,2);
    set(f, 'InvertHardCopy', 'off');
    print(f, fullfile(obj.Folder,['polColor_w_30-90-150',cropName{n}]) , '-dpng','-r450');
close(f)

f = figure('color',[1 1 1]);
frame = imcomplement(invFactor*X123b);
try
    if n == 2

frame = imcrop(frame,cropBox);
    end
catch
end
imshow(imrotate(frame,90))
% title('60/120/180','color','k')
axis equal
axis off
ax = gca;
ax.Position = [0 0 1 1];
f.Position(3) = 6*size(frame,1);
f.Position(4) = 6*size(frame,2);

set(f, 'InvertHardCopy', 'off');
print(f, fullfile(obj.Folder,['polColor_w_60-120-180',cropName{n}]) , '-dpng','-r450');
close(f)

f = figure('color',[1 1 1]);
frame = imcomplement(invFactor*X123ab);
try
 if n == 2   
frame = imcrop(frame,cropBox);
 end
catch
end
imshow(imrotate(frame,90))
% title('blended','color','k')
axis equal
axis off
ax = gca;
ax.Position = [0 0 1 1];
f.Position(3) = 6*size(frame,1);
f.Position(4) = 6*size(frame,2);

    set(f, 'InvertHardCopy', 'off');
    print(f, fullfile(obj.Folder,['polColor_w_combined_30+60_etc_',cropName{n}]) , '-dpng','-r450');
close(f)

end
% addExportFigToolbar(gcf)

%% non-inverted, blk bkgrnd

whiteThresh = 0;
x1scaled = x1.*(x1>0)./max(x1(:));
x1scaled = x1scaled.*(x1scaled>whiteThresh);
X1 = cat(3,x1scaled,zeros(size(x1)),x1scaled); 
% X1 = imcomplement(X1);

x15scaled = x15.*(x15>0)./max(x15(:));
x15scaled = x15scaled.*(x15scaled>whiteThresh);
X15 = cat(3,x15scaled,zeros(size(x15)),x15scaled); % Purple
% X15 = imcomplement(X15);
x15scaled = x15scaled.*(x15scaled>0.1);

x2scaled = x2.*(x2>0)./max(x2(:));
x2scaled = x2scaled.*(x2scaled>whiteThresh);
X2 = cat(3,zeros(size(x2)),x2scaled,zeros(size(x2))); % Green
% X2 = imcomplement(X2);

x25scaled = x25.*(x25>0)./max(x25(:));
x25scaled = x25scaled.*(x25scaled>whiteThresh);

X25 = cat(3,zeros(size(x25)),x25scaled,zeros(size(x25))); %  Green
% X25 = imcomplement(X25);


x3scaled = x3.*(x3>0)./max(x3(:));
x3scaled = x3scaled.*(x3scaled>whiteThresh);
X3 = cat(3,zeros(size(x3)),zeros(size(x3)),x3scaled); %  blue
% X3 = imcomplement(X3);

x35scaled = x35.*(x35>0)./max(x35(:));
x35scaled = x35scaled.*(x35scaled>whiteThresh);

X35 = cat(3,zeros(size(x35)),zeros(size(x35)),x35scaled); % BLue
% X35 = imcomplement(X35);

X12 = imfuse(X1,X2,'blend','Scaling','joint');
X123 = imfuse(X12,X3,'blend','Scaling','joint');

X12b = imfuse(X15,X25,'blend','Scaling','joint');
X123b = imfuse(X12b,X35,'blend','Scaling','joint');

X123ab = imfuse(X123,X123b,'blend','Scaling','joint');

invFactor = 6;

cropName{1} = 'nocrop';
cropName{2} = 'cropped';
for n = 1:2


f = figure('color',[1 1 1]);
frame = (invFactor*X123);
% frame = (2*invFactor*X123);

try
    if n == 2

frame = imcrop(frame,cropBox);
    end
catch
end
imshow(imrotate(frame,90))
% title('30/90/150','color','k')
axis equal
axis off
ax = gca;
ax.Position = [0 0 1 1];
f.Position(3) = 6*size(frame,1);
f.Position(4) = 6*size(frame,2);
    set(f, 'InvertHardCopy', 'off');
    print(f, fullfile(obj.Folder,['polColor_b_30-90-150',cropName{n}]) , '-dpng','-r450');
 close(f)

f = figure('color',[1 1 1]);
frame = (invFactor*X123b);
try
    if n == 2

frame = imcrop(frame,cropBox);
    end
catch
end
imshow(imrotate(frame,90))
% title('60/120/180','color','k')
axis equal
axis off
ax = gca;
ax.Position = [0 0 1 1];
f.Position(3) = 6*size(frame,1);
f.Position(4) = 6*size(frame,2);

set(f, 'InvertHardCopy', 'off');
print(f, fullfile(obj.Folder,['polColor_b_60-120-180',cropName{n}]) , '-dpng','-r450');
close(f)

f = figure('color',[1 1 1]);
frame = (invFactor*X123ab);
try
 if n == 2   
frame = imcrop(frame,cropBox);
 end
catch
end
imshow(imrotate(frame,90))
% title('blended','color','k')
axis equal
axis off
ax = gca;
ax.Position = [0 0 1 1];
f.Position(3) = 6*size(frame,1);
f.Position(4) = 6*size(frame,2);

    set(f, 'InvertHardCopy', 'off');
    print(f, fullfile(obj.Folder,['polColor_b_combined_30+60_etc_',cropName{n}]) , '-dpng','-r450');
close(f)

end
% addExportFigToolbar(gcf)



