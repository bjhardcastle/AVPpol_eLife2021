function makePolBarComparison(obj)



% figure,
% for n = 1:6
% fields = [];
% fields.TrialSeqNum = 5;
% fields.TrialPatNum = 2*n - 1;
% X = getActivityFrame(g4,fields);
% subplot(1,6,n)
% imagesc(X)
% end
% subplot(1,6,6)
% n =6;
% fields = [];
% fields.TrialSeqNum = 4;
% X = getActivityFrame(g4,fields);
% subplot(1,6,n)
% imagesc(X)
% colormap(hot)
%% g
% fields = [];
% fields.TrialSeqNum = 6;
% X = getActivityFrame(g,fields);
% subplot(1,2,1)
% imagesc(X)
% title('OF')

% OF
f= figure('color',[0 0 0]);

fields = [];
fields.TrialSeqNum = 6;
X6 = imrotate(getActivityFrameMax(obj,fields),90);
imagesc(X6)
colormap(hot)
axis equal
axis off
set(gca,'position',[0 0 1 1]);
f.Position = get(0, 'Screensize');

set(f, 'InvertHardCopy', 'off');
print(f, fullfile(obj.Folder,['heat6']) , '-dpng','-r450');
close(f)



% BAR
f= figure('color',[0 0 0]);

fields = [];
fields.TrialSeqNum = 5;
X5 = imrotate(getActivityFrameMax(obj,fields),90);
imagesc(X5)
colormap(hot)
axis equal
axis off
set(gca,'position',[0 0 1 1]);
f.Position = get(0, 'Screensize');

set(f, 'InvertHardCopy', 'off');
print(f, fullfile(obj.Folder,['heat5']) , '-dpng','-r450');
close(f)

% POL
f= figure('color',[0 0 0]);

fields = [];
fields.TrialSeqNum = 4;
X4 = imrotate(getActivityFrameMax(obj,fields),90);
imagesc(X4)
colormap(hot)
axis equal
axis off
set(gca,'position',[0 0 1 1]);
f.Position = get(0, 'Screensize');

set(f, 'InvertHardCopy', 'off');
print(f, fullfile(obj.Folder,['heat4']) , '-dpng','-r450');
close(f)

f= figure('color',[0 0 0]);
imshowpair(X5,X4,'ColorChannels','red-cyan','scaling','independent')
f= gcf;
axis equal
axis off
set(gca,'position',[0 0 1 1]);
f.Position = get(0, 'Screensize');

set(f, 'InvertHardCopy', 'off');
print(f, fullfile(obj.Folder,['diff4-5']) , '-dpng','-r450');
close(f)

f= figure('color',[0 0 0]);
imshowpair(X6,X4,'ColorChannels','red-cyan','scaling','independent')
f= gcf;
axis equal
axis off
set(gca,'position',[0 0 1 1]);
f.Position = get(0, 'Screensize');

set(f, 'InvertHardCopy', 'off');
print(f, fullfile(obj.Folder,['diff4-6']) , '-dpng','-r450');
close(f)

f= figure('color',[0 0 0]);
imshowpair(X6,X5,'ColorChannels','red-cyan','scaling','independent')
f= gcf;
axis equal
axis off
set(gca,'position',[0 0 1 1]);
f.Position = get(0, 'Screensize');

set(f, 'InvertHardCopy', 'off');
print(f, fullfile(obj.Folder,['diff5-6']) , '-dpng','-r450');
close(f)


