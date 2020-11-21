function axH = pbgridplot(tuning,psi,datetimestr,sortpsi,pubangles,varargin)
assert(nargin>0 && ~isempty(tuning),'Must provide array of tuning vectors')
% tuning should be [n x 8] array (8 PB glomeruli pairs)
if size(tuning,1)==8 && size(tuning,2)~=8
    tuning = tuning';
end

if nargin<2 || isempty(psi)
    % no psi speficified, all vectors will be unit length
    psi = ones(size(tuning));
end
if nargin < 5 || isempty(pubangles)
    % guess if pub angles used
    if all(tuning(:)<=180) && all(tuning(:)>=0)
        pubangles = 0;
    elseif all(tuning(:)<=90) && all(tuning(:)>=-90)
        pubangles = 1;
    else
        error('not sure which angle format was used')
    end
    
end

if pubangles
    % express in original angle format
    tuning = wrapTo360(-tuning-270);
end

if nargin < 3 || isempty(datetimestr)
    datetimestr = cellstr(num2str([1:size(tuning,1)]'));
elseif ~iscell(datetimestr)
    datetimestr = cellstr(datetimestr);
end

if nargin < 4 || isempty(sortpsi)
    sortpsi = 0;
end

if sortpsi
    % sort according to mean psi. plot vertically, highest will be at the
    % bottom
    [~,sortidx] = sort(mean(psi,2),'descend');
    psi = psi(sortidx,:);
    tuning = tuning(sortidx,:);
    datetimestr = datetimestr(sortidx);
end

figure('color','w')
axH = gca;
hold(axH,'on')
getAVPplotParams

numflies = size(tuning,1);
xlim([0 8])
ylim([0 numflies])
axH.XTick = [0:8];
axH.YTick = [0:numflies];
axH.XTickLabelRotation = 90;
axH.XAxis.Label.String = 'paired glomeruli';
axH.XAxis.Visible = 'on';
axH.XTickLabelRotation = 45;
axH.XTick = [0.5:1:7.5];
axH.XTickLabels = {'8L/2R';'7L/3R';'6L/4R';'5L/5R';'4L/6R';'3L/7R';'2L/8R';'1L/1R'};
% ax.XTickLabels = {'8L/2R';'';'';'5L/5R';'';'';'';'1L/1R'};
daspect(axH,[1,1,1])
% axH.Color = 'none';

% setup grid of centers for vectors:
xcenters = linspace(0,8,17);
xcenters = xcenters(2:2:end);
ycenters = linspace(0,numflies,numflies*2+1);
ycenters = ycenters(2:2:end);



for y = 1:numflies
    for x = 1:8
        p = plotPSIpatch(axH,xcenters(x),ycenters(y),psi(y,x));
        
        h = plotVector(axH,xcenters(x),ycenters(y),tuning(y,x),psi(y,x));
        h.Color = 'k';
        h.LineWidth = 0.5;     
    end
    
    s = plotCircCorrStar(axH,xcenters(x),ycenters(y),tuning(y,:),psi(y,:),numflies);
    s.FontSize = axisLabelFontSize;
end

axH.YTick = ycenters;
axH.YTickLabels = {repmat(' ',numflies,1)};

% axH.YTickLabels{1} = 'fly 1';
% axH.YTickLabels{end} = ['fly ' num2str(numflies)];

for yt= 1:length(datetimestr)
    axH.YTickLabels{yt} = ['fly ' num2str(datetimestr{yt})];
end

setAVPaxes(axH,[],defaultAxisHeight_cm)
tightfig(gcf)
axH.Layer = 'top';
uistack(axH,'top');
addExportFigToolbar(gcf)


end

function p = plotPSIpatch(ax,xc,yc,psi)
colors = [1 1 1;flipud(magma(255))];
psicol = colors(round(255*psi+1),:);
x = [xc-0.5  xc+0.5 xc+0.5 xc-0.5];
y = [yc-0.5 yc-0.5 yc+0.5 yc+0.5];
p = patch(ax,x,y,psicol,'EdgeColor','none');
uistack(p,'bottom');
end

function s = plotCircCorrStar(ax,xc,yc,tuning,psi,numflies)
% apply bonferroni correction
% alpha = 0.05/numflies;

alpha = 0.05;

% get glomeruli positions as a circular variable, [0:2pi]
pos = linspace(0,2*pi,17);
pos = pos(2:2:end);

% double tuning values, since they're axial
tun = deg2rad(2.*tuning);

% get lin-circ correlation
w=ones(1,8);
[rho, pval1] = circlin_wcorr(tun,pos,[-2 2],w');

% text will display rho value
txtstr1 = num2str(rho,'%.2f');

% add star to text if it passes significance
if pval1 < alpha
    txtstr1 = [txtstr1 ' *'];
end



% get circ-circ correlation
% may not need to double angles here
tun = deg2rad(tuning);
pos = linspace(0,2*pi,17);
pos = pos(2:2:end);

[rho,pval2] = circ_corrcc(tun,pos);
% [rho,pval] = circ_wcorrcc(tun,pos,w');
    
% text will display rho value
txtstr2 = num2str(rho,'%.2f');

% add star to text if it passes significance
if pval2 < alpha
    txtstr2 = [txtstr2 ' *'];
end

%{
% add text to plot
s = text(ax,xc+1,yc,[txtstr1 ' / ' txtstr2]);
%}
%only add circ corr, and if significant:
% add text to plot
if pval1 >= alpha
    txtstr1 = '';
end
if pval2 >= alpha
    txtstr2 = '';
end
s = text(ax,xc+1,yc,[txtstr1 '/' txtstr2]);


end

function h = plotVector(ax,xc,yc,angle,length)
% use original angles, not published.
% [0:180]
% 0 deg is horizontal
angle = 0.5.*wrapTo360(2.*angle);
% get x-component
x = length*cosd(angle);
y = length*sind(angle);

x1 = xc - x;
y1 = yc + y;
x2 = xc + x;
y2 = yc - y;

h = plot(ax,[x1,x2],[y1,y2]);

end

function [rho pval] = circ_wcorrcc(alpha1, alpha2, W)
%
% [rho pval ts] = circ_corrcc(alpha1, alpha2)
%   Circular correlation coefficient for two circular random variables.
%
%   Input:
%     alpha1	sample of angles in radians
%     alpha2	sample of angles in radians
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%
% References:
%   Topics in circular statistics, S.R. Jammalamadaka et al., p. 176
%
% PHB 6/7/2008
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
W = W./sum(W);

if size(alpha1,2) > size(alpha1,1)
	alpha1 = alpha1';
end

if size(alpha2,2) > size(alpha2,1)
	alpha2 = alpha2';
end

if size(W,2) > size(W,1)
	W = W';
end

if length(alpha1)~=length(alpha2)
  error('Input dimensions do not match.')
end

% compute mean directions
n = length(alpha1);
alpha1_bar = circ_mean(alpha1,W);
alpha2_bar = circ_mean(alpha2,W);

% compute correlation coeffcient from p. 176
num = sum(W.*sin(alpha1 - alpha1_bar) .* sin(alpha2 - alpha2_bar))./sum(W);
den = sqrt( (sum(W.*sin(alpha1 - alpha1_bar).^2)./sum(W)) .* (sum(W.*sin(alpha2 - alpha2_bar).^2)./sum(W)) );
rho = num / den;	

% compute pvalue
l20 = mean(W.*sin(alpha1 - alpha1_bar).^2);
l02 = mean(W.*sin(alpha2 - alpha2_bar).^2);
l22 = mean(W.*(sin(alpha1 - alpha1_bar).^2) .* W.*(sin(alpha2 - alpha2_bar).^2));

ts = sqrt((n * l20 * l02)/l22) * rho;
pval = 2 * (1 - normcdf(abs(ts)));

end
