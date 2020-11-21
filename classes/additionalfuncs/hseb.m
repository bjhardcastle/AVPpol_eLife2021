function H = hseb(x,y,errBar,lineProps,transparent)
% H = GSEB(x,y,errBar,lineProps,transparent)
%
% HISTOGRAM Shaded Error Bars (MSEB) makes a 2-d plot containing A HISTOGRAM
% with pretty shaded error bars. ben modded 2020
%
% This is an extension of the popular shadedErrorBar, by Rob Campbell,
% enabling plotting of multiple data lines with overlapping errorbars, as 
% well as turning off legends for all elements but the main lines. MSEB is
% directed at using the default renderer instead of openGL, which is known
% to cause problems, e.g., with logarithmic axes and not being able to save
% figures as vector graphics in the eps-format.
% To avoid the error bar patches concealing previously plotted main lines, 
% the different elements are plotted in a suitable order. The default
% setting plot edges of overshadowed patches in an non-obtrusive manner so 
% all error bars can be seen but avoids cluttering of the plot (see the 
% third example on how to customise this). The patches are plotted in
% reverse order to make sure the first entry is "on top".
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, Cognitive systems - January 2015
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of C x N, where C is the number of
%     lines to be plotted and N is the number of samples in each line and
%     should be equal to length(X)
% errBar - if a vector we draw symmetric errorbars. If it has a
%          size of [C,length(x),2] then asymmetric error bars are drawn,
%          with row 1 being the upper bar and row 2 being the lower
%          bar. In the present version errBar does not support two function
%          handles.
% lineProps  - [optional. Can also be set as "[]" for default settings].
%              Struct containing fields that define lineproperties for the 
%              plot function. It is possible to only define some of the
%              fields.
%       .col - cell array, where each element defines the colour of each
%              line. This can be done using either strings or three-element
%              RGB vectors e.g., either 'b' or [0 0 1] for blue.
%     .style - linestyle of the lines from y-data. Default is '-'.
%     .width - linewidth of the lines from y-data. Default is 2.
% .edgestyle - linestyle of edges that are overlapped by errorbars from
%              other lines. 
%             
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image. openGL does not support having logarithmic axes.
%
% Outputs
% H - structure with an element for each line entry containing handles to 
%     the generated plot objects (e.g. H(c) contains the handles to the
%     c'th line entry. 
%
%
%% Examples:
%
% x = 1:100; y = randn(1,100,30)*10;
% t = (1:100) - 50; y(2,:,:) = repmat(t,[1,1,30]);
% y(2,:,:) = y(2,:,:) + 0.06.*y(2,:,:).^2 + randn(1,100,30)*10;
% y(3,:,:) = repmat(t,[1,1,30]);
% y(3,:,:) = 60 - abs(y(3,:,:))+ randn(1,100,30)*10;
% y_mean = mean(y,3); y_std = std(y,[],3);
% figure; title('Default renderer'),
% mseb(x,y_mean,y_std);ylim([-50 150])
% legend('Line 1','Line 2','Line 3')
% figure; title('openGL');
% mseb(x,y_mean,y_std,[],1);ylim([-50 150]) 
% legend('Line 1','Line 2','Line 3')
% 
% % openGL issue with logarithmic axes:
% figure; title('Default renderer')
% mseb(x,y_mean,y_std);ylim([-50 150])
% set(gca,'xScale','log')
% figure; title('openGL')
% mseb(x,y_mean,y_std,[],1);ylim([-50 150])
% set(gca,'xScale','log')
%
% % Defing only some of the line parameters
% lineProps.width = 1;
% lineProps.edgestyle = ':';
% figure; title('Custom line properties')
% mseb([],y_mean,y_std,lineProps);ylim([-50 150])
% transparent = 0;

%% Error checking
error(nargchk(3,5,nargin))

% Cheking the y data
[C, N] = size(y);
if N==1
    C = 1;
    N = length(y);
    y = y';
end

% Cheking the x data
if isempty(x)
    x=repmat(1:N,[C,1]);
elseif length(x(:))==N && C>1
    x=repmat(x(:)',[C,1]);
end

if (size(x,1) ~= C) || (size(x,2)~=N)
    error('inputs x and y do not have same dimensions')
end

% ben: convert counts at bin centers into flat lines, like a regular
% histogram plot 
assert(C==1,'Needs work to plot multiple histograms at once')
% x should be bin centers
x(end+1) = x(end)+mode(diff(x));
edges = x-0.5*mode(diff(x));
binedges = reshape([edges;edges],1,[]);
binedges(1) = [];
binedges(end) = [];
y_binned = reshape([y;y],1,[]);

y = y_binned;
x = binedges;


% Checking ErrBar dimensions. If only one error bar is specified then we
% will mirror it, turning it into both upper and lower bars.
% % if (size(errBar,1) ~= C) || (size(errBar,2)~=N)
% %     if size(errBar,1)==N && size(errBar,2)==2 && size(errBar,3)~=1 && C==1
% %         errBar = permute(errBar,[3,1,2]);
% %     elseif size(errBar,1)==N && size(errBar,2)==1 && size(errBar,3)~=1 && C==1
% %         errBar=repmat(errBar,[1,2]);
% %         errBar = permute(errBar,[3,1,2]);
% %     else
% %         error('inputs errBar and y do not have same dimensions')
% %     end
% % end
% % if size(errBar,3)==1
% %     errBar = repmat(errBar,[1,1,2]);
% % end

% ben: make blocky errorbar to match histogram
if ~isempty(errBar)
errBar_binned = reshape([errBar;errBar],1,[]);
errBar = errBar_binned;
else
    errBar = nan(size(y));
end
%% Setting default line options
if nargin<4
    lineProps = [];
end
if nargin<5 || ~isnumeric(transparent)
    transparent=0;
end
if ~isfield(lineProps,'col')
    for c = 1:C
        colours = 'brkgmcy';
        col_ind = rem(c-1,length(colours))+1;
        lineProps.col{c} = colours(col_ind);
    end
end
if ~isfield(lineProps,'style'),	lineProps.style = '-'; end
if ~isfield(lineProps,'width'),	lineProps.width = 2; end
if ~isfield(lineProps,'edgestyle'),	lineProps.edgestyle = '--'; end

%% First loop over the number of lines to be plotted
% (More or less contains the original code from shadedErrorBar)
holdStatus=ishold;
if ~holdStatus, hold on,  end

for c = C:-1:1
    %% Plotting patches
    % Plot the main line. We plot this first in order to extract the RGB values
    % for the line colour. I am not aware of a function that does this.
    H(c).mainLine=plot(x(c,:),y(c,:),'color',lineProps.col{c});
    
    
    % Work out the color of the shaded region and associated lines
    % Using alpha requires the render to be openGL and so you can't
    % save a vector image. We therefore provide the option of choosing alpha
    % or a de-saturated solid colour for the patch surface.
    
    col=get(H(c).mainLine,'color');
    edgeColor=col+(1-col)*0.55;
    patchSaturation=0.8; %How de-saturated or transparent to make the patch
    if transparent
        faceAlpha=patchSaturation;
        patchColor=col;
        set(gcf,'renderer','openGL')
    else
        faceAlpha=1;
        patchColor=col+(1-col)*(1-patchSaturation);
        set(gcf,'renderer','painters')
    end
    
    %Calculate the y values at which we will place the error bars
    uE=y(c,:)+errBar(c,:,1);
%     lE=y(c,:)-errBar(c,:,2);
      lE=y(c,:);  % % ben: positive patch only, above histogram
     
    %Make the cordinats for the patch
    yP=[lE,fliplr(uE)];
    xP=[x(c,:),fliplr(x(c,:))];
    
    %remove any nans otherwise patch won't work
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];
    
    if isempty(errBar) || all(isnan(errBar)) % data were all nans: make an invisible patch to match the invisible line which is plot
        H(c).errPatch=patch(0,0,1,'facecolor','none',...
            'edgecolor','none',...
            'facealpha',0);
        
    else
        
        H(c).errPatch=patch(xP,yP,1,'facecolor',patchColor,...
            'edgecolor','none',...
            'facealpha',faceAlpha);
        
    end
    
    %Make nice edges around the patch.
%     H(c).edge(1)=plot(x(c,:),lE,'-','color',edgeColor);
%     H(c).edge(2)=plot(x(c,:),uE,'-','color',edgeColor);
end

%% Second loop over the number of lines to be plotted
for c = C:-1:1
    %% Plot egdes for patchoverlaps
    col=get(H(c).mainLine,'color');
    edgeColor=col+(1-col)*0.25;
%     lE = get(H(c).edge(1), 'ydata');
%     uE = get(H(c).edge(2), 'ydata');
%     
%     H(c).edgeoverlap(1)=plot(x(c,:),lE,lineProps.edgestyle,'color',edgeColor);
%     H(c).edgeoverlap(2)=plot(x(c,:),uE,lineProps.edgestyle,'color',edgeColor);
    
end

%% Third loop over the number of lines to be plotted
for c = 1:C
    %% Plot mainlines
    %The main line is now covered by the patch object and was plotted first to
    %extract the RGB value of the main plot line. I am not aware of an easy way
    %to change the order of plot elements on the graph so we'll just remove it
    %and put it back (yuk!)
     col=lineProps.col{c};
     mainLineColor = col*0.8;
    delete(H(c).mainLine)
    H(c).mainLine=plot(x(c,:),y(c,:),lineProps.style,'color',mainLineColor,...
        'linewidth',lineProps.width);
    
   
% ben: add patch instead of line

   
% delete(H(c).mainLine)

    % Work out the color of the shaded region and associated lines
    % Using alpha requires the render to be openGL and so you can't
    % save a vector image. We therefore provide the option of choosing alpha
    % or a de-saturated solid colour for the patch surface.
    
    patchSaturation=0.8; %How de-saturated or transparent to make the patch
    if transparent
        faceAlpha=patchSaturation;
        patchColor=col;
        set(gcf,'renderer','openGL')
    else
        faceAlpha=1;
        patchColor=col+(1-col)*(1-patchSaturation);
        set(gcf,'renderer','painters')
    end
    
    %Calculate the y values at which we will place the error bars
    uE=y(c,:);
    lE=0.*y(c,:);  % % ben: assumes histogram counts start at 0 on yaxis
     
    %Make the cordinats for the patch
    yP=[lE,fliplr(uE)];
    xP=[x(c,:),fliplr(x(c,:))];
    
    %remove any nans otherwise patch won't work
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];
    
    if isempty(xP) && isempty(yP) % data were all nans: make an invisible patch to match the invisible line which is plot
        H(c).mainPatch=patch(0,0,1,'facecolor','none',...
            'edgecolor','none',...
            'facealpha',0);
        
    else
        
        H(c).mainPatch=patch(xP,yP,1,'facecolor',patchColor,...
            'edgecolor','none',...
            'facealpha',faceAlpha);
        
    end
    
    %Make nice edges around the patch.
%     H(c).edge(1)=plot(x(c,:),lE,'-','color',edgeColor);
%     H(c).edge(2)=plot(x(c,:),uE,'-','color',edgeColor);    

    
    %% Turn legendinformation off for all but the main lines
    set(get(get(H(c).errPatch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
     set(get(get(H(c).mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
%     set(get(get(H(c).edge(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%     set(get(get(H(c).edge(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%     set(get(get(H(c).edgeoverlap(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
%     set(get(get(H(c).edgeoverlap(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(H(c).mainPatch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    
end
if ~holdStatus, hold off, end

