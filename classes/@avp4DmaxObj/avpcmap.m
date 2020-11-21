function avpcmap(obj)
% Shortcut to apply pol colormap map to tuning map image 
f = gcf;
ax  = gca;
set(f,'color','w');
pbaspect([1,1,1])

colormap([[1 1 1];flipud(hsv)])
set(ax,'CLim',[0 180])
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

if isempty(obj.ZPlane)
    titleSuffix = ' MIP';
else
    titleSuffix = [' layer' num2str(obj.ZPlane)];
end
ax.Title.String = [obj.Line obj.DateStr ' fly' num2str(obj.fly) titleSuffix];

end