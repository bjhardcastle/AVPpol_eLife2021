function superTimes(obj)
for oidx = 1:length(obj)
   h = str2double(obj(oidx).TimeStr(1:2));
   m =str2double(obj(oidx).TimeStr(3:4));
   d = m/60;
   
   t(oidx) = h+d;
   
end

figure
scatter(t,zeros(size(t)),20,'filled')
xlim([6 20])
