figure(1)
triplot(DT,'k')
xlim([min(x)-2, max(x)+2])
ylim([min(y)-2, max(y)+2])
hold on
rowx = min(x)-1:0.25:max(x)+ 1;
rowy = min(y)*ones(size(rowx));
plot(rowx,rowy,'k','LineWidth',3')
hold on
rowx = min(x)-1:0.25:max(x)+ 1;
rowy = max(y)*ones(size(rowx));
plot(rowx,rowy,'k','LineWidth',3')
axis off
print(gcf,'triangular.png','-dpng','-r1000')