figure(2)
h = plot(Strain, abs(Sigma(:,2)));
set(h, 'LineWidth', 1.5)
xlabel('\delta \rightarrow','Fontsize',16,'FontName',...
    'Arial','FontWeight', 'bold'), 
ylabel('Fy \rightarrow','Fontsize',16,'FontName', ...
    'Arial','FontWeight', 'bold'),
set(gca,'FontWeight', 'bold');
dim = [.65 .5 .3 .3];
str = ['200X200 Lattice'];
t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
t.FontSize = 12;
dim = [.65 .4 .3 .3];
str = ['p = 30%'];
t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
t.FontSize = 12;
hold on
h = plot(StrainExp, StressExp,'LineWidth', 1.5);
% h = plot(CompData(:,1),CompData(:,2));
set(h, 'LineWidth', 1.5)
legend({'Numerical','Experiment'}, ...
    'Location','northeast','Fontsize',10),



% L = 25;
% A = L^2;
% StressExp = CompData(:,2)/A;
% StrainExp = CompData(:,1)/L;
% h = plot(StrainExp, StressExp,'LineWidth', 1.5);

Eav(1) = (StressExp(150)-StressExp(100))/(StrainExp(150)-StrainExp(100))
Eav(2) = (StressExp(500)-StressExp(400))/(StrainExp(500)-StrainExp(400))
Eav(3) = (StressExp(1000)-StressExp(900))/(StrainExp(1000)-StrainExp(900))
% Eef = mean(Eav)
% 
E = (Sigma(10,2) - Sigma(5,2))/(Strain(10)-Strain(5))

% Eav(1) = (CompData(150,2)-CompData(100,2))/(CompData(150,1)-CompData(100,1))
% Eav(2) = (CompData(500,2)-CompData(400,2))/(CompData(500,1)-CompData(400,1))
% Eav(3) = (CompData(1000,2)-CompData(900,2))/(CompData(1000,1)-CompData(900,1))
% % Eef = mean(Eav)
% % 
% % E = (Force(16,2) - Force(5,2))/(delta(16)-delta(5))


%%
