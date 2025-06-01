%% Extracting stress-strain data for each porosity
% clear all
% load('250.mat')
% 
% figure(2)
% h = plot(Strain, abs(Sigma(:,2)));
% set(h, 'LineWidth', 1.5)
% xlabel('\delta \rightarrow','Fontsize',16,'FontName',...
%     'Arial','FontWeight', 'bold'), 
% ylabel('Fy \rightarrow','Fontsize',16,'FontName', ...
%     'Arial','FontWeight', 'bold'),
% set(gca,'FontWeight', 'bold');
% dim = [.65 .5 .3 .3];
% str = ['200X200 Lattice'];
% t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
% t.FontSize = 12;
% dim = [.65 .4 .3 .3];
% str = ['p = 30%'];
% t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
% t.FontSize = 12;
% 
% load('p40.mat')
% i = 3;
% x(:,i) = Strain(1:200);
% y(:,i) = abs(Sigma(1:200,2));
% save('p40.mat','x','y')
% 

%% Storing Every data in one variable
% F = zeros(500,3,7);
% d = zeros(500,3,7);

% load('p70.mat')
% i = 7;
% 
% F(:,1,i) = y(:,1);
% F(:,2,i) = y(:,2);
% F(:,3,i) = y(:,3);
% 
% 
% d(:,1,i) = x(:,1);
% d(:,2,i) = x(:,2);
% d(:,3,i) = x(:,3);
% 
% save('StressStrainData.mat','d','F')

%% Extracting stress field data at different step
clear all
load('300.mat')
NodalStress
FieldVisualisation
% SF = zeros(Ny,Nx,6);
load('SF_30R1.mat')
i = 6;
SF(:,:,i) = Fieldq;

save('SF_30R1.mat','SF')

% At the end
x = X;
y = Y;
SigO = zeros(1,10);
for i = 1:6
    SigO(i) = abs(Sigma(i*50,2));
end
save('SF_30R1.mat','x','y','SigO','-append')




