%% Program to estimate the following:
%% Number of bonds with applied displacement increment
%% Energy released with applied loading

%% Note: The variales - 'BrokenTensionBond' stores the information of the 
%% broken bonds during the simulations. Since, each bond has two identities
%% hence, total number of bonds broken are half of the length of this variable
%% The even or odd rows looping will remove redundancy

ReleaseDisp = unique(BrokenTensionBond(:,3));
NBrokenBonds = zeros(length(ReleaseDisp),1);
EnergyReleased = zeros(length(ReleaseDisp),1);
for i = 1:length(ReleaseDisp)
    for k = 1:2:TensionBondCount-1
        if(BrokenTensionBond(k,3) == ReleaseDisp(i))
            NBrokenBonds(i) = NBrokenBonds(i) + 1;
            EnergyReleased(i) = EnergyReleased(i) + 0.5*Kn*(BrokenTensionBond(k,4)^2 ...
                + TNRatio*(BrokenTensionBond(k,5) + a*BrokenTensionBond(k,7))^2 ...
                + TNRatio*a*BrokenTensionBond(k,6)*(BrokenTensionBond(k,5) + a*BrokenTensionBond(k,7)) ...
                + RTRatio*TNRatio*BrokenTensionBond(k,6)^2);
        end
    end
end

figure(3)
h = plot(ReleaseDisp, EnergyReleased,'s');
set(h, 'LineWidth', 1.5)
xlabel('\epsilon \rightarrow','Fontsize',16,'FontName',...
    'Arial','FontWeight', 'bold'), 
ylabel('\sigma \rightarrow','Fontsize',16,'FontName', ...
    'Arial','FontWeight', 'bold'),
set(gca,'FontWeight', 'bold');
dim = [.15 .6 .3 .3];
str = ['200X200 Lattice'];
t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
t.FontSize = 12;
dim = [.15 .5 .3 .3];
str = ['p = 30%'];
t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
t.FontSize = 12;
% legend({'Top','Bottom'}, ...
%     'Location','south','Fontsize',10),
    
        