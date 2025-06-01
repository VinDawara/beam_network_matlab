%% BEAM MODEL (TRIANGULAR LATTICE)
%% Dated: 11/08/2021
%% Problem: Central crack/Edge Crack in an plate with free boundary conditions
%% Test type: Tensile/Compression
tic
clear all
close all
disp('SIMULATION BEGINS')
disp('----------------------------------------------------------------')

% Lattice points
M = 100; % M should be even (Vertical points)
N = 100; % Horizontal points
porosity = 0.2;      % porosity level 
ao = 20;

% Creating folder to save data
MyFolder = [num2str(M),'X',num2str(N),'p',num2str(porosity*100),'ao',num2str(ao)];
mkdir(sprintf('%s',MyFolder))
mkdir(sprintf('%s/Data',MyFolder))
mkdir(sprintf('%s/Snapshots',MyFolder))


% Bond properties
lambda = 125;       % scaling factor
Kn = 1;             % bond normal stiffness
TNRatio = 0.1;      % Tangential to Normal stiffness ratio
RTRatio = 1/3;      % Rotational to tangential stiffness

% Creating triangular mesh
TriangularMesh
disp('STAGE 1: Mesh created')


% Introducing percolation porosity  
PorousConfg
disp(['STAGE 2: Porosity ',num2str(porosity),' introduced'])

% Defining voronoi cell at each node
VoronoiEdges
disp('STAGE 3: Voronoi cells are created at each node')

% Define simulation variables
index = 1;              % Remote load control variable
TensionBondCount = 1;   % Counts number of bond broken due to tension
BrokenTensionBond = []; % Stores the broken bond IDs 

% % Estimating initial pore distribution
if(~porosity == 0)
    PoreDistribution
    disp('STAGE 4: Pore size distribution is estimated')
else
    disp('STAGE 4: No porosity in the lattice')
end

% Intrdrocing central crack of length a_o
CentralCrack

% Initial Lattice
LatticeView
disp('STAGE 5: Displaying and saving initial lattice')
print(gcf,sprintf('%s/InitialLattice.png',MyFolder),'-dpng','-r1000')
disp('         Lattice saved.')

% Loop variables
delta = zeros(5,1);
Force = zeros(5,4);
Sigma = zeros(5,4);
Strain = zeros(5,1);
Energy = zeros(5,1);
Energy_in = zeros(5,1);
Energy_Axial = zeros(5,1);
Energy_Bending = zeros(5,1);
if(~porosity == 0)
    CN = zeros(5,7);
else
    CN = zeros(5,6);
end


loop = 300;
disp('STAGE 6: Applying loading and solving the lattice')
for step = 1:loop
    % BC on bottom and top nodes
    yB = 0;
    yT = -2*0.05*a*index;  
    disp('---------------------------------------------------------------')
    disp(['         Step =  ',num2str(step),' Remote disp. = ', num2str(yT)])
    if(step == 1)
        MatrixEquationSparse
    else
        MatrixAUpdate
        UpdateB
    end
    Solver
    disp(['         Solution is obtained at ', num2str(itr),' iterations'])
    disp(['         Relative residual = ',num2str(relresi)])

    ForceComponents
    ForceEstimation
    disp('         Forces and moment are computed')
   
    % Storing force-displacemnt or stress-strain values
    delta(step+1) = abs(yT);
    Strain(step+1) = delta(step+1)/((M-1)*(sqrt(3)/2)*a);
    Force(step+1,:) = [FT, FB];
    Sigma(step+1,:) = [SigT, SigB];
    
    % Applying fracture criterion
    FractureCriterion
      
    SystemEnergy

    Energy_in(step+1) = trapz(delta(1:step+1),abs(Force(1:step+1,2)));
    Energy(step+1) = SysEnergy;
    Energy_Axial(step+1) = SysEnergy_Axial;
    Energy_Bending(step+1) = SysEnergy_Bending;
    

    if(BondBreak == 0)  % if no critial bond found
        index = index + 1;
        disp('         No Critical bond found. Increasing the load')
    else               % if critial bond found 
        disp('         Critical bond found. Most critial ones are broken.')
        disp('         Re-equilibrating at same load.')
    end
    
    % Lattice visualisation with crack path
    if(mod(step,10) == 0)
    disp('----------SAVING LATTICE IMAGE-----------------------')    
    clf
    CrackPath
    print(gcf,sprintf('%s/Snapshots/%d.png',MyFolder, step),'-dpng','-r800');
    disp('----------LATTICE IMAGE SAVED-----------------------')
    end
    
    % Nodal coordination number distribution
    NodeCN
    CN(step,:) = NodeDis; 
    
    % When the loop will break
    if(abs(FT(2)) <1e-6 || abs(FB(2)) <1e-6)
        disp('         Simulation stops because remote forces are below tolerance limit.')
        break;
    elseif(step == loop)
        disp('         Simulation stops because maximum iteration limit is reached.')
    end 
    
     

    for ID = 1:TotalNodes
        node(ID).r = node(ID).r_initial;
    end
    if(mod(step,50) == 0)
        save(sprintf('%s/Data/%d.mat',MyFolder, step))
        disp('         SAVING WORKSPACE')
    end
end

% Data reterival
disp('Stage 7: Retriving data')


figure(2)
h = plot(delta, Energy_in, delta, Energy, delta, Energy_Axial, delta, Energy_Bending);
set(h, 'LineWidth', 1.5)
xlabel('\delta \rightarrow','Fontsize',16,'FontName',...
    'Arial','FontWeight', 'bold'), 
ylabel('\Phi \rightarrow','Fontsize',16,'FontName', ...
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
legend({'Input', 'Total','Axial','Bending'}, ...
    'Location','north','Fontsize',10),

figure(3)
h = plot(Strain, abs(Sigma(:,2)),'-',Strain, abs(Sigma(:,4)),'-');
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
legend({'Top','Bottom'}, ...
    'Location','south','Fontsize',10),


figure(4)
timestep = 1:step;
for i = 1:6
    h = plot(timestep, CN(:,i));
    set(h, 'LineWidth', 1.5)
    hold on
end
xlim([0,step+10])
xlabel('iterations \rightarrow','Fontsize',16,'FontName',...
    'Arial','FontWeight', 'bold'), 
ylabel('CN \rightarrow','Fontsize',16,'FontName', ...
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
legend({'1','2','3','4','5','6'}, ...
    'Location','northeast','Fontsize',10),
hold off

Load_Step = 2*0.05*a*(1:index);
BrokenBondAnalyzer
figure(5)
h = plot(Load_Step, NumBondBroken,'-o');
set(h, 'LineWidth', 1.5)
xlabel('\delta \rightarrow','Fontsize',16,'FontName',...
    'Arial','FontWeight', 'bold'), 
ylabel('Fy \rightarrow','Fontsize',16,'FontName', ...
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
%------------------------------------------------------------------------
disp('         Data plotted. Simulation ends successfully!')


% Removing unneccessary variables
clrvars = {'i', 'I', 'J', 'Entry','k','ID','s','l','v','CriticalNode','SigmaT','maxStress','BondDeflect',...
     'Pos', 'Network','rowx','rowy','index','t','dim','str','j'};
clear(clrvars{:})
clrvars = {'SysEnergy', 'BondBreak','FT','FB','flag','relresi','itr','maxit','tol','A','B','SigB','SigT','NoBondNode'};
clear(clrvars{:})

% Saving output data
save(sprintf('%s/Data/Output.mat',MyFolder))


