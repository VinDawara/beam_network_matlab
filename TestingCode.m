clear all
M = 10; % M should be even
N = 10;
porosity = 0.1;
lambda = 48*100;
Kn = 1;
TNRatio = 0.1;    % Tangential to Normal stiffness ratio
RTRatio = 1/3;    % Rotational to tangential stiffness
TriangularMesh
PorousConfg
VoronoiEdges

TensionBondCount = 1;

index = 10;
step = 1;
yB = 0.005*a*index;
yT = -0.005*a*index;
disp('---------------------------------------------------------------')
disp(['Computing at steps =  ',num2str(step),' For remote disp. = ', num2str(2*yT)])
MatrixEquationSparse
Solver
disp(['Solution is obtained at ', num2str(itr),' iterations'])
disp(['Relative residual = ',num2str(relresi)])
ForceComponents
FractureCriterion
% ValidEquilibrium
ForceEstimation
NodalStress
FieldVisualisation
% hold on
% LatticeView3D



% x = zeros(TotalNodes,1);
% y = x;
% for ID = 1:TotalNodes
%     x(ID) = node(ID).r(1);
%     y(ID) = node(ID).r(2);
% end
% 
hold on
vxlabels = arrayfun(@(n) {sprintf('P%d', n)}, (1:TotalNodes)'); 
Hpl = text(NodePos(1:TotalNodes,1),NodePos(1:TotalNodes,2),vxlabels,'FontWeight','bold','HorizontalAlignment',...
   'center','BackgroundColor','none');
hold off