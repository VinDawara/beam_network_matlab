%% Field interpolation within the domain using delaunay triangulation


x = zeros(TotalNodes,1);
y = zeros(TotalNodes,1);
Field = zeros(TotalNodes,1);

for ID = 1:TotalNodes
    x(ID) = node(ID).r_initial(1);
    y(ID) = node(ID).r_initial(2);
    Eigen = eig(node(ID).Stress);
    Field(ID) = max(Eigen);
%     Field(ID) = node(ID).Stress(1,2);
end
syT = abs(Sigma(step,2));
Field = Field/syT; % divide this by remote loading
DT = delaunayTriangulation(x,y);
Nx = 5*N;
Ny = 5*M;
xq = linspace(min(x(2)),max(x(TotalNodes-1)),Nx);
yq = linspace(min(y),max(y),Ny);
Fieldq = zeros(Ny,Nx);

% Crack domain
% CrackPoly = zeros(4,2);
% CrackPoly(1,:) = node(CrackNodes(1)-1).r;
% CrackPoly(2,:) = node(CrackNodes(end)+1).r;
% CrackPoly(3,:) = node(node(CrackNodes(end)+1).neighbors((node(CrackNodes(end)+1).aph == 120))).r;
% CrackPoly(4,:) = node(node(CrackNodes(1)-1).neighbors((node(CrackNodes(1)-1).aph == 120))).r;
% CrackPoly = polyshape(CrackPoly);

for i = 1:Ny
    for j = 1:Nx
%         if(isinterior(CrackPoly, xq(j), yq(i)))
%             Fieldq(i,j) = 0;
%         else    
            TriangleID = pointLocation(DT,xq(j),yq(i));
            if(~isnan(TriangleID))
                VertexID = DT.ConnectivityList(TriangleID,:);
                Vertex = DT.Points(VertexID,:);
                d  = zeros(3,1);
                for k = 1:3
                    d(k) = norm([(Vertex(k,1) - xq(j));(Vertex(k,2) - yq(i))]);
                end
                w = (1./d)/sum(1./d);
                Fieldq(i,j) = w(1)*Field(VertexID(1)) + w(2)*Field(VertexID(2))...
                    + w(3)*Field(VertexID(3));
            end
%         end  
    end
end
        
% figure(1)
[X,Y] = meshgrid(xq,yq);

surfc(X,Y,Fieldq)
view([0 90])
% caxis([0 1])
colorbar
shading interp
axis tight
title('\sigma_{yy}/\sigma_{o}')
xlabel('x \rightarrow','Fontsize',16,'FontName',...
    'Arial','FontWeight', 'bold'), 
ylabel('y \rightarrow','Fontsize',16,'FontName', ...
    'Arial','FontWeight', 'bold'),
zlabel('\sigma_{xx} \rightarrow','Fontsize',16,'FontName', ...
    'Arial','FontWeight', 'bold'),
set(gca,'FontWeight', 'bold');
% xlim([node(2,1).r(1),node(1,N).r(1)])
% ylim([node(1,1).r(2),node(M,1).r(2)])

% Plotting the quarter image
% First quarter
% X = X(0.5*Ny:Ny,0.5*Nx:Nx) - X(0.5*Ny,0.5*Nx);
% Y = Y(0.5*Ny:Ny,0.5*Nx:Nx) - Y(0.5*Ny,0.5*Nx);
% Fieldq = Fieldq(0.5*Ny:Ny,0.5*Nx:Nx);




