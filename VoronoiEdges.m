%% DEFINING NEIGHBORING EDGES
%% Date: 04/04/2022

%% STEP 1: Generating voronoi set

NodePos = zeros(TotalNodes,2);
for ID = 1:TotalNodes
    NodePos(ID,:) = node(ID).r;
end

% Artificial points to define the rectangular boundary
k = 1;
Arti_ID = 0;
dx = a;
dy = 0.5*sqrt(3)*a;
ID = ID + 1;

% Bottom of the lattice
i = 0;
for j = 1:N+1
    NodePos(ID,1) = (j-1)*dx - 0.5*dx;
    NodePos(ID,2) = (i-1)*dy;
    Arti_ID(k) = ID;
    k = k + 1;
    ID = ID + 1;
end

% Left of the lattice
j = -0.5;
for i = 1:2:M+1
    NodePos(ID,1) = (j-1)*dx;
    NodePos(ID,2) = (i-1)*dy;
    Arti_ID(k) = ID;
    k = k + 1;
    ID = ID + 1;
end

% Right of the lattice
j = N+1;
for i = 2:2:M-1
    NodePos(ID,1) = (j-1)*dx - 0.5*dx;
    NodePos(ID,2) = (i-1)*dy;
    Arti_ID(k) = ID;
    k = k + 1;
    ID = ID + 1;
end

% Top of the lattice
i = M +1;
for j = 1:N
    Arti_ID(k) = ID;
    k = k + 1;
    if(mod(M,2) == 0)
       NodePos(ID,1) = (j-1)*dx ;
       NodePos(ID,2) = (i-1)*dy;
       ID = ID + 1; 
    else
        NodePos(ID,1) = (j-1)*dx - 0.5*dx;
        NodePos(ID,2) = (i-1)*dy;
        ID = ID + 1;
    end
end

% One in the row of last node (Auxillary node)
i = M;
j = N+1;
Arti_ID(k) = ID;
NodePos(ID,1) = (j-1)*dx ;
NodePos(ID,2) = (i-1)*dy;

DT = delaunayTriangulation(NodePos(:,1),NodePos(:,2));

[Voro_Edges,Voro_Region] = voronoiDiagram(DT);
% figure(1)
% voronoi(DT)




%% STEP 2: DEFINING THE VORONOI POLYGON AT EACH NUCLEUS

L = (N-1)*a;
W = (M-1)*0.5*sqrt(3)*a;
for ID = 1:TotalNodes
    PolyVertices = Voro_Edges(Voro_Region{ID},:);
    poly = polyshape(PolyVertices);
    Linesegment = [-10 0;(L+10) 0];
    in = intersect(poly,Linesegment);
    in(:,2) = 0*(abs(in(:,2))<= 1e-9);
    if(~isempty(in))
        PolyVertices = PolyVertices((PolyVertices(:,2)>0),:);
        PolyVertices = [PolyVertices; in];
        PolyVertices = AntiClockwise(PolyVertices);
    end
    poly = polyshape(PolyVertices);
    Linesegment = [-0.5 -10;-0.5 (W+10)];
    in = intersect(poly,Linesegment);
    in(:,1) = -0.5*(abs(in(:,1)+0.5)<=1e-9);
    if(~isempty(in))
        PolyVertices = PolyVertices((PolyVertices(:,1)>-0.5),:);
        PolyVertices = [PolyVertices; in];
        PolyVertices = AntiClockwise(PolyVertices);
    end
    poly = polyshape(PolyVertices);
    Linesegment = [L -10;L (W+10)];
    in = intersect(poly,Linesegment);
    in(:,1) = L*(abs(in(:,1)-L)<=1e-9);
    if(~isempty(in))
        PolyVertices = PolyVertices((PolyVertices(:,1)<L),:);
        PolyVertices = [PolyVertices; in];
        PolyVertices = AntiClockwise(PolyVertices);
    end
    poly = polyshape(PolyVertices);
    Linesegment = [-10 W;(L+10) W];
    in = intersect(poly,Linesegment);
    in(:,2) = W*(abs(in(:,2)-W)<=1e-9);
    if(~isempty(in))
        PolyVertices = PolyVertices((PolyVertices(:,2)<W),:);
        PolyVertices = [PolyVertices; in];
        PolyVertices = AntiClockwise(PolyVertices);
    end
    
    
    node(ID).Voro_Polygon = round(PolyVertices,6);
%     figure(1)
%     plot(polyshape(node(ID).Voro_Polygon))
%     hold on
%     plot(node(ID).Voro_Polygon(:,1),node(ID).Voro_Polygon(:,2),'k')
%     axis([-1 L+1 -1 W+1])
%     hold on

end

%% STEP 3: Computing the voronoi edges from the 'Voro_Polygon'
for ID = 2:TotalNodes-1
    s = 1;
    for k = 1:length(node(ID).neighbors)
        [X,Y] = polyxpoly(node(ID).Voro_Polygon(:,1),node(ID).Voro_Polygon(:,2), ...
        node(node(ID).neighbors(s)).Voro_Polygon(:,1),node(node(ID).neighbors(s)).Voro_Polygon(:,2));
        
        if(isempty(X))                  % No contact
            node(ID).neighbors(s) = [];
        elseif(length(X) == 1)          % single point contact
            node(ID).XEdge(s,1:2) = X;
            node(ID).YEdge(s,1:2) = Y;
        else                            % line contact
            node(ID).XEdge(s,:) = X;
            node(ID).YEdge(s,:) = Y;
        % Calculating intercept
            P1 = [(node(ID).XEdge(s,1) - node(ID).r(1)), (node(ID).YEdge(s,1) - node(ID).r(2))];
            P2 = [(node(ID).XEdge(s,2) - node(ID).r(1)), (node(ID).YEdge(s,2) - node(ID).r(2))];
            intercept = P1(2) - P1(1)*((P2(2) - P1(2))/(P2(1) - P1(1)));
        
        % Storing edges in the anticlockwise sense by checking the
        % intercept
            temp1 = [node(ID).XEdge(s,1), node(ID).YEdge(s,1)];
            temp2 = [node(ID).XEdge(s,2), node(ID).YEdge(s,2)];
            if(isinf(intercept))
                if(P1(1) > 0) 
                    node(ID).YEdge(s,1) = min(temp1(2),temp2(2));
                    node(ID).YEdge(s,2) = max(temp1(2), temp2(2));
                else
                    node(ID).YEdge(s,1) = max(temp1(2),temp2(2));
                    node(ID).YEdge(s,2) = min(temp1(2), temp2(2));
                end
            elseif(intercept > 0)
                if(temp1(1) <= temp2(1))
                    node(ID).XEdge(s,1) = temp2(1);
                    node(ID).YEdge(s,1) = temp2(2);
                    node(ID).XEdge(s,2) = temp1(1);
                    node(ID).YEdge(s,2) = temp1(2);
                end
            else
                if(temp1(1) >= temp2(1))
                    node(ID).XEdge(s,1) = temp2(1);
                    node(ID).YEdge(s,1) = temp2(2);
                    node(ID).XEdge(s,2) = temp1(1);
                    node(ID).YEdge(s,2) = temp1(2);
                end
            end
            
%         % Calculating angles made by the outward normal
%             t = [(node(ID).XEdge(s,2)-node(ID).XEdge(s,1)), (node(ID).YEdge(s,2)-node(ID).YEdge(s,1))];
%             if(t(1) >= 0 && t(2) >= 0)
%                 node(ID).aph(s) = atand(t(2)/t(1)) -90;     % 90 + theta = alpha
%             elseif(t(1) < 0 && t(2) > 0)
%                 node(ID).aph(s) = 180 - atand(t(2)/abs(t(1))) - 90;
%             elseif(t(1) < 0 && t(2) <= 0)
%                 node(ID).aph(s) = 180 + atand(abs(t(2))/abs(t(1))) - 90;
%             else
%                 node(ID).aph(s) = -atand(abs(t(2))/(t(1))) - 90;
%             end
             
            s = s + 1;
        end
    end        
end


clrvars = {'dx','dy', 'Arti_ID', 'DT', 'Voro_Edges', 'Voro_Region','PolyVertices', 'X', 'Y', 'P1', 'P2', ...
    'intercept', 'temp1', 'temp2', 'poly','in','Linesegment'};
clear(clrvars{:})
