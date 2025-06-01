%% Tracing crack path with the pore spatial distribution
%% Date: 04/04/2022

L = (N-1)*a;
W = (M-1)*0.5*sqrt(3)*a;

figure('visible','off')
% [E,R] = voronoi(NodePos(:,1),NodePos(:,2));
% plot(E,R, 'k','Marker','none')
% axis([-0.5 L 0 W])
% hold on


k = 1;
l = 1;
Pos = zeros(TotalNodes,2);
Network = [];
for ID = 1:TotalNodes
    if(~any(ID == UpperBoundary))
        for s = 1:length(node(ID).neighbors)
            if(any(node(ID).aph(s) == [0, 60, 120]))
                if(node(ID).Kn(s) ~= 0)
                    Network(k,1) = ID;
                    Network(k,2) = node(ID).neighbors(s);
                    k = k +1;
                end
            end
        end
    end
    l = l + 1;
end
Network(k,:) = [UpperBoundary(end-1),UpperBoundary(end)];
ID = TotalNodes - N;
rowx = node(ID).r_initial(1)-1:0.25:node(ID + N).r(1)+ 1;
rowy = node(ID).r_initial(2)*ones(size(rowx));
plot(rowx,rowy,'k','LineWidth',3')
hold on
ID = 1;
rowx = node(ID).r_initial(1)- 1:0.25:node(ID+N).r(1)+ 1;
rowy = node(ID).r_initial(2)*ones(size(rowx));
plot(rowx,rowy,'k','LineWidth',3')
hold on
xlim([(node(1).r_initial(1)-a) (node(TotalNodes).r(1)+a)]);
ylim([(node(1).r_initial(2)-a) (node(TotalNodes).r(2)+a)]); 
axis equal
g = graph(Network(:,1),Network(:,2));
plot(g,'k','Marker','none','EdgeAlpha',0.4, 'XData', NodePos(1:TotalNodes,1), 'YData', NodePos(1:TotalNodes,2));
hold on

% Plotting the Broken Bonds (Crack path) in the output
Vx = zeros(0.5*(TensionBondCount-1),2);
Vy = Vx;
if(TensionBondCount ~=1)
    for i = 1:2: length(BrokenTensionBond(:,1))
        Vx(0.5*(i+1),:) = [node(BrokenTensionBond(i,1)).XEdge(BrokenTensionBond(i,2),1), node(BrokenTensionBond(i,1)).XEdge(BrokenTensionBond(i,2),2)];
        Vy(0.5*(i+1),:) = [node(BrokenTensionBond(i,1)).YEdge(BrokenTensionBond(i,2),1), node(BrokenTensionBond(i,1)).YEdge(BrokenTensionBond(i,2),2)];
    end
    plot([Vx(:,1)';Vx(:,2)'],[Vy(:,1)'; Vy(:,2)'],'r','LineWidth',2)
    hold on
end


% Plotting pores as individual polygons
% for i = 1:length(Poresize(:,1))
%     for j = 1:nnz(Poresize(i,:))
%         if (j == 1)
%             poly = polyshape(node(Poresize(i,j)).Voro_Polygon);
%         else
%             poly = union(poly,polyshape(node(Poresize(i,j)).Voro_Polygon));
%         end
%     end
%     plot(poly)
%     hold on
% end

% Plotting pores as one single polygon
% for i = 1:length(PoreID)
%     ID = PoreID(i);
%     if(i == 1)
%         PoreStructure = polyshape(node(ID).Voro_Polygon);
%     else
%         PoreStructure = union(PoreStructure,polyshape(node(ID).Voro_Polygon));
%     end
%     FinalVolume =  area(PoreStructure);
% end

if(~porosity == 0)
    plot(PoreStructure,'Facecolor',[0.3010 0.7450 0.9330], 'LineWidth',0.3, 'EdgeAlpha',0.5)
end
% title(['STEPS = ',num2str(step)]);
axis([-0.5 L 0 W])
axis off
hold off