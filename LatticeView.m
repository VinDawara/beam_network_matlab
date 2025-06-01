%% Lattice visualisation
%% Date: 22/09/2021

L = (N-1)*a;
W = (M-1)*0.5*sqrt(3)*a;
k = 1;
l = 1;
Pos = zeros(TotalNodes,2);
Network = [];
for ID = 1:TotalNodes
    Pos(l,:) = node(ID).r;
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
rowx = node(ID).r(1)-1:0.25:node(ID + N).r(1)+ 1;
rowy = node(ID).r(2)*ones(size(rowx));
% figure(2)
% figure('WindowState','maximized')
figure('visible','off')
plot(rowx,rowy,'k','LineWidth',3')
hold on
ID = 1;
rowx = node(ID).r(1)- 1:0.25:node(ID+N).r(1)+ 1;
rowy = node(ID).r(2)*ones(size(rowx));
plot(rowx,rowy,'k','LineWidth',3')
hold on
xlim([(node(1).r(1)-a) (node(TotalNodes).r(1)+a)]);
ylim([(node(1).r(2)-a) (node(TotalNodes).r(2)+a)]); 
axis equal
g = graph(Network(:,1),Network(:,2));
plot(g,'k','Marker','.','EdgeAlpha',0.5, 'XData', Pos(:,1), 'YData', Pos(:,2));
hold on
if(~porosity == 0)
    plot(PoreStructure,'Facecolor',[0.3010 0.7450 0.9330])
end
axis([-0.5 L 0 W])
axis off
hold off
