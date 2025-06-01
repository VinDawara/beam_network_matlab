%% Lattice visualisation
%% Date: 22/09/2021

k = 1;
l = 1;
Pos = zeros(TotalNodes,3);
Pos(:,3) = 40*ones(TotalNodes,1);
Network = [];
for ID = 1:TotalNodes
    Pos(l,1:2) = node(ID).r;
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

ID = TotalNodes - N;
rowx = node(ID).r(1)-1:0.25:node(ID + N).r(1)+ 1;
rowy = node(ID).r(2)*ones(size(rowx));
rowz = 40*ones(size(rowx));
% figure('WindowState','maximized')
% figure('visible','on')
plot3(rowx,rowy, rowz,'k','LineWidth',3')
hold on
ID = 1;
rowx = node(ID).r(1)- 1:0.25:node(ID+N).r(1)+ 1;
rowy = node(ID).r(2)*ones(size(rowx));
rowz = 40*ones(size(rowx));
plot3(rowx,rowy,rowz,'k','LineWidth',3')
hold on
xlim([(node(1).r(1)-a) (node(TotalNodes).r(1)+a)]);
ylim([(node(1).r(2)-a) (node(TotalNodes).r(2)+a)]); 
axis equal
g = graph(Network(:,1),Network(:,2));
plot(g,'k','Marker','none','EdgeAlpha',1, 'XData', Pos(:,1), 'YData', Pos(:,2),'ZData', Pos(:,3));
view([0 90])
axis equal
hold off
