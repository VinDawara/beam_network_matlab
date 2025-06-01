%% Computing Node degree
%% Date: 08/04/2022

% STEP 1: Create the graph of bonded network
k = 1;
l = 1;
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
g = graph(Network(:,1),Network(:,2));
NodeDeg = degree(g);
% figure(3)
% histogram(NodeDeg)
NodeDis = histcounts(NodeDeg);