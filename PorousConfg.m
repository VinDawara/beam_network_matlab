%% Introducing porosity by removing nodes
%% Date: 20/10/2021


TotalPorousNode = ceil(porosity*(M-1)*(N-1));   


RandNum = UpperBoundary(1) - LowerBoundary(end)+1;
r = randi([LowerBoundary(end)+1,UpperBoundary(1)-1],RandNum,1);
r = unique(r,'stable');
for i = 1:TotalPorousNode
    for s = 1:length(node(r(i)).neighbors)
        node(r(i)).Kn(s) = 0;
        node(r(i)).Kt(s) = 0;
        node(r(i)).Kphi(s) = 0;
        node(node(r(i)).neighbors(s)).Kn(node(node(r(i)).neighbors(s)).neighbors == r(i)) = node(r(i)).Kn(s);
        node(node(r(i)).neighbors(s)).Kt(node(node(r(i)).neighbors(s)).neighbors == r(i)) = node(r(i)).Kt(s);
        node(node(r(i)).neighbors(s)).Kphi(node(node(r(i)).neighbors(s)).neighbors == r(i)) = node(r(i)).Kphi(s);
    end   
end
PoreID = r(1:TotalPorousNode);


clrvars = {'RandNum','r'};
clear(clrvars{:})
