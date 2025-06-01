%% Computing pore size distribution
%% Date: 06/04/22

% STEP 1: Construct the graph of the pore nodes
k = 1;
Link = [];
temp = PoreID;
for i = 1:(length(PoreID)-sum(PoreID==0))
    if(PoreID(i)~=0)
        ID = PoreID(i);
        Poreneighbor = intersect(node(ID).neighbors,PoreID);
        PoreID(PoreID == ID) = 0;
        if(~isempty(Poreneighbor))
            for s = 1:length(Poreneighbor)
                Link(k,1) = ID;
                Link(k,2) = Poreneighbor(s);
                k = k + 1;
            end
        end
    end
end
PoreID = temp;

% for i = 1:length(PoreID)
%     ID = PoreID(i);
%     Poreneighbor = intersect(node(ID).neighbors,PoreID);
%     if(~isempty(Poreneighbor))
%         for s = 1:length(Poreneighbor)
%             Link(k,1) = ID;
%             Link(k,2) = Poreneighbor(s);
%             k = k + 1;
%         end
%     end    
% end



G = graph;
G = addnode(G,TotalNodes);
G = addedge(G,Link(:,1),Link(:,2));
% plot(G,'k','Marker','.','EdgeAlpha',1, 'XData', NodePos(1:TotalNodes,1), 'YData', NodePos(1:TotalNodes,2));
% hold on

% STEP 2: Estimate the pore size
Poresize = [];
temp = PoreID;
k = 1;
for i = 1:length(PoreID)
    if(PoreID(i)~= 0)
        ConnNode = nearest(G,PoreID(i),Inf);
        if(isempty(ConnNode))
            Poresize(k,1) = PoreID(i);
            PoreID(PoreID == PoreID(i)) = 0;
            k = k +1;
        else
            for j = 1:length(ConnNode)
                PoreID(PoreID == ConnNode(j)) = 0;
            end
            Poresize(k,1:(length(ConnNode)+1)) = [PoreID(i), ConnNode'];
            PoreID(PoreID == PoreID(i)) = 0;
            k = k + 1;
        end
    end
end
PoreID = temp;

% % STEP 3: Plot the pore size as a polygon
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

% STEP 4: Estimating pore size distribution
% Method 1: Counting number of nodes
PoreLen = zeros(length(Poresize(:,1)),1);
for i = 1:length(Poresize(:,1))
    PoreLen(i) = nnz(Poresize(i,:));
end
% figure(2)
% histogram(PoreLen)
% PoreDist = histcounts(PoreLen);

% Method 2: Calculating the perimeter of the pore polygon
PoreParameter = zeros(length(Poresize(:,1)),1);
for i = 1:length(Poresize(:,1))
    for j = 1:nnz(Poresize(i,:))
        if(j == 1)
            Poly = polyshape(node(Poresize(i,j)).Voro_Polygon);
        else
            Poly = union(Poly,polyshape(node(Poresize(i,j)).Voro_Polygon));
        end
    end
    PoreParameter(i) = perimeter(Poly)/(6*a/sqrt(3));
end
% figure(2)
% histogram(PoreParameter)
% PoreDist = histcounts(PoreLen);

% STEP 5: Forming all pore polygon as one polygon
for i = 1:length(PoreID)
    ID = PoreID(i);
    if(i == 1)
        PoreStructure = polyshape(node(ID).Voro_Polygon);
    else
        PoreStructure = union(PoreStructure,polyshape(node(ID).Voro_Polygon));
    end
    FinalVolume =  area(PoreStructure);
end
% plot(PoreStructure,'Facecolor',[0.3010 0.7450 0.9330],'Edgecolor','none')


clrvars = {'Link', 'temp', 'Poreneighbor','ConnNode','Poly'};
clear(clrvars{:})


            
            