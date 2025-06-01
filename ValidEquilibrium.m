%% VALIDATION OF THE EQUILIBRIUM CONDITIONS
%% DATE: 13/08/2021

F = zeros(TotalNodes,3);
for ID = 1:TotalNodes
    for s = 1:length(node(ID).neighbors)
        
        % X component of the force
        F(ID,1) = F(ID,1) + node(ID).Fn(s)*cosd(node(ID).aph(s));
        F(ID,1) = F(ID,1) - node(ID).Ft(s)*sind(node(ID).aph(s));
        
        % Y component of the force
        F(ID,2) = F(ID,2) + node(ID).Fn(s)*sind(node(ID).aph(s));
        F(ID,2) = F(ID,2) + node(ID).Ft(s)*cosd(node(ID).aph(s));
        
        % moment at the node
        F(ID,3) = F(ID,3) + node(ID).M(s);
        
    end
end