%% COMPUTATION OF NORMAL - SHEAR FORCE AND MOMENT COMPONENTS
%% DATE: 13/08/2021

%% Uncomment the codes, if breaking criterion is maximum normal or shear stress of the bond
% % Defining critical strength of the bonds
% CriticalShear = 0.03;
% CriticalTension = 0.01;
% 
% BondBreak = 0;     % Indicator which tells whether critical bonds are 
%                     % present(1) or not(0)
% 
% % Variable to store critical bond information
% CriticalNode = zeros(1,3);
% k = 1;

for ID = 1:TotalNodes
    node(ID).Fn = zeros(size(node(ID).neighbors));
    node(ID).Ft = zeros(size(node(ID).neighbors));
    node(ID).M = zeros(size(node(ID).neighbors));
    for s = 1:length(node(ID).neighbors)
        % Calculating normal and shear forces
        node(ID).Fn(s) = node(ID).Kn(s)*((node(ID).u(1) - node(node(ID).neighbors(s)).u(1))*cosd(node(ID).aph(s)) + ...
            (node(ID).u(2) - node(node(ID).neighbors(s)).u(2))*sind(node(ID).aph(s)));
        node(ID).Ft(s) = node(ID).Kt(s)*(-(node(ID).u(1) - node(node(ID).neighbors(s)).u(1))*sind(node(ID).aph(s)) + ...
            (node(ID).u(2) - node(node(ID).neighbors(s)).u(2))*cosd(node(ID).aph(s))) + ...
            0.5*node(ID).Kt(s)*a*(node(ID).u(3) + node(node(ID).neighbors(s)).u(3));
        node(ID).M(s) = 0.5*a*node(ID).Kt(s)*(-(node(ID).u(1) - node(node(ID).neighbors(s)).u(1))*sind(node(ID).aph(s)) + ...
            (node(ID).u(2) - node(node(ID).neighbors(s)).u(2))*cosd(node(ID).aph(s))) + ...
            node(ID).Kphi(s) *a^2*(node(ID).u(3) + 0.5*node(node(ID).neighbors(s)).u(3));
        
        % Determining the stress by dividing the force with the respective voronoi face area
        FaceArea = 1;
        node(ID).Fn(s) = lambda*node(ID).Fn(s)/FaceArea;
        node(ID).Ft(s) = lambda*node(ID).Ft(s)/FaceArea;
        node(ID).M(s) = lambda*node(ID).M(s);
        
%         % Checking whether the bond is critical and then storing in
%         % CriticalNode
%         if(~any(ID == FixedNodes))
%             if(abs(node(ID).Ft(s))>=CriticalShear || node(ID).Fn(s) >= CriticalTension)
%                 Eff_stress = sqrt(node(ID).Fn(s)^2 + node(ID).Ft(s)^2);
%                 CriticalNode(k,:) = [ID s Eff_stress];
%                 k = k+1;
%                 BondBreak = 1;
%             end
%         end
    end 
end

% Removing the bond which have maximum effective stresses
% if(BondBreak == 1)
%     maxStress = max(CriticalNode(:,3));
%     for i = 1:length(CriticalNode(:,3))
%         if(abs(maxStress - CriticalNode(i,3))<=1e-6)
%             node(CriticalNode(i,1)).Kn(CriticalNode(i,2)) = 0;
%             node(CriticalNode(i,1)).Kt(CriticalNode(i,2)) = 0;
%             node(CriticalNode(i,1)).Kphi(CriticalNode(i,2)) = 0;
%             BrokenTensionBond(TensionBondCount,:) = CriticalNode(i,1:2);
%             TensionBondCount = TensionBondCount + 1; 
%         end
%     end
% end  



