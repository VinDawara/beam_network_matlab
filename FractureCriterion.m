%% BREAKING CRITERION
%% BASED ON AXIAL AND BENDING STRESS
%% DATE: 15/08/2021

alpha = 0.5;
h = sqrt(TNRatio)*a;
Area = h^2;
W = h^3/6;

% Defining critical strength of the bonds
CriticalTension =   0.52*lambda;

BondBreak = 0;     % Indicator which tells whether critical bonds are 
                    % present(1) or not(0)

% Variable to store critical bond information
CriticalNode = zeros(1,3);
k = 1;

for ID = 1:TotalNodes
    for s = 1:length(node(ID).neighbors)
        SigmaT = -node(ID).Fn(s)/Area + (alpha/W)*max(abs(node(ID).M(s)), ...
            abs(node(node(ID).neighbors(s)).M(node(node(ID).neighbors(s)).neighbors == ID)));

        % Checking whether the bond is critical and then storing in
        % CriticalNode
%         if(~any(ID == BoundaryNodes) || ~any(node(ID).neighbors(s) == BoundaryNodes))
            if(~any(ID == FixedNodes))
                if(SigmaT >= CriticalTension)
                    CriticalNode(k,:) = [ID s SigmaT];
                    k = k+1;
                    BondBreak = 1;
                end
            end
%         end
    end
end

% Removing the bond which have maximum effective stresses
if(BondBreak == 1)
    maxStress = max(CriticalNode(:,3));
    for i = 1:length(CriticalNode(:,3))
        if(abs(maxStress - CriticalNode(i,3))<=1e-3)
            node(CriticalNode(i,1)).Kn(CriticalNode(i,2)) = 0;
            node(CriticalNode(i,1)).Kt(CriticalNode(i,2)) = 0;
            node(CriticalNode(i,1)).Kphi(CriticalNode(i,2)) = 0;
            node(node(CriticalNode(i,1)).neighbors(CriticalNode(i,2))).Kn(node(node(CriticalNode(i,1)).neighbors(CriticalNode(i,2))).neighbors == CriticalNode(i,1)) = 0;
            node(node(CriticalNode(i,1)).neighbors(CriticalNode(i,2))).Kt(node(node(CriticalNode(i,1)).neighbors(CriticalNode(i,2))).neighbors == CriticalNode(i,1)) = 0;
            node(node(CriticalNode(i,1)).neighbors(CriticalNode(i,2))).Kphi(node(node(CriticalNode(i,1)).neighbors(CriticalNode(i,2))).neighbors == CriticalNode(i,1)) = 0;
            BondDeflect = node(CriticalNode(i,1)).u - node(node(CriticalNode(i,1)).neighbors(CriticalNode(i,2))).u;
            AvgDeflect = 0.5*(node(CriticalNode(i,1)).u + node(node(CriticalNode(i,1)).neighbors(CriticalNode(i,2))).u);
            BrokenTensionBond(TensionBondCount,:) = [CriticalNode(i,1:2),abs(yT), BondDeflect, node(CriticalNode(i,1)).u(3),AvgDeflect];
            TensionBondCount = TensionBondCount + 1; 
        end
    end
end  
            
        
    