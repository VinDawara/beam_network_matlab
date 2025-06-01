%% Energy of the system with load
%% Date: 02/01/2022

% SysEnergy = 0;
SysEnergy_Axial = 0;
SysEnergy_Bending = 0;

for ID = 1:TotalNodes
    for s = 1:length(node(ID).neighbors)
        theta = node(ID).aph(s);
        Q = [cosd(theta), sind(theta), 0;-sind(theta), cosd(theta), 0;0, 0, 1];
        ui = Q'*node(ID).u';
        uj = Q'*node(node(ID).neighbors(s)).u';
        
        SysEnergy_Axial = SysEnergy_Axial + 0.5*(node(ID).Kn(s)*(ui(1) - uj(1))^2);
        SysEnergy_Bending = SysEnergy_Bending + 0.5*(node(ID).Kt(s)*(ui(2) - uj(2) + a*ui(3))^2 ...
            + node(ID).Kt(s)*a*(ui(3) - uj(3))*(ui(2) - uj(2) + a*ui(3)) ...
            + node(ID).Kphi(s)*(ui(3) - uj(3))^2);
        
%         SysEnergy_Axial = SysEnergy_Axial + 0.5*(node(ID).Kn(s)*(node(ID).u(1) - node(node(ID).neighbors(s)).u(1))^2);
%         SysEnergy_Bending = SysEnergy_Bending + 0.5*(node(ID).Kt(s)*(node(ID).u(2) - node(node(ID).neighbors(s)).u(2) + a*node(ID).u(3))^2 ...
%             + node(ID).Kt(s)*a*(node(ID).u(3) - node(node(ID).neighbors(s)).u(3))*(node(ID).u(2) - node(node(ID).neighbors(s)).u(2) + a*node(ID).u(3)) ...
%             + node(ID).Kphi(s)*(node(ID).u(3) - node(node(ID).neighbors(s)).u(3))^2);
     
    end
end

% SysEnergy = 0.5*SysEnergy;
SysEnergy_Axial = 0.5*SysEnergy_Axial ;
SysEnergy_Bending = 0.5*SysEnergy_Bending;
SysEnergy = SysEnergy_Axial + SysEnergy_Bending;