%% NODAL STRESS COMPUTATION 
%% DATE: 21/10/2021

for ID = 1:TotalNodes
    node(ID).Stress = zeros(2,2);
    % Construct horizontal plane for syx and syy
    if(any(ID == LowerBoundary))
        theta = 180;
    else
        theta = 0;
    end
    RN = [];        % Neighbors that lie right side of the plane
    Fr = [];        % Radial force on bonds
    Ft = [];        % Tangential force
    A = [];         % Projected area of the beam on the plane
    
    % Variable to store angle of the bonds measured anticlockwise sense from bond 1
    AntiNodeAngle = zeros(length(node(ID).neighbors),1);
    k = 1;  % index updater
    
    % Tangential and normal vector of the plane
    t = [cosd(theta), sind(theta)];
    n = [-sind(theta), cosd(theta)];
    
    An = 0; % Total projected area in the plane
    
    syx = 0;
    syy = 0;
    %%-----LOOP FOR EACH BOND--------%%
    for s = 1:length(node(ID).neighbors)
        
        % Measuring angle anticlockwise w.r.t bond 1
        if(node(ID).aph(s) >= 0)
            AntiNodeAngle(s) = node(ID).aph(s);
        else
            AntiNodeAngle(s) = 360 + node(ID).aph(s);
        end
        
        if(theta <= 180)
            % Checking whether bond lies inside the cut plane
            if(AntiNodeAngle(s) <= theta || AntiNodeAngle(s) >= theta + 180)
                RN(k) = node(ID).neighbors(s);
                er = [cosd(AntiNodeAngle(s)), sind(AntiNodeAngle(s))];
                et = [-sind(AntiNodeAngle(s)), cosd(AntiNodeAngle(s))];
                uij = node(ID).u(1:2) - node(node(ID).neighbors(s)).u(1:2);
                Fr(k,:) = lambda*(node(ID).Kn(s)*uij*er')*er;
                phi = node(node(ID).neighbors(s)).u(3) + node(ID).u(3);
                Ft(k,:) = lambda*node(ID).Kt(s)*(uij*et' + 0.5*a*phi)*et;
                A(k) = Area*er*n';
                if(A(k) ~= 0)
                    syx = syx + Fr(k,:)*t' + Ft(k,:)*t';
                    syy= syy + Fr(k,:)*n' + Ft(k,:)*n';
                    An = An + abs(A(k));
                end
                k = k + 1;
            end
        end
    end
    if(An ~= 0)
        syx = syx/An;
        syy = syy/An;
    end
    
    % Construct horizontal plane for sxy and sxx
    sxy = 0;
    sxx = 0;
    theta = 90;
    
    RN = [];        % Neighbors that lie right side of the plane
    Fr = [];        % Radial force on bonds
    Ft = [];        % Tangential force
    A = [];         % Projected area of the beam on the plane
    
    % Variable to store angle of the bonds measured anticlockwise sense from bond 1
    AntiNodeAngle = zeros(length(node(ID).neighbors),1);
    k = 1;  % index updater
    
    % Tangential and normal vector of the plane
    t = [cosd(theta), sind(theta)];
    n = [-sind(theta), cosd(theta)];
    
    An = 0; % Total projected area in the plane
    
    for s = 1:length(node(ID).neighbors)
        
        % Measuring angle anticlockwise w.r.t bond 1
        if(node(ID).aph(s) >= 0)
            AntiNodeAngle(s) = node(ID).aph(s);
        else
            AntiNodeAngle(s) = 360 + node(ID).aph(s);
        end
        
        if(theta <= 180)
            % Checking whether bond lies inside the cut plane
            if(AntiNodeAngle(s) <= theta || AntiNodeAngle(s) >= theta + 180)
                RN(k) = node(ID).neighbors(s);
                er = [cosd(AntiNodeAngle(s)), sind(AntiNodeAngle(s))];
                et = [-sind(AntiNodeAngle(s)), cosd(AntiNodeAngle(s))];
                uij = node(ID).u(1:2) - node(node(ID).neighbors(s)).u(1:2);
                Fr(k,:) = (node(ID).Kn(s)*uij*er')*er;
                phi = node(node(ID).neighbors(s)).u(3) + node(ID).u(3);
                Ft(k,:) = node(ID).Kt(s)*(uij*et' + 0.5*a*phi)*et;
                A(k) = Area*er*n';
                if(A(k) ~= 0)
                    sxy = sxy + Fr(k,:)*t' + Ft(k,:)*t';
                    sxx= sxx + Fr(k,:)*n' + Ft(k,:)*n';
                    An = An + abs(A(k));
                end
                k = k + 1;
            end
        end
    end
    if(An ~= 0)
        sxy = sxy/An;
        sxx = sxx/An;
    end
    if(any(ID == BoundaryNodes))
        node(ID).Stress = [sxx, syx; syx, syy];
    else
        node(ID).Stress = [sxx, syx; syx, syy];
    end
end





    
    
    