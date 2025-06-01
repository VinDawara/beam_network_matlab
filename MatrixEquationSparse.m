%% MATRIX EQUATION
%% DATE: 11/08/2021

NoBondNode = [];    % Store the ID of the node that either have no bonds
                    % or prescribed all the displacment components (u, v,
                    % phi)
% k = 1;              % Increment the NoBondNode variable

% Apply the boundary conditions to upper and lower boundary nodes
for i = 1:length(LowerBoundary)
    node(LowerBoundary(i)).u(2) = yB;
end

for i = 1:length(UpperBoundary)
    node(UpperBoundary(i)).u(2) = yT;
end

% A = sparse(3*TotalNodes,3*TotalNodes,TotalNodes);
I = zeros(TotalNodes,1);
J = zeros(TotalNodes,1);
Entry = zeros(TotalNodes,1);
k = 0;
B = zeros(3*TotalNodes,1);

for ID = 1:TotalNodes
    if(~any(ID == BoundaryNodes))
        % Coefficient of u_i, v_i and phi_i in the X force equilibrium
        % equation
        k = k + 1;
        I(k) = 3*ID-2;
        J(k) = 3*ID-2;
        Entry(k) = sum(node(ID).Kn.*(cosd(node(ID).aph).^2) + node(ID).Kt.*(sind(node(ID).aph).^2));
        k = k + 1;
        I(k) = 3*ID-2;
        J(k) = 3*ID-1;
        Entry(k) = sum((node(ID).Kn - node(ID).Kt).*sind(node(ID).aph).*cosd(node(ID).aph));
        k = k + 1;
        I(k) = 3*ID-2;
        J(k) = 3*ID;
        Entry(k) = -0.5*a*sum(node(ID).Kt.*sind(node(ID).aph));
        
        
        % Coefficient of u_i, v_i and phi_i in the Y force equilibrium
        % equation
        k = k + 1;
        I(k) = 3*ID-1;
        J(k) = 3*ID-2;
        Entry(k) = sum((node(ID).Kn - node(ID).Kt).*sind(node(ID).aph).*cosd(node(ID).aph));
        k = k + 1;
        I(k) = 3*ID-1;
        J(k) = 3*ID-1;
        Entry(k) = sum(node(ID).Kn.*(sind(node(ID).aph).^2) + node(ID).Kt.*(cosd(node(ID).aph).^2));
        k = k + 1;
        I(k) = 3*ID-1;
        J(k) = 3*ID;
        Entry(k) = 0.5*a*sum(node(ID).Kt.*cosd(node(ID).aph));
        
        
        % Coefficient of u_i, v_i and phi_i in the Moment equilibrium
        % equation
        k = k + 1;
        I(k) = 3*ID;
        J(k) = 3*ID-2;
        Entry(k) = -0.5*a*sum(node(ID).Kt.*sind(node(ID).aph));
        k = k + 1;
        I(k) = 3*ID;
        J(k) = 3*ID-1;
        Entry(k) = 0.5*a*sum(node(ID).Kt.*cosd(node(ID).aph));
        k = k + 1;
        I(k) = 3*ID;
        J(k) = 3*ID;
        Entry(k) = a^2*sum(node(ID).Kphi);
        
        
       % Cofficients of u_j, v_j and phi_j in X,Y force and moment
       % equilibrium equations 
        for s = 1:length(node(ID).neighbors)
            if(~any(node(ID).neighbors(s) == BoundaryNodes))
                % Coefficient of u_j, v_j and phi_j in the X force equilibrium
                % equation
                k = k + 1;
                I(k) = 3*ID-2;
                J(k) = 3*node(ID).neighbors(s)-2;
                Entry(k) = - node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) - node(ID).Kt(s)*(sind(node(ID).aph(s))^2);
                Entry(k) = round(Entry(k),12);
                k = k + 1;
                I(k) = 3*ID-2;
                J(k) = 3*node(ID).neighbors(s)-1;
                Entry(k) = -(node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s));
                Entry(k) = round(Entry(k),12);
                k = k + 1;
                I(k) = 3*ID-2;
                J(k) = 3*node(ID).neighbors(s);
                Entry(k) = -0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                Entry(k) = round(Entry(k),12);
                
                
                % Coefficient of u_j, v_j and phi_j in the Y force equilibrium
                % equation
                k = k + 1;
                I(k) = 3*ID-1;
                J(k) = 3*node(ID).neighbors(s)-2;
                Entry(k) = - (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s));
                Entry(k) = round(Entry(k),12);
                k = k + 1;
                I(k) = 3*ID-1;
                J(k) = 3*node(ID).neighbors(s)-1;
                Entry(k) = - node(ID).Kn(s)*(sind(node(ID).aph(s))^2) - node(ID).Kt(s)*(cosd(node(ID).aph(s))^2);
                Entry(k) = round(Entry(k),12);
                k = k + 1;
                I(k) = 3*ID-1;
                J(k) = 3*node(ID).neighbors(s);
                Entry(k) =  0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s));
                Entry(k) = round(Entry(k),12);
                
                
                % Coefficient of u_j, v_j and phi_j in the moment equilibrium
                % equation
                k = k + 1;
                I(k) = 3*ID;
                J(k) = 3*node(ID).neighbors(s)-2;
                Entry(k) = 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                Entry(k) = round(Entry(k),12);
                k = k + 1;
                I(k) = 3*ID;
                J(k) = 3*node(ID).neighbors(s)-1;
                Entry(k) = - 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s));
                Entry(k) = round(Entry(k),12);
                k = k + 1;
                I(k) = 3*ID;
                J(k) = 3*node(ID).neighbors(s);
                Entry(k) =  0.5*a^2*node(ID).Kphi(s);
                Entry(k) = round(Entry(k),12);
                

            else
                if(any(node(ID).neighbors(s) == FixedNodes))  
                    % RHS values for X equations
                    B(3*ID-2) = B(3*ID-2) + (node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) + node(ID).Kt(s)*(sind(node(ID).aph(s))^2))*node(node(ID).neighbors(s)).u(1);
                    B(3*ID-2) = B(3*ID-2) + (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                    B(3*ID-2) = B(3*ID-2) + 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s))*node(node(ID).neighbors(s)).u(3);
                    
                    % RHS values for Y equations
                    B(3*ID-1) = B(3*ID-1) + (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(1);
                    B(3*ID-1) = B(3*ID-1) + (node(ID).Kn(s)*(sind(node(ID).aph(s))^2) + node(ID).Kt(s)*(cosd(node(ID).aph(s))^2))*node(node(ID).neighbors(s)).u(2);
                    B(3*ID-1) = B(3*ID-1) - 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(3);
                    
                    % RHS values for moment equations
                    B(3*ID) = B(3*ID) - 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s))*node(node(ID).neighbors(s)).u(1);
                    B(3*ID) = B(3*ID) + 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                    B(3*ID) = B(3*ID) - 0.5*a^2*node(ID).Kphi(s)*node(node(ID).neighbors(s)).u(3);
                    
                else
                    % Coefficient of u_j and phi_j in the X force equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    Entry(k) = - node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) - node(ID).Kt(s)*(sind(node(ID).aph(s))^2);
                    Entry(k) = round(Entry(k),12);
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s);
                    Entry(k) =  -0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                    Entry(k) = round(Entry(k),12);
                    
                
                    % Coefficient of u_j and phi_j in the Y force equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID-1;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    Entry(k) = - (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s));
                    Entry(k) = round(Entry(k),12);
                    k = k + 1;
                    I(k) = 3*ID-1;
                    J(k) = 3*node(ID).neighbors(s);
                    Entry(k) =  0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s));
                    Entry(k) = round(Entry(k),12);
                    
                
                    % Coefficient of u_j and phi_j in the moment equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    Entry(k) = 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                    Entry(k) = round(Entry(k),12);
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s);
                    Entry(k) =  0.5*a^2*node(ID).Kphi(s);
                    Entry(k) = round(Entry(k),12);
 
                    % RHS v_j values for X equations
                    B(3*ID-2) = B(3*ID-2) + (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                    
                    % RHS v_j values for Y equations
                    B(3*ID-1) = B(3*ID-1) + (node(ID).Kn(s)*(sind(node(ID).aph(s))^2) + node(ID).Kt(s)*(cosd(node(ID).aph(s))^2))*node(node(ID).neighbors(s)).u(2);
                    
                    % RHS values for moment equations
                    B(3*ID) = B(3*ID) + 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);

                end
            end
        end
    else
        if(~any(ID == FixedNodes))
            % Coefficient of u_i, and phi_i in the X force equilibrium
            % equation
            k = k + 1;
            I(k) = 3*ID-2;
            J(k) = 3*ID-2;
            Entry(k) = sum(node(ID).Kn.*(cosd(node(ID).aph).^2) + node(ID).Kt.*(sind(node(ID).aph).^2));
            k = k + 1;
            I(k) = 3*ID-2;
            J(k) = 3*ID;
            Entry(k) =  -0.5*a*sum(node(ID).Kt.*sind(node(ID).aph));
        
            % Coefficient of u_i, and phi_i in the Moment equilibrium
            % equation
            k = k + 1;
            I(k) = 3*ID;
            J(k) = 3*ID-2;
            Entry(k) = -0.5*a*sum(node(ID).Kt.*sind(node(ID).aph));
            k = k + 1;
            I(k) = 3*ID;
            J(k) = 3*ID;
            Entry(k) =  a^2*sum(node(ID).Kphi);
           
            
            % RHS (Term with v_i) in the X equations
            B(3*ID-2) = B(3*ID-2) - sum((node(ID).Kn - node(ID).Kt).*sind(node(ID).aph).*cosd(node(ID).aph))*node(ID).u(2);
            
            % RHS (Term with v_i) in the moment equations
            B(3*ID) = B(3*ID) - 0.5*a*sum(node(ID).Kt.*cosd(node(ID).aph))*node(ID).u(2);
        
            for s = 1:length(node(ID).neighbors)
                if(~any(node(ID).neighbors(s) == BoundaryNodes))
                    % Coefficient of u_j, v_j and phi_j in the X force equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    Entry(k) = - node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) - node(ID).Kt(s)*(sind(node(ID).aph(s))^2);
                    Entry(k) = round(Entry(k),12);
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s)-1;
                    Entry(k) = - (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s));
                    Entry(k) = round(Entry(k),12);
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s);
                    Entry(k) =  -0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                    Entry(k) = round(Entry(k),12);
                    
                
                    % Coefficient of u_j, v_j and phi_j in the moment equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    Entry(k) = 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                    Entry(k) = round(Entry(k),12);
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s)-1;
                    Entry(k) = - 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s));
                    Entry(k) = round(Entry(k),12);
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s);
                    Entry(k) = 0.5*a^2*node(ID).Kphi(s);
                    Entry(k) = round(Entry(k),12);
                     
                else
                    if(any(node(ID).neighbors(s) == FixedNodes))
                        % RHS values for X equations
                        B(3*ID-2) = B(3*ID-2) + (node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) + node(ID).Kt(s)*(sind(node(ID).aph(s))^2))*node(node(ID).neighbors(s)).u(1);
                        B(3*ID-2) = B(3*ID-2) + (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                        B(3*ID-2) = B(3*ID-2) + 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s))*node(node(ID).neighbors(s)).u(3);
                    
                        % RHS values for moment equations
                        B(3*ID) = B(3*ID) - 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s))*node(node(ID).neighbors(s)).u(1);
                        B(3*ID) = B(3*ID) + 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                        B(3*ID) = B(3*ID) - 0.5*a^2*node(ID).Kphi(s)*node(node(ID).neighbors(s)).u(3);
                                
                    else
                        % Coefficient of u_j and phi_j in the X force equilibrium
                        % equation
                        k = k + 1;
                        I(k) = 3*ID-2;
                        J(k) = 3*node(ID).neighbors(s)-2;
                        Entry(k) = - node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) - node(ID).Kt(s)*(sind(node(ID).aph(s))^2);
                        Entry(k) = round(Entry(k),12);
                        k = k + 1;
                        I(k) = 3*ID-2;
                        J(k) = 3*node(ID).neighbors(s);
                        Entry(k) = -0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                        Entry(k) = round(Entry(k),12);
                        
                
                        % Coefficient of u_j and phi_j in the moment equilibrium
                        % equation
                        k = k + 1;
                        I(k) = 3*ID;
                        J(k) = 3*node(ID).neighbors(s)-2;
                        Entry(k) = 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                        Entry(k) = round(Entry(k),12);
                        k = k + 1;
                        I(k) = 3*ID;
                        J(k) = 3*node(ID).neighbors(s);
                        Entry(k) = 0.5*a^2*node(ID).Kphi(s);
                        Entry(k) = round(Entry(k),12);
 
                        % RHS v_j values for X equations
                        B(3*ID-2) = B(3*ID-2) + (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                    
                        % RHS values for moment equations
                        B(3*ID) = B(3*ID) + 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                    end       
                end
            end      
        end    
    end
%     if(sum(abs(A(3*ID-2,:))) == 0 && sum(abs(A(:,3*ID-2))) == 0)
%        NoBondNode(k) = ID;
%        k = k + 1;
%     end
end

A = sparse(I,J,Entry,3*TotalNodes,3*TotalNodes);


    