%% PROGRAM FOR UPDATING ONLY BROKEN ENTRIES IN THE MATRIX 
%% DATE: 23/09/2021

if(~isempty(BrokenTensionBond))
NoBondNode = [];    % Store the ID of the node that either have no bonds
                    % or prescribed all the displacment components (u, v,
                    % phi)
% k = 1;              % Increment the NoBondNode variable

l = length(BrokenTensionBond(:,1));
I = zeros(l,1);
J = zeros(l,1);
k = 0;

for i = 1:l
    ID = BrokenTensionBond(i,1);
    if(~any(ID == BoundaryNodes))
        % Coefficient of u_i, v_i and phi_i in the X force equilibrium
        % equation
        k = k + 1;
        I(k) = 3*ID-2;
        J(k) = 3*ID-2;
        A(I(k),J(k)) = sum(node(ID).Kn.*(cosd(node(ID).aph).^2) + node(ID).Kt.*(sind(node(ID).aph).^2));
        k = k + 1;
        I(k) = 3*ID-2;
        J(k) = 3*ID-1;
        A(I(k),J(k)) = sum((node(ID).Kn - node(ID).Kt).*sind(node(ID).aph).*cosd(node(ID).aph));
        k = k + 1;
        I(k) = 3*ID-2;
        J(k) = 3*ID;
        A(I(k),J(k)) = -0.5*a*sum(node(ID).Kt.*sind(node(ID).aph));
        
        
        % Coefficient of u_i, v_i and phi_i in the Y force equilibrium
        % equation
        k = k + 1;
        I(k) = 3*ID-1;
        J(k) = 3*ID-2;
        A(I(k),J(k)) = sum((node(ID).Kn - node(ID).Kt).*sind(node(ID).aph).*cosd(node(ID).aph));
        k = k + 1;
        I(k) = 3*ID-1;
        J(k) = 3*ID-1;
        A(I(k),J(k)) = sum(node(ID).Kn.*(sind(node(ID).aph).^2) + node(ID).Kt.*(cosd(node(ID).aph).^2));
        k = k + 1;
        I(k) = 3*ID-1;
        J(k) = 3*ID;
        A(I(k),J(k)) = 0.5*a*sum(node(ID).Kt.*cosd(node(ID).aph));
        
        
        % Coefficient of u_i, v_i and phi_i in the Moment equilibrium
        % equation
        k = k + 1;
        I(k) = 3*ID;
        J(k) = 3*ID-2;
        A(I(k),J(k)) = -0.5*a*sum(node(ID).Kt.*sind(node(ID).aph));
        k = k + 1;
        I(k) = 3*ID;
        J(k) = 3*ID-1;
        A(I(k),J(k)) = 0.5*a*sum(node(ID).Kt.*cosd(node(ID).aph));
        k = k + 1;
        I(k) = 3*ID;
        J(k) = 3*ID;
        A(I(k),J(k)) = a^2*sum(node(ID).Kphi);
        
        
       % Cofficients of u_j, v_j and phi_j in X,Y force and moment
       % equilibrium equations 
        for s = 1:length(node(ID).neighbors)
            if(~any(node(ID).neighbors(s) == BoundaryNodes))
                % Coefficient of u_j, v_j and phi_j in the X force equilibrium
                % equation
                k = k + 1;
                I(k) = 3*ID-2;
                J(k) = 3*node(ID).neighbors(s)-2;
                A(I(k),J(k)) = - node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) - node(ID).Kt(s)*(sind(node(ID).aph(s))^2);
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                k = k + 1;
                I(k) = 3*ID-2;
                J(k) = 3*node(ID).neighbors(s)-1;
                A(I(k),J(k)) = -(node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s));
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                k = k + 1;
                I(k) = 3*ID-2;
                J(k) = 3*node(ID).neighbors(s);
                A(I(k),J(k)) = -0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                
                
                % Coefficient of u_j, v_j and phi_j in the Y force equilibrium
                % equation
                k = k + 1;
                I(k) = 3*ID-1;
                J(k) = 3*node(ID).neighbors(s)-2;
                A(I(k),J(k)) = - (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s));
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                k = k + 1;
                I(k) = 3*ID-1;
                J(k) = 3*node(ID).neighbors(s)-1;
                A(I(k),J(k)) = - node(ID).Kn(s)*(sind(node(ID).aph(s))^2) - node(ID).Kt(s)*(cosd(node(ID).aph(s))^2);
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                k = k + 1;
                I(k) = 3*ID-1;
                J(k) = 3*node(ID).neighbors(s);
                A(I(k),J(k)) =  0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s));
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                
                
                % Coefficient of u_j, v_j and phi_j in the moment equilibrium
                % equation
                k = k + 1;
                I(k) = 3*ID;
                J(k) = 3*node(ID).neighbors(s)-2;
                A(I(k),J(k)) = 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                k = k + 1;
                I(k) = 3*ID;
                J(k) = 3*node(ID).neighbors(s)-1;
                A(I(k),J(k)) = - 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s));
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                k = k + 1;
                I(k) = 3*ID;
                J(k) = 3*node(ID).neighbors(s);
                A(I(k),J(k)) =  0.5*a^2*node(ID).Kphi(s);
                A(I(k),J(k)) = round(A(I(k),J(k)),12);
                

            else
                if(~any(node(ID).neighbors(s) == FixedNodes))  
                    % Coefficient of u_j and phi_j in the X force equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    A(I(k),J(k)) = - node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) - node(ID).Kt(s)*(sind(node(ID).aph(s))^2);
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s);
                    A(I(k),J(k)) =  -0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    
                
                    % Coefficient of u_j and phi_j in the Y force equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID-1;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    A(I(k),J(k)) = - (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s));
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    k = k + 1;
                    I(k) = 3*ID-1;
                    J(k) = 3*node(ID).neighbors(s);
                    A(I(k),J(k)) =  0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s));
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    
                
                    % Coefficient of u_j and phi_j in the moment equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    A(I(k),J(k)) = 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s);
                    A(I(k),J(k)) =  0.5*a^2*node(ID).Kphi(s);
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);

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
            A(I(k),J(k)) = sum(node(ID).Kn.*(cosd(node(ID).aph).^2) + node(ID).Kt.*(sind(node(ID).aph).^2));
            k = k + 1;
            I(k) = 3*ID-2;
            J(k) = 3*ID;
            A(I(k),J(k)) =  -0.5*a*sum(node(ID).Kt.*sind(node(ID).aph));
        
            % Coefficient of u_i, and phi_i in the Moment equilibrium
            % equation
            k = k + 1;
            I(k) = 3*ID;
            J(k) = 3*ID-2;
            A(I(k),J(k)) = -0.5*a*sum(node(ID).Kt.*sind(node(ID).aph));
            k = k + 1;
            I(k) = 3*ID;
            J(k) = 3*ID;
            A(I(k),J(k)) =  a^2*sum(node(ID).Kphi);
           
        
            for s = 1:length(node(ID).neighbors)
                if(~any(node(ID).neighbors(s) == BoundaryNodes))
                    % Coefficient of u_j, v_j and phi_j in the X force equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    A(I(k),J(k)) = - node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) - node(ID).Kt(s)*(sind(node(ID).aph(s))^2);
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s)-1;
                    A(I(k),J(k)) = - (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s));
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    k = k + 1;
                    I(k) = 3*ID-2;
                    J(k) = 3*node(ID).neighbors(s);
                    A(I(k),J(k)) =  -0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    
                
                    % Coefficient of u_j, v_j and phi_j in the moment equilibrium
                    % equation
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s)-2;
                    A(I(k),J(k)) = 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s)-1;
                    A(I(k),J(k)) = - 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s));
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                    k = k + 1;
                    I(k) = 3*ID;
                    J(k) = 3*node(ID).neighbors(s);
                    A(I(k),J(k)) = 0.5*a^2*node(ID).Kphi(s);
                    A(I(k),J(k)) = round(A(I(k),J(k)),12);
                     
                else
                    if(~any(node(ID).neighbors(s) == FixedNodes))         
                        % Coefficient of u_j and phi_j in the X force equilibrium
                        % equation
                        k = k + 1;
                        I(k) = 3*ID-2;
                        J(k) = 3*node(ID).neighbors(s)-2;
                        A(I(k),J(k)) = - node(ID).Kn(s)*(cosd(node(ID).aph(s))^2) - node(ID).Kt(s)*(sind(node(ID).aph(s))^2);
                        A(I(k),J(k)) = round(A(I(k),J(k)),12);
                        k = k + 1;
                        I(k) = 3*ID-2;
                        J(k) = 3*node(ID).neighbors(s);
                        A(I(k),J(k)) = -0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                        A(I(k),J(k)) = round(A(I(k),J(k)),12);
                        
                
                        % Coefficient of u_j and phi_j in the moment equilibrium
                        % equation
                        k = k + 1;
                        I(k) = 3*ID;
                        J(k) = 3*node(ID).neighbors(s)-2;
                        A(I(k),J(k)) = 0.5*a*node(ID).Kt(s)*sind(node(ID).aph(s));
                        A(I(k),J(k)) = round(A(I(k),J(k)),12);
                        k = k + 1;
                        I(k) = 3*ID;
                        J(k) = 3*node(ID).neighbors(s);
                        A(I(k),J(k)) = 0.5*a^2*node(ID).Kphi(s);
                        A(I(k),J(k)) = round(A(I(k),J(k)),12);
 
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
end






    