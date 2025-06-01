%% Update B
%% Date: 23/09/2021

% Apply the boundary conditions to upper and lower boundary nodes
for i = 1:length(LowerBoundary)
    node(LowerBoundary(i)).u(2) = yB;
end

for i = 1:length(UpperBoundary)
    node(UpperBoundary(i)).u(2) = yT;
end

B = zeros(3*TotalNodes,1);

for ID = 1:TotalNodes
    if(~any(ID == BoundaryNodes))     
        for s = 1:length(node(ID).neighbors)
            if(any(node(ID).neighbors(s) == BoundaryNodes))
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
            % RHS (Term with v_i) in the X equations
            B(3*ID-2) = B(3*ID-2) - sum((node(ID).Kn - node(ID).Kt).*sind(node(ID).aph).*cosd(node(ID).aph))*node(ID).u(2);
            
            % RHS (Term with v_i) in the moment equations
            B(3*ID) = B(3*ID) - 0.5*a*sum(node(ID).Kt.*cosd(node(ID).aph))*node(ID).u(2);
        
            for s = 1:length(node(ID).neighbors)
                if(any(node(ID).neighbors(s) == BoundaryNodes))
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
                        % RHS v_j values for X equations
                        B(3*ID-2) = B(3*ID-2) + (node(ID).Kn(s) - node(ID).Kt(s))*sind(node(ID).aph(s))*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                    
                        % RHS values for moment equations
                        B(3*ID) = B(3*ID) + 0.5*a*node(ID).Kt(s)*cosd(node(ID).aph(s))*node(node(ID).neighbors(s)).u(2);
                    end       
                end
            end      
        end    
    end
end

