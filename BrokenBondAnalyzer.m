%% Analyzing broken bonds
%% Date: 11/04/2022

% 1. How many bonds brokens at given remote displacement?
% 2. What is nodal displacement distribution for the broken bond?

% About variable - 'BrokenTensionBond'
% Index 1 --> Node ID of the node to which broken bond is connected
% Index 2 --> Bond number of broken bond with respect to that node 
% Index 3 --> Remote load at which bond broke
% Index 4,5, 6 --> change in bond displacements and rotation
% Index 7 --> Node rotation
% Index 8,9,10 --> Average displacements and rotation of the nodes
% connecting the bond

NumBondBroken = zeros(length(Load_Step),1);
for i = 1:length(Load_Step)
    NumBondBroken(i) = nnz(BrokenTensionBond(:,3) == Load_Step(i));
end


u_values = sqrt(BrokenTensionBond(:,8).^2 + BrokenTensionBond(:,9).^2);
% Plot the histogram of 'u_values' to get the distribution

u_relative = sqrt(BrokenTensionBond(:,4).^2 + BrokenTensionBond(:,5).^2);
% Plot the histogram of 'u_relative' to get the distribution




    
    
