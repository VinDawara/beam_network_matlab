%% ESTIMATION OF THE REMOTE LOADING
%% DATE: 01/05/2021

FT = [0 0];
FB = [0 0];
FaceArea = 0;
h = sqrt(TNRatio)*a;
Area = h^2;
% W = h^3/6;
for i = 1:length(LowerBoundary)
    for s = 1: length(node(LowerBoundary(i)).neighbors)
        FB(1) = FB(1) - (node(LowerBoundary(i)).Fn(s)*cosd(node(LowerBoundary(i)).aph(s)) - ...
            node(LowerBoundary(i)).Ft(s)*sind(node(LowerBoundary(i)).aph(s)));
        FB(2) = FB(2) - (node(LowerBoundary(i)).Fn(s)*sind(node(LowerBoundary(i)).aph(s)) + ...
            node(LowerBoundary(i)).Ft(s)*cosd(node(LowerBoundary(i)).aph(s)));
    end
end
for i = 1:length(UpperBoundary)
    for s = 1: length(node(UpperBoundary(i)).neighbors)
        FaceArea = FaceArea + Area*abs(sind(node(UpperBoundary(i)).aph(s)));
        FT(1) = FT(1) - (node(UpperBoundary(i)).Fn(s)*cosd(node(UpperBoundary(i)).aph(s)) - ...
            node(UpperBoundary(i)).Ft(s)*sind(node(UpperBoundary(i)).aph(s)));
        FT(2) = FT(2) - (node(UpperBoundary(i)).Fn(s)*sind(node(UpperBoundary(i)).aph(s)) + ...
            node(UpperBoundary(i)).Ft(s)*cosd(node(UpperBoundary(i)).aph(s)));
    end
end
SigT = FT/FaceArea;
SigB = FB/FaceArea;
