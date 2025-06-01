%% FUNCTION TO ARRANGE THE VERTICES OF THE POLYGON IN ANIT-CLOCKWISE SENSE
%% FROM 0 TO 360 DEGREE
%% DATE: 17/04/2021

function [Vertices] = AntiClockwise(Vertices)
theta = zeros(length(Vertices(:,1)),1);
xc = mean([min(Vertices(:,1)), max(Vertices(:,1))]);
yc = mean([min(Vertices(:,2)), max(Vertices(:,2))]);

for i = 1:length(Vertices)
    vector = [(Vertices(i,1) - xc),(Vertices(i,2) - yc)];
    if(vector(1) >= 0 && vector(2) >= 0)
        theta(i) = atand(vector(2)/vector(1));
    elseif(vector(1) < 0 && vector(2) > 0)
        theta(i) = 180 - atand(vector(2)/abs(vector(1)));
    elseif(vector(1) < 0 && vector(2) <= 0)
        theta(i) = 180 + atand(abs(vector(2))/abs(vector(1)));
    else
        theta(i) = 360 - atand(abs(vector(2))/vector(1));
    end
end

A = [theta, Vertices];
A = sortrows(A);
Vertices = A(:,2:3);
end
    
     