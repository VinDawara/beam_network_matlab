%% STORING LATTICE INFORMATION
%% DATE: 30/04/2021

% Lattice spacing
a = 1;

% Node index
ID = 1;

% Left and Right boundary nodes for nodal stress calculation
LeftBoundary = [];
RightBoundary = [];
k = 1;
l = 1;
% First Auxillary node
node(ID).r = [-a, 0];
node(ID).r_initial = node(ID).r;
node(ID).neighbors = [2, N+2];
node(ID).Kn = [Kn, Kn];
node(ID).aph = [0 60];
node(ID).Kt = TNRatio*node(ID).Kn;
node(ID).Kphi = RTRatio*node(ID).Kt;
ID = ID + 1;
% Generating triangular lattice and defining its connectivity
for i = 1:M
    for j = 1:N
        if( i == 1)
            if(j ~= N)
                s = 1;
                Lexi = i;
                Lexj = j + 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                s = 2;
                Lexi = i + 1;
                Lexj = j + 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                s = 3;
                Lexi = i + 1;
                Lexj = j;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                s = 4;
                Lexi = i;
                Lexj = j - 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                node(ID).Kn = [Kn Kn Kn Kn];
                node(ID).aph = [0 60 120 180];
            else
                s = 1;
                Lexi = i + 1;
                Lexj = j;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                s = 2;
                Lexi = i;
                Lexj = j - 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
    
                node(ID).Kn = [Kn Kn];
                node(ID).aph = [120 180];
                RightBoundary(l) = ID;
                l = l + 1;
            end
        elseif(i == M) % Note: M should be even
            if(j == 1)
                s = 1;
                Lexi = i;
                Lexj = j + 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                s = 2;
                Lexi = i - 1;
                Lexj = j;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                
                node(ID).Kn = [Kn Kn];
                node(ID).aph = [0 -60];
                LeftBoundary(k) = ID;
                k = k + 1;
            else
                s = 1;
                Lexi = i;
                Lexj = j + 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                s = 2;
                Lexi = i;
                Lexj = j - 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;  
                s = 3;
                Lexi = i - 1;
                Lexj = j - 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                s = 4;
                Lexi = i - 1;
                Lexj = j;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;                
                node(ID).Kn = [Kn Kn Kn Kn];
                node(ID).aph = [0 180 -120 -60];
            end
        else
            if(j == 1)
                s = 1;
                Lexi = i;
                Lexj = j + 1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                if(mod(i,2) == 0)
                    s = 2;
                    Lexi = i + 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 3;
                    Lexi = i - 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    
                    node(ID).Kn = [Kn Kn Kn];
                    node(ID).aph = [0 60 -60]; 
                    LeftBoundary(k) = ID;
                    k = k + 1;
                else
                    s = 2;
                    Lexi = i + 1;
                    Lexj = j + 1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 3;
                    Lexi = i + 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 4;
                    Lexi = i - 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 5;
                    Lexi = i - 1;
                    Lexj = j + 1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    node(ID).Kn = [Kn Kn Kn Kn Kn];
                    node(ID).aph = [0 60 120 -120 -60];
                end
            elseif(j == N)
                if(mod(i,2) == 0)
                    s = 1;
                    Lexi = i + 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 2;
                    Lexi = i + 1;
                    Lexj = j - 1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 3;
                    Lexi = i;
                    Lexj = j-1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 4;
                    Lexi = i - 1;
                    Lexj = j - 1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 5;
                    Lexi = i - 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    node(ID).Kn = [Kn Kn Kn Kn Kn];
                    node(ID).aph = [60 120 180 -120 -60];
                else
                    s = 1;
                    Lexi = i + 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 2;
                    Lexi = i;
                    Lexj = j-1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 3;
                    Lexi = i - 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    
                    node(ID).Kn = [Kn Kn Kn];
                    node(ID).aph = [120 180 -120];
                    RightBoundary(l) = ID;
                    l = l + 1;
                end
            else
                s = 1;
                Lexi = i;
                Lexj = j+1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                s = 4;
                Lexi = i;
                Lexj = j-1;
                Lex = (Lexi-1)*N + Lexj + 1;
                node(ID).neighbors(s) = Lex;
                if(mod(i,2) == 0)
                    s = 2;
                    Lexi = i + 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 3;
                    Lexi = i + 1;
                    Lexj = j - 1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 5;
                    Lexi = i - 1;
                    Lexj = j - 1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 6;
                    Lexi = i - 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                else
                    s = 2;
                    Lexi = i + 1;
                    Lexj = j + 1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 3;
                    Lexi = i + 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 5;
                    Lexi = i - 1;
                    Lexj = j;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                    s = 6;
                    Lexi = i - 1;
                    Lexj = j + 1;
                    Lex = (Lexi-1)*N + Lexj + 1;
                    node(ID).neighbors(s) = Lex;
                end   
                node(ID).Kn = [Kn Kn Kn Kn Kn Kn];
                node(ID).aph = [0 60 120 180 -120 -60];
            end
        end
       node(ID).Kt = TNRatio*node(ID).Kn;
       node(ID).Kphi = RTRatio*node(ID).Kt;
       if(mod(i,2) == 0)
            node(ID).r = [((j-1)*a-0.5*a), (i-1)*a*0.5*sqrt(3)];
            node(ID).r_initial = node(ID).r;
        else
            node(ID).r = [(j-1)*a, (i-1)*a*0.5*sqrt(3)];
            node(ID).r_initial = node(ID).r;
        end
       ID = ID + 1;
    end                                                        
end
% End Auxillary node
s = 1;
Lexi = M;
Lexj = N;
Lex = (Lexi-1)*N + Lexj + 1;
node(ID).neighbors(s) = Lex;
s = 2;
Lexi = M-1;
Lexj = N;
Lex = (Lexi-1)*N + Lexj + 1;
node(ID).neighbors(s) = Lex;
node(ID).Kn = [Kn Kn];
node(ID).aph = [180 -120];
node(ID).Kt = TNRatio*node(ID).Kn;
node(ID).Kphi = RTRatio*node(ID).Kt;
i = M;
j = N + 1;
node(ID).r = [((j-1)*a-0.5*a), (i-1)*a*0.5*sqrt(3)];
node(ID).r_initial = node(ID).r;
TotalNodes = ID;

% Connecting first node in second row to the first Auxillary node
ID = N + 2;
s = 3;
Lexi = 1;
Lexj = 0;
Lex = (Lexi-1)*N + Lexj + 1;
node(ID).neighbors = [node(ID).neighbors(1:2), Lex, node(ID).neighbors(end)];
node(ID).aph = [node(ID).aph(1:2), -120, node(ID).aph(end)];
node(ID).Kn = [node(ID).Kn(1:2), Kn, node(ID).Kn(end)];
node(ID).Kt = TNRatio*node(ID).Kn;
node(ID).Kphi = RTRatio*node(ID).Kt;
% Connnecting last node in M-1 row with the end Auxillary node
Lexi = M-1;
Lexj = N;
Lex = (Lexi-1)*N + Lexj + 1;
ID = Lex;
Lexi = M;
Lexj = N + 1;
Lex = (Lexi-1)*N + Lexj + 1;
node(ID).neighbors = [Lex, node(ID).neighbors];
node(ID).aph = [60, node(ID).aph];
node(ID).Kn = [Kn,  node(ID).Kn];
node(ID).Kt = TNRatio*node(ID).Kn;
node(ID).Kphi = RTRatio*node(ID).Kt;

% Defining boundary nodes
BoundaryNodes = [1:N+1, TotalNodes-N:TotalNodes];
LowerBoundary = 1:N+1;
UpperBoundary = TotalNodes-N:TotalNodes;
FixedNodes = [];
AuxillaryNodes = [1 TotalNodes];

clrvars = {'Lexi','Lexj','Lex'};
clear(clrvars{:})


