%% SOLVING Av = B using CGS
%% DATE = 13/08/2021


if(~isempty(NoBondNode))
    ZeroIndex = [3*NoBondNode-2; 3*NoBondNode-2; 3*NoBondNode];
    A(ZeroIndex,:) = [];
    A(:,ZeroIndex) = [];
    B(ZeroIndex,:) = [];
end
maxit = 10000;
tol = 1e-8;
[v, flag, relresi, itr] = pcg(A,B,tol,maxit);

k = 1;
for ID = 1:TotalNodes
    if(~any(ID == BoundaryNodes))
        if(~any(ID == NoBondNode))
            node(ID).u = [v(k), v(k+1), v(k+2)];
            k = k + 3;
        else
            node(ID).u = [0 0 0];
        end
    else
        if(~any(ID == FixedNodes))
            if(~any(ID == NoBondNode))
                node(ID).u = [v(k), node(ID).u(2), v(k+2)];
                k = k + 3;
            else
                node(ID).u = [0 0 0];
            end
        end
    end
    node(ID).r = node(ID).r + node(ID).u(1:2);
end


