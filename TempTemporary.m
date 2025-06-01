Temporary variables
-------------------

1. index
2. t, dim, str, i


Program 1 - TriangularMesh
1. ID --> node referncing number
2. k,l,i,j,l --> referencing index
3. Lexi, Lexj, Lex --> RI

clrvars = {'ID','i','j','l','s','Lexi','Lexj','Lex'};
clear(clrvars{:})

Program 2 - PorousConfg
1. RandNum
2. r
3. i,s

clrvars = {'RandNum','r','i','s'};
clear(clrvars{:})

Program 3 - VoronoiEdges
1. k
2. dx, dy
3. Arti_ID
4. ID, i, j
5. DT
6. Voro_Edges, Voro_Region
7. PolyVertices, poly, Linesegment, in
8. s, k
9. X,Y
10. P1, P2, intercept
11. temp1, temp2

clrvars = {'k','dx','dy','s', 'Arti_ID', 'ID', 'DT', 'Voro_Edges', 'Voro_Region','PolyVertices', 'X', 'Y', 'P1', 'P2', ...
    'intercept', 'temp1', 'temp2', 'poly'};
clear(clrvars{:})

Program 4 - PoreDistribution

clrvars = {'k', 'Link', 'temp', 'Poreneighbor','i'}
clear(clrvars{:})

Program 5 - LatticeView
clrvars = {'k', 'l', 'Pos', 'Network','ID','s','rowx','rowy'}

Program 6 - MatrixEquationSparse/MatrixAUpdate
clrvars = {'i', 'I', 'J', 'Entry','k','ID','s','l'}
clear(clrvars{:})

Program 7 - Solver
clrvars = {'k', 'ID', 'v'}

Program 8 - ForceComponents
clrvars = {'ID'}

Program 9 - ForceEstimation
clrvars = {'i'}

Program 10 - FractureCriterion
clrvars = {'CriticalNode','k','s','SigmaT','ID','maxStress','BondDeflect'}

Program 11 - SystemEnergy
clrvars = {'ID','s'}

Program 12 - CrackPath
clrvars = {'k', 'l', 'Pos', 'Network','ID','s','rowx','rowy','i'}

Program 13 - NodeCN
clrvars = {'k', 'l', 'Pos', 'Network','ID','s'}













