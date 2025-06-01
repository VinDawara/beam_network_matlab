Permanent variables
1. M --> Number of lattice points in vertical direction
2. N --> Number of lattice points in horizontal direction
3. lamba --> Scaling paramter
4. Kn --> Normal stiffness of the bond
5. TNRatio --> Tangential to normal stiffness' ratio
6. RNRatio --> Rotational to tangential siffness' ratio
7. node
8. porosity
9. TensionBondCount
10. BrokenTensionBond
11. delta
12. Force
13. Sigma
14. Strain
15. Energy
16. CN
17. loop
18. yB
19. yT



Program 1 - TriangularMesh
1. a --> lattice spacing
2. LeftBoundary --> Stores the IDs of the left boundaries of the lattice
3. RightBoundary --> Stores the IDs of the right boundaries of the lattice
4. TotalNodes --> Stores the number of total nodes
5. LowerBoundary
6. UpperBoundary
7. FixedNodes
8. AuxillaryNodes

Program 2 - PorousConfg
1. TotalPorousNode
2. PoreID

Program 3 - VoronoiEdges
1. NodePos
2. L, W

Program 4 - PoreDistribution
1. PoreLen
2. PoreDist
3. Poresize
4. PoreStructure
5. FinalVolume

Program 5 - LatticeView
1. g
2. PoreStructure

Program 6 - MatrixEquationSparse/MatrixAUpdate/Update
1. NoBondNode --> This variable is further used in solver
2. A, B

Program 7 - Solver
1. flag, relresi, itr, maxit, tol --> at the end, you can delete them

Program 8 - ForceComponent
1. FT, FB --> you can delete at the end
2. FaceArea

Program 9 - FractureCriterion
1. alpha, h, Area, W, CriticalTension
2. BondBreak --> at the end



Program 10 - SystemEnergy
1. SysEnergy --> you can delete at the end

Program 11 - CrackPath
1. L,W

Program 12 - NodeCN
1. NodeDeg
2. NodeDist










