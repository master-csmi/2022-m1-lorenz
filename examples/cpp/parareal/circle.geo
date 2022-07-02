h = DefineNumber[0.1, Name "Parameters/h"];

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, h};
Point(3) = {0, 1, 0, h};
Point(4) = {-1, 0, 0, h};
Point(5) = {0, -1, 0, h};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 2};

Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};

Physical Surface("Omega", 4) = {1};
Physical Curve("Dirichlet", 5) = {1};
Physical Curve("Neumann", 6) = {2};
Physical Curve("Fourier", 7) = {3};
