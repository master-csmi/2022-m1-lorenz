SetFactory("OpenCASCADE");
h = 0.1

Point(1) = {1, 0, 0, h};
Point(2) = {0, 1, 0, h};
Point(3) = {-1, 0, 0, h};
Point(4) = {0, -1, 0, h};

Circle(1) = {0, 0, 0, 1, 0, 2*Pi};

Curve Loop(1) = {1};
Plane Surface(1) = {1};

Physical Curve("Dirichlet", 2) = {1};
Physical Surface("Omega", 3) = {1};
