h = DefineNumber[ 0.1, Name "Parameters/h" ];
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {0, 0, 2.7, h};
//+
Point(3) = {0, 8, 2.7, h};
//+
Point(4) = {0, 8, 0, h};
//+
Point(5) = {-6, 0, 0, h};
//+
Point(6) = {-6, 0, 2.7, h};
//+
Point(7) = {-6, 8, 2.7, h};
//+
Point(8) = {-6, 8, 0, h};
//+
Point(9) = {0, 0.5, 0.2, h};
//+
Point(10) = {0, 0.5, 2.2, h};
//+
Point(11) = {0, 3.5, 2.2, h};
//+
Point(12) = {0, 3.5, 0.2, h};
//+
Point(13) = {0, 4.5, 0.2, h};
//+
Point(14) = {0, 4.5, 2.2, h};
//+
Point(15) = {0, 7.5, 2.2, h};
//+
Point(16) = {0, 7.5, 0.2, h};//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {11, 10};
//+
Line(6) = {10, 9};
//+
Line(7) = {9, 12};
//+
Line(8) = {12, 11};
//+
Line(9) = {14, 13};
//+
Line(10) = {13, 16};
//+
Line(11) = {16, 15};
//+
Line(12) = {15, 14};
//+
Line(13) = {3, 7};
//+
Line(14) = {7, 8};
//+
Line(15) = {8, 4};
//+
Line(16) = {7, 6};
//+
Line(17) = {6, 5};
//+
Line(18) = {5, 8};
//+
Line(19) = {2, 6};
//+
Line(20) = {5, 1};
//+
Curve Loop(1) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {12, 9, 10, 11};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 3, 4, 1};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {14, 15, -3, 13};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {16, -19, 2, 13};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {17, 20, 1, 19};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {17, 18, -14, 16};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {20, -4, -15, -18};
//+
Plane Surface(8) = {8};
//+
Surface Loop(1) = {5, 7, 6, 8, 3, 4};
//+
Volume(1) = {1};
