mesh_size = 0.02; // 0.02 is pretty good
//+
Point(1) = {-0.5, 0, 0, 1.0};
//+
Point(2) = {-0.5, -0.5, 0, 1.0};
//+
Point(3) = {0, -0.5, 0, 1.0};
//+
Point(4) = {0, -1.5, 0, 1.0};
//+
Point(5) = {3, -1.5, 0, 1.0};
//+
Point(6) = {3, -0.5, 0, 1.0};
//+
Point(7) = {3.5, -0.5, 0, 1.0};
//+
Point(8) = {3.5, 0, 0, 1.0};
//+
Point(9) = {3.0, 0, 0, 1.0};
//+
Point(10) = {0, 0, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 6};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 1};
//+
Line(11) = {1, 2};
//+
Line(12) = {3, 10};
//+
Curve Loop(1) = {10, 11, 1, 12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, -12, 2, 3, 4, -8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 8, 5, 6};
//+
Plane Surface(3) = {3};
//+
MeshSize{ PointsOf{ Surface{1, 2, 3}; } } = mesh_size;//+
Physical Curve(1) = {11};
//+
Physical Curve(2) = {6};
//+
Physical Curve(3) = {1, 2, 3, 4, 5};
//+
Physical Curve(4) = {10, 7, 9};
//+
Physical Surface(5) = {2};
Physical Surface(6) = {1};
Physical Surface(7) = {3};
