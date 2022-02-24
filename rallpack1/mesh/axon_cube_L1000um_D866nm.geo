Point(1) = {4.431135e-07, 4.431135e-07, 0, 5e-4};
Point(2) = {4.431135e-07, -4.431135e-07, 0, 5e-4};
Point(3) = {-4.431135e-07, 4.431135e-07, 0, 5e-4};
Point(4) = {-4.431135e-07, -4.431135e-07, 0, 5e-4};
Point(5) = {4.431135e-07, 4.431135e-07, 1e-3, 5e-4};
Point(6) = {4.431135e-07, -4.431135e-07, 1e-3, 5e-4};
Point(7) = {-4.431135e-07, 4.431135e-07, 1e-3, 5e-4};
Point(8) = {-4.431135e-07, -4.431135e-07, 1e-3, 5e-4};

Line(1) = {8, 6};
Line(2) = {6, 5};
Line(3) = {5, 7};
Line(4) = {7, 8};
Line(5) = {4, 2};
Line(6) = {2, 1};
Line(7) = {1, 3};
Line(8) = {3, 4};
Line(9) = {6, 2};
Line(10) = {5, 1};
Line(11) = {7, 3};
Line(12) = {8, 4};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {1, 9, -5, -12};
Plane Surface(3) = {3};
Curve Loop(4) = {2, 10, -6, -9};
Plane Surface(4) = {4};
Curve Loop(5) = {3, 11, -7, -10};
Plane Surface(5) = {5};
Curve Loop(6) = {12, -8, -11, 4};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 3, 4, 5, 6, 2};
Volume(1) = {1};

Physical Surface("z_max", 1) = {1};
Physical Surface("z_min", 2) = {2};
Physical Surface("memb", 3) = {6, 3, 4, 5};

Physical Volume("cyto", 4) = {1};
