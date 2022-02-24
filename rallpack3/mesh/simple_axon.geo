SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthFactor = 10;

Point(1) = {-4.431135e-7, -4.431135e-7, 0, 1.0};
Point(2) = {4.431135e-7, -4.431135e-7, 0, 1.0};
Point(3) = {4.431135e-7, 4.431135e-7, 0, 1.0};
Point(4) = {-4.431135e-7, 4.431135e-7, 0, 1.0};

Point(5) = {-4.431135e-7, -4.431135e-7, 1.0e-6, 1.0};
Point(6) = {4.431135e-7, -4.431135e-7, 1.0e-6, 1.0};
Point(7) = {4.431135e-7, 4.431135e-7, 1.0e-6, 1.0};
Point(8) = {-4.431135e-7, 4.431135e-7, 1.0e-6, 1.0};

Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {2, 6};
Line(10) = {7, 3};
Line(11) = {4, 8};
Line(12) = {5, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {9, 6, 10, 4};
Plane Surface(3) = {3};
Curve Loop(4) = {10, -3, 11, -7};
Plane Surface(4) = {4};
Curve Loop(5) = {12, 2, 11, 8};
Plane Surface(5) = {5};
Curve Loop(6) = {12, -1, 9, -5};
Plane Surface(6) = {6};

Surface Loop(1) = {2, 6, 5, 1, 4, 3};
Volume(1) = {1};

Physical Surface("memb", 1) = {3, 4, 5, 6};
Physical Surface("z_max", 2) = {2};
Physical Surface("z_min", 3) = {1};
Physical Volume("cyto", 4) = {1};
