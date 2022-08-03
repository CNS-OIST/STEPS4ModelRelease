SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthFactor = 100;

Point(1) = {-4, -4, 0, 100};
Point(2) = {4, -4, 0, 100};
Point(3) = {4, 4, 0, 100};
Point(4) = {-4, 4, 0, 100};

Point(5) = {-4, -4, 10, 100};
Point(6) = {4, -4, 10, 100};
Point(7) = {4, 4, 10, 100};
Point(8) = {-4, 4, 10, 100};

Point(9) = {-4, -4, 20, 100};
Point(10) = {4, -4, 20, 100};
Point(11) = {4, 4, 20, 100};
Point(12) = {-4, 4, 20, 100};

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
Line(13) = {11, 10};
Line(14) = {10, 9};
Line(15) = {9, 12};
Line(16) = {12, 11};
Line(17) = {11, 7};
Line(18) = {10, 6};
Line(19) = {9, 5};
Line(20) = {12, 8};

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
Curve Loop(7) = {19, -8, -20, -15};
Plane Surface(7) = {7};
Curve Loop(8) = {15, 16, 13, 14};
Plane Surface(8) = {8};
Curve Loop(9) = {14, 19, 5, -18};
Plane Surface(9) = {9};
Curve Loop(10) = {18, 6, -17, 13};
Plane Surface(10) = {10};
Curve Loop(11) = {17, 7, -20, 16};
Plane Surface(11) = {11};

Surface Loop(1) = {2, 6, 5, 1, 4, 3};
Volume(1) = {1};
Surface Loop(2) = {8, 7, 9, 10, 11, 2};
Volume(2) = {2};

Physical Volume("spiny", 1) = {1};
Physical Volume("smooth", 2) = {2};

