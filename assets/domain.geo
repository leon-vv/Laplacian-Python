//+
SetFactory("OpenCASCADE");
Disk(1) = {0, 0, 0, 2, 2};
//+
Rectangle(2) = {-0.5, -0.5, 0, 1, 1, 0};
//+
Curve Loop(3) = {1};
//+
Curve Loop(4) = {4, 5, 2, 3};
//+
Plane Surface(3) = {3, 4};
//+
Physical Curve("outer-border") = {1};
//+
Physical Curve("inner-border") = {5, 4, 3, 2};
//+
Physical Surface("domain") = {3};
