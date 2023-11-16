// Gmsh project created on Wed Feb 15 10:52:30 2023
SetFactory("OpenCASCADE");
//+
y_rem=101.5662;
x_rem=26;
z_rem=40;
//+
Sphere(1) = {30.48, 5, 42.88155, 3.7938, -Pi/2, Pi/2, 2*Pi};
//+
Cylinder(2) = {30.48, 118.8814-y_rem, 28.5877, 0, 0, 33.5153-28.5877, 10.16, 2*Pi};
//+
Box(3) = {30.48-11, 114.1062-y_rem, 28.5877, 22, -10, 6};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{3}; Delete; }
//+
Cylinder(3) = {30.48, 118.8814-y_rem, 33.5153, 0, 0, 35.5753-33.5153, 1.1938, 2*Pi};
//+
Cylinder(4) = {30.48, 118.8814-y_rem, 35.5753, 0, 0, 45.7353-35.5753, 1.1938, 2*Pi};
//+
Cylinder(5) = {30.48, 118.8814-y_rem, 45.7353, 0, 0, 52.2478-45.7353, 1.1938, 2*Pi};
//+
Cylinder(6) = {30.48, 118.8814-y_rem, 33.5153, 0, 0, 52.2478-33.5153, 3.75666, 2*Pi};
//+
Cylinder(7) = {30.48, 118.8814-y_rem, 33.5153, 0, 0, 52.2478-33.5153, 3.91414, 2*Pi};
//+
Cylinder(8) = {30.48, 118.8814-y_rem, 33.5153, 0, 0, 52.2478-33.5153, 10.16, 2*Pi};
//+
Box(9) = {30.48-11, 118.8814-y_rem, 33.5153, 22, -15, 20};
//+
BooleanDifference{ Volume{8}; Delete; }{ Volume{9}; Delete; }
//+
Cylinder(9) = {30.48, 118.8814-y_rem, 33.5153, 0, 0, 52.2478-33.5153, 10.16, 2*Pi};
//+
Box(10) = {30.48-11, 118.8814-y_rem, 33.5153, 22, 15, 20};
//+
Box(11) = {30.48-11, 114.1062-y_rem, 33.5153, 22, -15, 20};
//+
Box(12) = {34.39414, 118.8814-11-y_rem, 33.5153, -22, 22, 20};
//+
BooleanDifference{ Volume{9}; Delete; }{ Volume{10}; Volume{11}; Volume{12}; Delete; }
//+
Cylinder(10) = {30.48, 118.8814-y_rem, 33.5153, 0, 0, 52.2478-33.5153, 10.16, 2*Pi};
//+
Box(11) = {30.48-11, 118.8814-y_rem, 33.5153, 22, 15, 20};
//+
Box(12) = {30.48-11, 114.1062-y_rem, 33.5153, 22, -15, 20};
//+
Box(13) = {26.56586, 118.8814-11-y_rem, 33.5153, 22, 22, 20};
//+
BooleanDifference{ Volume{10}; Delete; }{ Volume{11}; Volume{12}; Volume{13}; Delete; }
//+
Box(11) = {23.1648, 111.5662-y_rem, 28.5877, 37.7952-23.1648, 114.1062-111.5662, 57.1754-28.5877};
//+
Cylinder(12) = {30.48, 118.8814-y_rem, 52.2478, 0, 0, 57.1754-52.2478, 10.16, 2*Pi};
//+
Box(13) = {30.48-11, 114.1062-y_rem, 52.2478, 22, -10, 6};
//+
BooleanDifference{ Volume{12}; Delete; }{ Volume{13}; Delete; }
//+
Box(13) = {0+x_rem/2, 0, 0+z_rem/2, 60.96-x_rem, 135-y_rem, 85.7631-z_rem};
//+
Coherence;
//+
Field[1] = Constant;
//+
Field[1].SurfacesList = {1};
//+
Field[1].VIn = 3;
//+
Field[2] = Constant;
//+
Field[2].VIn = 1;
//+
Field[2].VolumesList = {21, 22};
//+
Field[3] = Constant;
//+
Field[3].VIn = 1;
//+
Field[3].SurfacesList = {2, 6, 9, 12, 14, 17};
//+
Field[4] = Min;
//+
Field[4].FieldsList = {1, 2, 3};
//+
Background Field = 4;
