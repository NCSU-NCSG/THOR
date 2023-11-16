// Gmsh project created on Wed Oct 26 17:08:27 2022
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.39218, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 0.40005, 0, 2*Pi};
//+
Circle(3) = {0, 0, 0, 0.45720, 0, 2*Pi};
//+
Point(4) = {-0.62992, -0.62992, 0, 1.0};
//+
Point(5) = {0.62992, -0.62992, 0, 1.0};
//+
Point(6) = {0.62992, 0.62992, 0, 1.0};
//+
Point(7) = {-0.62992, 0.62992, 0, 1.0};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 4};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2};
//+
Curve Loop(3) = {1};
//+
Plane Surface(2) = {2, 3};
//+
Curve Loop(4) = {3};
//+
Curve Loop(5) = {2};
//+
Plane Surface(3) = {4, 5};
//+
Curve Loop(6) = {5, 6, 7, 4};
//+
Curve Loop(7) = {3};
//+
Plane Surface(4) = {6, 7};
//+
h=100;
runsize=2;
//+
Extrude {0, 0, h/2} {
  Surface{1};
  Layers{ h/(runsize*12) };
}
//+
Extrude {0, 0, h/2} {
  Surface{2};
  Layers{ h/(runsize*12) };
}
//+
Extrude {0, 0, h/2} {
  Surface{3};
  Layers{ h/(runsize*12) };
}
//+
Extrude {0, 0, h/2} {
  Surface{4};
  Layers{ h/(runsize*12) };
}
//+
tn=10/runsize;
Transfinite Curve {1} = tn Using Progression 1;
//+
Transfinite Curve {2} = tn Using Progression 1;
//+
Transfinite Curve {3} = tn Using Progression 1;
//+
Transfinite Curve {6} = tn/3 Using Progression 1;
//+
Transfinite Curve {5} = tn/3 Using Progression 1;
//+
Transfinite Curve {7} = tn/3 Using Progression 1;
//+
Transfinite Curve {4} = tn/3 Using Progression 1;
