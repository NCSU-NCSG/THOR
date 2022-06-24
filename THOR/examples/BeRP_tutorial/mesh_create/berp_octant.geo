// -----------------------------------------------------------------------------
//
//  BeRP in Poly Full Sphere
//
//  Constructive Solid Geometry, OpenCASCADE geometry kernel
//
// -----------------------------------------------------------------------------

// Instead of constructing a model in a bottom-up fashion with Gmsh's built-in
// geometry kernel, starting with version 3 Gmsh allows you to directly use
// alternative geometry kernels. Here we use the OpenCASCADE kernel:

SetFactory("OpenCASCADE");

// BeRP is just a sphere of radius 3.7938 cm
Sphere(6) = {0,0,0,3.7938};
// The shell has outer radius of 11.4138 cm
Sphere(7) = {0,0,0,11.4138};

// boxes to remove
Box(3) = {-12,-12,-12, 12,24,24};
Box(4) = {0,-12,-12, 12,24,12};
Box(5) = {0,-12,0, 12,12,12};

// remove bottom half and negative side
BooleanDifference(1) = { Volume{6}; Delete; }{ Volume{3:5}; };
Physical Volume(1) = {1};
BooleanDifference(2) = { Volume{7}; Delete; }{ Volume{3:5}; Delete; };
Physical Volume(2) = {2};

// now we want a consistent mesh using fragments
v() = BooleanFragments{ Volume{2}; Delete; }{ Volume{1}; Delete; };

// restores tag to outer sphere
Physical Volume(9) = v(#v()-1);

// Assign a mesh size to all the points of all the volumes:
MeshSize{ PointsOf{ Volume{:}; } } = 6;
// should be smaller for berp
MeshSize{ PointsOf{ Volume{1}; } } = 2;
