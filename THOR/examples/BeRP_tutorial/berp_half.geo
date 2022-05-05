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
Sphere(1) = {0,0,0,3.7938};
// The shell has outer radius of 11.4138 cm
Sphere(2) = {0,0,0,11.4138};

// boxes to remove
Box(3) = {-12,-12,-12, 12,24,24};

// remove bottom half
BooleanDifference(4) = { Volume{1}; Delete; }{ Volume{3}; };
Physical Volume(4) = {4};
BooleanDifference(5) = { Volume{2}; Delete; }{ Volume{3}; Delete; };
Physical Volume(5) = {5};

// now we want a consistent mesh using fragments
v() = BooleanFragments{ Volume{5}; Delete; }{ Volume{4}; Delete; };

// restores tag to outer sphere
Physical Volume(7) = v(#v()-1);

// Assign a mesh size to all the points of all the volumes:
MeshSize{ PointsOf{ Volume{:}; } } = 6;
// should be smaller for berp
MeshSize{ PointsOf{ Volume{4}; } } = 2;
