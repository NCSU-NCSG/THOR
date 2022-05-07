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
Physical Volume(1) = {1};
// The shell has outer radius of 11.4138 cm
Sphere(2) = {0,0,0,11.4138};
Physical Volume(2) = {2};
// now we want a consistent mesh using fragments
v() = BooleanFragments{ Volume{2}; Delete; }{ Volume{1}; Delete; };

// restores tag to outer sphere
Physical Volume(4) = v(#v()-1);

// Assign a mesh size to all the points of all the volumes:
MeshSize{ PointsOf{ Volume{:}; } } = 6;
// should be smaller for berp
MeshSize{ PointsOf{ Volume{1}; } } = 2;