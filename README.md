# THOR

THOR is a discrete ordinates radiation transport code using the AHOT-C method on unstructured meshes (can be generated using [Gmsh](https://gmsh.info/) and converted using [OpenMeshConverter](https://github.com/nfherrin/OpenMeshConverter)) and multigroup XS (can be generated using [OpenMC](https://github.com/openmc-dev/openmc) and converted using [OpenXSConverter](https://github.com/nfherrin/OpenXSConverter)).

For usage, see the user's manual in `<THOR_dir>/docs/usermanual/CurrentVersion/usermanual.pdf`.
This can be found online at https://github.com/NCSU-NCSG/THOR/raw/v1.1.0/docs/usermanual/CurrentVersion/usermanual.pdf .
The theory manual can be found in `<THOR_dir>/docs/theorymanual/CurrentVersion/theorymanual.pdf`.
The programmer's manual is currently just an outline and is not yet considered released.

NOTE: The current release version of THOR (the default branch in this repository) is not necessarily the most up to date version of THOR.
The most up to date tested version of THOR can be found in the `master` branch which can be checked out after cloning the repository.
All features merged into `master` have been tested via the regression tests, however they may not be complete if they are part of a larger over-arching feature OR they may not yet have adequate documentation in the THOR manuals.
