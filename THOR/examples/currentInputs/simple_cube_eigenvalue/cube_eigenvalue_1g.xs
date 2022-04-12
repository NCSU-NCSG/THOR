!this is current THOR XS format.
!Lines starting with ! will be ignored, comments starting with ! will be ignored
!blank lines will also be ignored
!on the same line directly following the format indicator will be the number of materials, then the
! number of energy groups, and finally the spherical harmonics scattering order

THOR_XS_V1 1 1 0
!each xs can have an integer id index and a name (which cannot contain spaces)
!if no name is given then the name mat_<id> will be assigned
!names, ids, and indicators must be separated by spaces, ids do not need to be in order
!and do not need to be a complete integer sequence
!i.e. a file can contain ids 3, 5, and 6 without ids 1, 2, and 4
id 0 name cube_mat

!1st row of data is fission spectrum chi
1.0
!2nd row of data is energy bounds, typically this can be all 0s
0.0
!3rd row of data is fission cross section SigmaF
0.39
!4th row of data is fission production nu
1.0
!5th row of data is total or transport cross section
1.0
!6th row onward is scattering matrix. Format is traditional with downscatter being literally
! downward, i.e. row 1 is group 1->1, 2->1, 3->1, ...
! row 2 is group 1->2, 2->2, 3->2, ... etc.
.7
