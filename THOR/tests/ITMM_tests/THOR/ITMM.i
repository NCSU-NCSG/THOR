ITMM_testing

  start problem
    itmm = yes
    execution = yes
    lambda = 0
    type = keig
    sweep = precomp
    page_refl = save
    innerconv = 1E-6 ; outerconv = 1E-4
    maxinner = 100 ; maxouter = 100
  end problem

  start inout
    mesh_file     =  ../create_mesh/from_CUBIT/ITMM2cell.msh
    xs_file       =  ../cross_sections/ITMM.xs
    flux_file     =  pi_hard_eig_1grp.flux
    print_xs = yes
    source_file = ../src/ITMM.src
  end inout

  start cross_sections
     ngroups = 2 ; upscattering = no
  end cross_sections

  start quadrature
    qdtype = levelsym
    qdorder = 4
  end quadrature

  start regionmap
     1
  end regionmap

end file
