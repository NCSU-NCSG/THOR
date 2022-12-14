Godiva 6 group example

  start problem
    execution = yes
    lambda = 0
    type = keig ; keigsolver = pi ; piacc = none
    sweep = precomp
    page_refl = save
    kconv = 1E-8; innerconv = 1E-12 ; outerconv = 1E-7
    maxinner = 4 ; maxouter = 5000
  end problem

  start inout
    mesh_file     =  ../create_mesh/from_CUBIT/godiva_1_o.thrm
    xs_file       =  ../cross_sections/godiva.xs
    flux_file     =  pi_hard_eig_1grp.flux
    print_xs = yes
  end inout

  start cross_sections
     ngroups = 6 ; upscattering = no
  end cross_sections

  start quadrature
    qdtype = levelsym
    qdorder = 4
  end quadrature

  start regionmap
     1
  end regionmap

end file
