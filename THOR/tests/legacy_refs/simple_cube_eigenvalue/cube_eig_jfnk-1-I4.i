infinite medium 1 group w/ upscattering

  start problem
    execution = yes
    lambda = 0
    type=keig ; keigsolver = jfnk ; jfnk_method = outer
    sweep=precomp
    page_refl=save
    kconv=1E-8; innerconv=1E-12 ; outerconv=1E-7
    maxinner = 4 ; maxouter = 5000
  end problem

  start inout
    mesh_file     =  cube_level1.mesh
    xs_file       =  cube_eigenvalue_1g.xs
    flux_file     =  pi_hard_eig_1grp.flux
    print_xs=yes
  end inout

  start cross_sections
     ngroups=1 ; upscattering = yes
  end cross_sections

  start quadrature
qdtype=levelsym ; qdorder = 12
  end quadrature

  start regionmap
     0
  end regionmap

end file
