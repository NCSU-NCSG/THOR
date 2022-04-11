[Mesh]
  type = FileMesh
  #file = c5g7-2d_r1.e
[]

[MeshModifiers]
  [./extrude]
    type = MeshExtruder
    num_layers = 16
    extrusion_vector = '0 0 57.12'
    existing_subdomains = '1 2 3 4 5 6 7'
    layers = '6 7 8 9'
    new_ids = '1 1 8 4 1 1 1
               1 1 8 4 1 1 1
               1 1 8 4 1 1 1
               1 1 8 4 1 1 1'
    bottom_sideset = '2'
    top_sideset = '4'
  [../]
[]

[UserObjects]
  [./dump_mesh]
    type = DumpMesh
    #filename = 'c5g7_3d_hex.dat'
  [../]
[]

[Variables]
  active = 'u'

  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'diff'

  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [./bottom]
    type = DirichletBC
    variable = u
    boundary = 2
    value = 0
  [../]

  [./top]
    type = DirichletBC
    variable = u
    boundary = 4
    value = 1
  [../]
[]

[Executioner]
  type = Steady

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  nl_max_its = 2
  l_max_its = 10
  nl_rel_tol = 0.1
[]

[Outputs]
  file_base = extrusion_out
  exodus = true
[]
