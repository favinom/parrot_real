[Mesh]
 type = GeneratedMesh
 dim = 3
 nx = 15
 ny = 15
 nz = 15
 xmin = 0
 xmax = 1
 ymin = 0
 ymax = 1
 zmin = 0
 zmax = 1
#If you want to use quadratic elements, regenerate the mesh with the following line uncommented.
#elem_type = QUAD9
partitioner = parmetis
[]

[MeshModifiers]
  [./subdomains]
    type = SubdomainBoundingBox
    bottom_left = '0.5 0.5 0.5'
    block_id = 2
    top_right = '1.0 1.0 1.0'
    location = INSIDE
  [../]
[]

[Variables]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x]
 type = Reaction
 variable = pressure
[../]
[]

[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
[]

 
 
[Problem]
 type = FEProblem
 solve = false
[]

[Outputs]
 file_base = mesh0_reg
 print_linear_residuals = true
 print_perf_log = true
 xda = true
 exodus=true
[]

