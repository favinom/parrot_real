[Mesh]
 type = GeneratedMesh
 dim = 3
 nx = 20
 ny = 20
 nz = 20
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
 file_base = uniformOut
 print_linear_residuals = true
 print_perf_log = true
 xda = true
[]

