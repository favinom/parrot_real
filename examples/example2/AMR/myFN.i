 [Mesh]
file = adapt1.xda
# uniform_refine = 1
#partitioner = parmetis
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

[Adaptivity]
  marker = simplemark
  steps = 2
  [./Markers]
    [./simplemark]
      type = SimpleMarker3D
fn = 9
fx_string = '0.5,0.5,0.5,
0.749975,0.75,0.749975,
0.625,0.625,0.625'
fy_string = '0.5,0.5,0.5,
0.749975,0.749975,0.75,
0.625,0.625,0.625'
fz_string = '0.5,0.5,0.5,
0.75,0.749975,0.749975,
0.625,0.625,0.625'
fa1_string = '0.0,0.0,0.0,
0.0,0.0,0.0,
0.0,0.0,0.0'
fa2_string = '0.0,90.0,0.0,
0.0,90.0,0.0,
0.0,90.0,0.0'
fa3_string = '0.0,0.0,90.0,
0.0,0.0,90.0,
0.0,0.0,90.0'
fd1_string = '1.0,1.0,1.0,
0.50005,0.50005,0.50005,
0.2501,0.2501,0.2501'
fd2_string = '1.0,1.0,1.0,
0.50005,0.50005,0.50005,
0.2501,0.2501,0.2501'
fd3_string = '0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001'
    [../]
  [../]
[]

 
 
 

 
 
[Problem]
 type = FEProblem
 solve = false
[]

[Outputs]
 file_base = adaptedAMR1
 exodus = true
 print_linear_residuals = true
 print_perf_log = true
 xda = true
[]

