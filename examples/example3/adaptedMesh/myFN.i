 [Mesh]
type = GeneratedMesh
dim = 3
xmin = 0.0
xmax = 1.0000000001
ymin = 0.0
ymax = 2.25
zmin = 0.0
zmax = 1.0
nx = 9
ny = 19
nz = 9
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
  steps = 5
  [./Markers]
    [./simplemark]
      type = SimpleMarker3D
fn = 8
fx_string = '0.5  ,0.5  ,0.77,0.83, 0.2,0.2  , 0.5,0.5'
fy_string = '1.125,0.175,2.05,2.05, 2.05,2.05 , 1.6,1.6'
fz_string = '0.5  ,0.5  ,0.5 ,0.5, 0.5,0.5 , 0.675,0.31'
fa1_string = '0,90,90,90,78.6901,-78.6901,0,0'
fa2_string = '0, 0, 0, 0,0,0,0,0'
fa3_string = '0,90,90,90,90,90,16.2602,-15.8192'
fd1_string = '0.9,0.25,0.3,0.3,0.3059,0.3059,0.9,0.9'
fd2_string = '1.75,0.9,0.4,0.4,0.4,0.4,1.25,1.2472'
fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
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

