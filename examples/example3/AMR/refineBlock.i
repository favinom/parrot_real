[Problem]
type = FEProblem
solve = false
[]

[Mesh]
file =  test_3_blocchi_c2.e
#second_order=true
[]

[MeshModifiers]
active=''
[./rotate]
type = Transform
transform = TRANSLATE
vector_value = '50 50 50'
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

#[BCs]
#[./Periodic]
#[./pressure_real_periodic] variable = pressure auto_direction = 'x y z' [../]
#[../]
#[]


[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
 nl_abs_tol = 1e-8
 []

[Outputs]
 file_base = refineMesh_1
 exodus = true
 print_linear_residuals = true
 print_perf_log = true
 xdr = true
[]

[Adaptivity]
 marker = simplemark
 steps = 2
 [./Markers]
 [./simplemark]
 type = BlockMarker
 blockNum = '1 2 3 4 5 6 7 8'
 [../]
 
 [../]
[]
