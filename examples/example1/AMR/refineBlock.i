[Problem]
type = FEProblem
solve = false
[]

[Mesh]
file = Mesh_level1.e
block_id = '2 4 5 6 7'
boundary_id = '1 2'
boundary_name = 'inflow outflow'
#second_order=true
[]

[MeshModifiers]
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
 file_base = refineMesh
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
 blockNum = 2
 [../]
 [../]
[]
