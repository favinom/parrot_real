[Mesh]
type = FileMesh
file = mesh_fracture_3.e
#2D_1_square_frac_in_matrix_multiple_insquares_frac_new.e
construct_side_list_from_node_list = true
uniform_refine=3
[]

[Variables]
[./p]
order = FIRST
family = LAGRANGE
[../]
[]

[AuxVariables]
[./lambda]
[../]
[./p_aux]
[../]
[]

[Materials]
[./conductivity1]
 type =  HydraulicConductivity
 conductivity = 1.0
 pressure = p_aux
 [../]
 []
 
[Kernels]
[./diff_1]
 type = MyDiffusion
 variable = p
 coef=1.0
 [../]
[]

[BCs]
 active='x_disp'
[./x_disp]
type = DirichletBC
variable = p
boundary = '1 3'
value=1
[../]
[./y_disp]
type = NeumannBC
variable = p
boundary = 2
value=1
[../]
[]


[Problem]
type = FEProblem
solve = false
kernel_coverage_check = false
[]



[Executioner]
type = Transient
dt = 1.0
start_time = 0.0
num_steps = 1.0
solve_type ='LINEAR'
[]


