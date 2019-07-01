[Mesh]
type = FileMesh
file = mesh_fracture_2d_non_conf.e
uniform_refine=1
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
 conductivity = 1e-3
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
 active=''
[./x_disp]
type = DirichletBC
variable = p
boundary = 1
value=1
[../]
[./y_disp]
type = DirichletBC
variable = p
boundary = 2
value=4
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


