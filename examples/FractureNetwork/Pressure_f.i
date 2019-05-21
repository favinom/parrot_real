[Mesh]
type = FileMesh
file = inclined_line.e
[]

#[Mesh]
# type = GeneratedMesh
# dim = 1
# ymin = 0
# ymax = 1.0
# ny = 10
# elem_type = EDGE3
#[]


[Variables]
[./p]
order = FIRST
family = LAGRANGE
[../]
[]

[AuxVariables]
[./lambda]
[../]
[./pp_slave]
[../]
[]

[Kernels]
[./diff1]
type = CoefDiffusion
variable = p
coef=0.1
[../]
[]

[BCs]
 active=''
[./x_disp]
type = DirichletBC
variable = p
boundary = 1
value=4
[../]
[./y_disp]
type = DirichletBC
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


