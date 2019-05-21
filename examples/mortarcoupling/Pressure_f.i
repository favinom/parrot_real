[Mesh]
type = FileMesh
file = mesh_f_tri.e
second_order=false
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
[./pp_slave]
[../]
[]

[Kernels]
[./diff1]
type = Diffusion
variable = p
[../]
[]

[BCs]
#active=''
[./_b_11]
type = DirichletBC
variable = p
boundary = 11
value=4
[../]
[./_b_12]
type = DirichletBC
variable = p
boundary = 12
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


