
[Mesh]
type = GeneratedMesh
x_min=0
x_max=1.0
y_min=0
y_max=1
z_min=0
z_max=1
nx = 20
ny = 20
nz = 20
dim=3
[]

[Variables]
[./CM]
[../]
[]

[AuxVariables]
[./Pvel]
[../]
[]

[ICs]
[./u_ic]
type = FunctionIC
variable = 'Pvel'
function = '5*x-5*z'
[../]
[]



[Kernels]
[./convection]
type = Advection
variable = CM
p=Pvel
#velocity='-1 0 0'
int_by_parts=true
[]

[./advectionsupg]
type = AdvectionSUPG
coef=0.5
use_h=true
variable = CM
p=Pvel
[../]


[./time]
type = TimeDerivative
variable = CM
[../]
[]


[Executioner]
type = Transient
solve_type = LINEAR
dt = 0.1
num_steps=10.0
l_tol = 1E-14
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

[BCs]
active='u_injection_left'
[./u_injection_left]
 type = DirichletBC
 boundary = 'right'
 variable = CM
 value='1'
[../]

[./u_injection_right]
 type = OutflowBC
 boundary = 'left top bottom front back'
 variable = CM
 p=Pvel
[../]
[]


[Problem]
type = FEProblem
solve = true
kernel_coverage_check = false
[]


[Functions]
[./ic_func]
type = ParsedFunction
value = '10*(x<0.9)*(x>0.7)*(y>0.7)*(y<0.9)'
[../]
[]

[Outputs]
file_base = cons
exodus = true
execute_on = 'timestep_end'
[]


