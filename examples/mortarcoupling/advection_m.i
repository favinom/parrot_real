
[Mesh]
type = FileMesh
file = test_1.e
[]

[Variables]
[./CM]
[../]
[]

[AuxVariables]
[./P_vel]
[../]
[./flux_x]
order = FIRST
family = MONOMIAL
[../]
[./flux_y]
order = FIRST
family = MONOMIAL
[../]
[]

[Kernels]
[./convection]
type = Advection_Pressure
variable = CM
coef=0.1
p = P_vel
int_by_parts=false
[../]

#[./time]
#type = TimeDerivative
#variable = CM
#[../]
[]


[Executioner]
type = Transient
solve_type = LINEAR
dt = 0.1
num_steps=100.0
l_tol = 1E-14
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]



[BCs]
[./u_injection_left]
 type = DirichletBC
 boundary = '1'
 variable = CM
 value='1'
[../]
[]

[AuxKernels]
[./flux_x]
type = DiffusionFluxAux
diffusivity = 1.0
variable = flux_x
diffusion_variable = P_vel
component = x
[../]
[./flux_y]
type = DiffusionFluxAux
diffusivity = 1.0
variable = flux_y
diffusion_variable = P_vel
component = y
[../]
[]


[Problem]
type = FEProblem
solve = false
kernel_coverage_check = false
[]


[AuxKernels]
#active=''
[./en]
type = SolutionAux
solution = soln
variable = P_vel
scale_factor = 1.0
execute_on = 'initial'
[../]
[]

[MultiApps]
#active=''
[./transfer_1]
type =  TransientMultiApp
app_type=WhaleApp
execute_on='timestep_end'
input_files=advection_f.i
positions='0.0 0.0 0.0'
[../]
[]


[UserObjects]
#active='soln'
[./soln]
type = SolutionUserObject
mesh = matrix.e
timestep = 1
system_variables = P
execute_on = 'timestep_begin'
[../]

[./operator]
type = StoreTransferOperators
[../]

[./u2_slave]
type=FractureNetworkUserObject
multi_app = transfer_1
fracture_variable = cm
matrix_variable = CM
lagrange_variable = lambda
operator_type=MONOLITHIC
operator_userobject = operator
execute_on = 'timestep_begin'
solve_cg = false
pressure = false
transport = true
constraint_m = true
constraint_f = false
porosity_m=0.2
porosity_f=2.0e-5
[../]

[]

[Functions]
[./ic_func]
type = ParsedFunction
value = '10*(x<0.9)*(x>0.7)*(y>0.7)*(y<0.9)'
[../]
[]

[Outputs]
exodus = true
execute_on = 'timestep_end'
[]


