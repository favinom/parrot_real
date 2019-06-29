[Mesh]
type = FileMesh
file = mesh_matrix_3.e
uniform_refine=3
[]

[Variables]
[./CM]
[../]
[]

[AuxVariables]
[./P_aux]
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

[Materials]
[./conductivity1]
type =  HydraulicConductivity
conductivity = 1.0
pressure = P_aux
[../]
[]

[Kernels]
[./convection]
type = Advection
variable = CM
p = P_vel
int_by_parts=false
[../]
[]


[Executioner]
type = Transient
solve_type = LINEAR
dt = 1e-3
num_steps=500.0
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
boundary = '2'
variable = CM
value='1'
[../]
[]

[AuxKernels]
[./flux_x]
type = DiffusionFluxAux
diffusivity = 1.0
variable = flux_x
diffusion_variable = P_aux
component = x
[../]
[./flux_y]
type = DiffusionFluxAux
diffusivity = 1.0
variable = flux_y
diffusion_variable = P_aux
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
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[]

[MultiApps]
#active=''
[./transfer_1]
type =  TransientMultiApp
app_type=Parrot_realApp
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
type=FractureAppNConforming
multi_app = transfer_1
fracture_variable = cm
matrix_variable = CM
lagrange_variable = lambda
operator_type=MONOLITHIC
operator_userobject = operator
execute_on = 'timestep_begin'
solve_cg = true
pressure = false
transport = true
constraint_m = true
constraint_f = false
block_id='101 1'
#  2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18'
value_p='1.0 1.0'
#1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4'
boundary=false
solve_mg=false
stabilize=true
biorth=true
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


