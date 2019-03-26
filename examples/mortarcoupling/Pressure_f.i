[Mesh]
type = FileMesh
file = mesh_f_tri.e
second_order=true
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
order = SECOND
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

[MultiApps]
 active='transfer_1'
[./transfer_1]
 type =  TransientMultiApp
 app_type=Parrot_realApp
 execute_on=timestep_end
 input_files=Pressure_m.i
 positions='0.0 0.0 0.0'
 [../]
 []


[UserObjects]
 
[./operator]
 type = StoreTransferOperators
 [../]
 
[./u2_slave]
 type=FractureNetworkUserObject
 multi_app = transfer_1
 fracture_variable = p
 matrix_variable = P
 operator_type=MONOLITHIC
 execute_on='timestep_begin'
 solve_cg=true
 solve_mg=false
 pressure=true
 operator_userobject = operator
 transport=false
 constraint_m = false
 constraint_f = false
 porosity_m=1
 porosity_f=1
 id_slave=1
 boundary=true
 [../]
 
 []
