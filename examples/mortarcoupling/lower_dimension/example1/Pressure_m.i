[Mesh]
type = FileMesh
file = mesh_matrix.e
[]

[Variables]
[./P]
order = FIRST
family = LAGRANGE
[../]
[]


[AuxVariables]
[./P_aux]
[../]
[]

[Materials]
[./conductivity1]
 block='1'
 type =  HydraulicConductivity
 conductivity = 1e-6
 pressure = P_aux
[../]
[./conductivity2]
 block='2'
 type =  HydraulicConductivity
 conductivity = 1e-5
 pressure = P_aux
[../]
[]
 
[Kernels]
[./diff_1]
type = MyDiffusion
variable = P
coef=1.0
[../]
[]

[BCs]
[./x_disp]
type = DirichletBC
variable = P
boundary = 1
value=1
[../]
[./y_disp]
type = DirichletBC
variable = P
boundary = 2
value=4
[../]
[]


[MultiApps]
active='transfer_1'
[./transfer_1]
type =  TransientMultiApp
app_type=Parrot_realApp
execute_on=timestep_begin
input_files=Pressure_f.i
positions='0.0 0.0 0.0'
[../]
[]


[Problem]
type = FEProblem
solve = false
kernel_coverage_check = false
[]


[Preconditioning]
[./SMP]
type = SMP
full = true
ksp_norm = default
[../]
[]

[Executioner]
type = Transient
dt = 1.0
start_time = 0.0
end_time = 1.0
solve_type ='LINEAR'
[]

[Outputs]
file_base = matrix
exodus = true
execute_on = 'timestep_end'
[]

[UserObjects]

[./operator]
type = StoreTransferOperators
[../]

[./u2_slave]
type=FractureAppNConforming
multi_app = transfer_1
fracture_variable = p
matrix_variable = P
lagrange_variable = lambda
operator_type=CG_DUAL
execute_on='timestep_begin'
solve_cg=true
pressure=true
operator_userobject = operator
transport=false
constraint_m = false
constraint_f = false
porosity_m=1
porosity_f=1
solve_mg=false
boundary=false
block_id='1 2 3 4'
value_p='0.25 0.2 0.2 0.4'
[../]

[]


