[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin =-0.5
  xmax = 0.5
  ymin =-0.5
  ymax = 0.5
  nx = 40
  ny = 40
  elem_type = QUAD4
[]

#[Mesh]
#type = FileMesh
#file = mesh_m_q3.e
#second_order=false
#[]


[Variables]
[./P]
order = FIRST
[../]
[]





[Kernels]
[./diff_1]
type = Diffusion
variable = P
[../]
[]

[BCs]
# active=''
[./_b_22]
 type = DirichletBC
 variable = P
 boundary = right
 value=4
[../]
 
[./_b_21]
 type = DirichletBC
 variable = P
 boundary = left
 value=1
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



[MultiApps]
active='transfer_1'
[./transfer_1]
type =  TransientMultiApp
app_type=Parrot_realApp
execute_on=timestep_end
input_files=Pressure_f.i
positions='0.0 0.0 0.0'
[../]
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
