
#[Mesh]
#type = FileMesh
#file = refineMesh_0001_mesh.xdr
#[]

[Mesh]
type = GeneratedMesh
dim = 3
xmin = 0
ymin = 0
zmin = 0
xmax = 100
ymax = 100
zmax = 100
nx = 50
ny = 50
nz = 50
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
active='conductivity2 porosity2'
[./conductivity2] type = HydraulicConductivityFracture block = 0 conductivity = '1e-1 1e-1 1e-1' theta = '0.0 30.9638 0.0' pressure = P_aux [../]
[./conductivity4] type = HydraulicConductivity         block = 4 conductivity = 1e-6                                       pressure = P_aux [../]
[./conductivity5] type = HydraulicConductivity         block = 5 conductivity = 1e-6                                       pressure = P_aux [../]
[./conductivity6] type = HydraulicConductivity         block = 6 conductivity = 1e-6                                       pressure = P_aux [../]
[./conductivity7] type = HydraulicConductivity         block = 7 conductivity = 1e-5                                       pressure = P_aux [../]

[./porosity2] type = Porosity block = 0 phi = 0.4  [../]
[./porosity4] type = Porosity block = 4 phi = 0.2  [../]
[./porosity5] type = Porosity block = 5 phi = 0.2  [../]
[./porosity6] type = Porosity block = 6 phi = 0.2  [../]
[./porosity7] type = Porosity block = 7 phi = 0.25 [../]

[./epsInt2] type = GenericConstantMaterial block = 0 prop_names = epsInt prop_values = 0.01  [../]
[./epsInt4] type = GenericConstantMaterial block = 4 prop_names = epsInt prop_values = 0.0   [../]
[./epsInt5] type = GenericConstantMaterial block = 5 prop_names = epsInt prop_values = 0.0   [../]
[./epsInt6] type = GenericConstantMaterial block = 6 prop_names = epsInt prop_values = 0.0   [../]
[./epsInt7] type = GenericConstantMaterial block = 7 prop_names = epsInt prop_values = 0.0   [../]

[./dummy2] type = GenericConstantMaterial block = 2 prop_names = dummy prop_values = 0.0  [../]
[./dummy4] type = GenericConstantMaterial block = 4 prop_names = dummy prop_values = 0.0  [../]
[./dummy5] type = GenericConstantMaterial block = 5 prop_names = dummy prop_values = 0.0  [../]
[./dummy6] type = GenericConstantMaterial block = 6 prop_names = dummy prop_values = 0.0  [../]
[./dummy7] type = GenericConstantMaterial block = 7 prop_names = dummy prop_values = 1.0 [../]

[]


[Kernels]
[./convection]
type = Advection
variable = CM
int_by_parts=false
[../]
[]


[Executioner]
type = Transient
solve_type = LINEAR
dt = 1.0e7
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
boundary = 'right'
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
active=''
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
mesh = Diffusion_out_sub0.e
timestep = LATEST
system_variables = v
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
solve_mg=false
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


