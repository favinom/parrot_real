
[Mesh]
type = FileMesh
file = refineMesh_0_0001_mesh.xdr
[]

[Variables]
[./CM]
[../]
[]

[MeshModifiers]
#active=''
[./rotate]
type = Transform
transform = TRANSLATE
vector_value = '50 50 50'
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
[./conductivity2] type = HydraulicConductivityFracture block = 2 conductivity = '1e-1 1e-1 1e-1' theta = '0.0 30.9638 0.0' pressure = P_aux [../]
[./conductivity4] type = HydraulicConductivity         block = 4 conductivity = 1e-6                                       pressure = P_aux [../]
[./conductivity5] type = HydraulicConductivity         block = 5 conductivity = 1e-6                                       pressure = P_aux [../]
[./conductivity6] type = HydraulicConductivity         block = 6 conductivity = 1e-6                                       pressure = P_aux [../]
[./conductivity7] type = HydraulicConductivity         block = 7 conductivity = 1e-5                                       pressure = P_aux [../]

[./porosity2] type = Porosity block = 2 phi = 0.4  [../]
[./porosity4] type = Porosity block = 4 phi = 0.2  [../]
[./porosity5] type = Porosity block = 5 phi = 0.2  [../]
[./porosity6] type = Porosity block = 6 phi = 0.2  [../]
[./porosity7] type = Porosity block = 7 phi = 0.25 [../]

[./epsInt2] type = GenericConstantMaterial block = 2 prop_names = epsInt prop_values = 0.01  [../]
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
boundary = '1'
variable = CM
value='0.01'
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
mesh = OutputBenchmark1.e
timestep = 2
system_variables = pressure
execute_on = 'timestep_begin'
[../]

[./operator]
type = StoreTransferOperators
[../]

[./u2_slave]
type=FractureAppConforming
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
csv = true
execute_on = 'timestep_end'
[]


[Postprocessors]

[./flux_left]
type = SideFlux
variable = CM
boundary = outflow
# PER FAVORE CONTROLLARE IL COEF
coef = 0.0
execute_on = 'timestep_end'
[../]

[./int3]
type = ElementIntegral_phi_c_MatProp
variable = CM
mat_prop = dummy
execute_on = 'timestep_end'
[../]

[./intFrac]
type = ElementIntegral_phi_c_MatProp
variable = CM
mat_prop = epsInt
execute_on = 'timestep_end'
[../]

[]
