[Problem]
type = ParrotProblem
[]

[Mesh]
 file = refineMesh_0003_mesh.xdr
 block_id = '2 4 5 6 7'
 boundary_id = '1 2'
 boundary_name = 'inflow outflow'
 uniform_refine = 0
#second_order = true
 []

[AuxVariables]
[./P_aux] [../]
[]


[MeshModifiers]
active=''
[./rotate]
type = Transform
transform = TRANSLATE
vector_value = '50 50 50'
[../]
[]


[AuxKernels]
[./en]
type = SolutionAux
solution = soln
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[]

[Variables]
[./CM] [../]
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
active='diff Stab_f convection stab  time'

[upwind]
type = AlgebraicDiffusion
#upwinding_type=full
variable = CM
#int_by_parts=false
[../]

[./Stab_f]
block='2'
type = MyDiffusion
variable = CM
coef=5.e-3
[../]

[./convection]
type = Advection
variable = CM
int_by_parts=false
[../]

[./stab]
type = AdvectionSUPG
variable = CM
coef=0.17  #0.13
use_h=true
[../]

[./timestab]
type = TimeAdvectionSUPG
variable = CM
coef=0.3
use_h=true
lumping=true
[../]

[./time]
type = PorosityTimeDerivative
variable = CM
lumping = true
[../]

[./diff]
block='4 5'
type = MyDiffusion
variable = CM
coef=1.0e-3
#tensor_coeff='1.e-10 0 0 0 1.0e-10 0 0 0 1.0e-10'
[../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = inflow variable = CM value='0.01' [../]
#[./u_injection_right] type = OutflowBC boundary = 'outflow 3' variable = CM [../]
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

[Executioner]

type = Transient
solve_type= LINEAR
line_search = none

petsc_options_iname=' -ksp_type            '
petsc_options_value='  ksp_parrot_preonly  '

dt = 1e7
num_steps=100

[./Quadrature]
order=SIXTH
[../]

[]

[Outputs]
exodus = true
csv=true
print_perf_log = true
[]


[UserObjects]
[./soln]
type = SolutionUserObject
mesh = OutputBenchmark1.e
timestep = 2
system_variables = pressure
execute_on = 'initial'
[../]
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





