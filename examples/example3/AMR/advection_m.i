[Problem]
type = ParrotProblem
[]

[Mesh]
file = test_3_blocchi.e
uniform_refine = 0
#second_order=true
[]

[AuxVariables]
[./P_aux] [../]
[./vel_x] order=CONSTANT  family=MONOMIAL [../]
[./vel_y] order=CONSTANT  family=MONOMIAL [../]
[./vel_z] order=CONSTANT  family=MONOMIAL [../]
[]

[AuxKernels]
[./en]
type = SolutionAux
solution = soln
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
 
[./velxout] type = MaterialRealVectorValueAux i = 0 variable = vel_x property = VelocityVector [../]
[./velyout] type = MaterialRealVectorValueAux i = 1 variable = vel_y property = VelocityVector [../]
[./velzout] type = MaterialRealVectorValueAux i = 2 variable = vel_z property = VelocityVector [../]
 
[]

[Variables]
[./CM] [../]
[]

[Materials]
[./conductivity1]
block='1 2 3 4 5 6 7 8'
type =  HydraulicConductivity
conductivity = 1.e4
 pressure = P_aux
[../]
 
[./conductivity2]
block ='11 12 13'
type =  HydraulicConductivity
conductivity = 1.0
 pressure = P_aux
[../]

[./porosity2] type = Porosity block = '1 2 3 4 5 6 7 8 11 12 13' phi = 0.02  [../]

#[./epsInt2] type = GenericConstantMaterial block = 2 prop_names = epsInt prop_values = 0.01  [../]
#[./epsInt4] type = GenericConstantMaterial block = 4 prop_names = epsInt prop_values = 0.0   [../]
#[./epsInt5] type = GenericConstantMaterial block = 5 prop_names = epsInt prop_values = 0.0   [../]
#[./epsInt6] type = GenericConstantMaterial block = 6 prop_names = epsInt prop_values = 0.0   [../]
#[./epsInt7] type = GenericConstantMaterial block = 7 prop_names = epsInt prop_values = 0.0   [../]

#[./dummy2] type = GenericConstantMaterial block = 2 prop_names = dummy prop_values = 0.0  [../]
#[./dummy4] type = GenericConstantMaterial block = 4 prop_names = dummy prop_values = 0.0  [../]
#[./dummy5] type = GenericConstantMaterial block = 5 prop_names = dummy prop_values = 0.0  [../]
#[./dummy6] type = GenericConstantMaterial block = 6 prop_names = dummy prop_values = 0.0  [../]
#[./dummy7] type = GenericConstantMaterial block = 7 prop_names = dummy prop_values = 1.0 [../]

[]


[Kernels]
active='convection stab time'
#active='time diff'
 
[upwind]
type = AdvectionUpwind
upwinding_type=full
variable = CM
[../]


[./convection]
type = Advection
variable = CM
int_by_parts=false
[../]

[./stab]
type = AdvectionSUPG
variable = CM
coef=0.3
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
type = AnisotropicDiffusion
variable = CM
tensor_coeff='1.e-8 0 0 0 1.0e-8 0 0 0 1.0e-8'
[../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = 21 variable = CM value= 1 [../]
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

dt = 0.01
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
mesh = DiffusionOutput_0.e
timestep = LATEST
system_variables = pressure
execute_on = 'initial'
[../]
[]

 [Postprocessors]
 [./integral1]
 type = ElementIntegralVariablePostprocessor
 variable = CM
 block = 1
 [../]
 
 [./integral2]
 type = ElementIntegralVariablePostprocessor
 variable = CM
 block = 2
 [../]

 [./integral3]
 type = ElementIntegralVariablePostprocessor
 variable = CM
 block = 3
 [../]
 
 [./integral4]
 type = ElementIntegralVariablePostprocessor
 variable = CM
 block = 4
 [../]

 [./integral5]
 type = ElementIntegralVariablePostprocessor
 variable = CM
 block = 5
 [../]
 
 [./integral6]
 type = ElementIntegralVariablePostprocessor
 variable = CM
 block = 6
 [../]
 
 [./integral7]
 type = ElementIntegralVariablePostprocessor
 variable = CM
 block = 7
 [../]
 
 [./integral8]
 type = ElementIntegralVariablePostprocessor
 variable = CM
 block = 8
 [../]
 
 [./volume1] type = VolumePostprocessor block = 1 [../] # execute_on = 'final'
 [./volume2] type = VolumePostprocessor block = 2 [../] # execute_on = 'final'
 [./volume3] type = VolumePostprocessor block = 3 [../] # execute_on = 'final
 [./volume4] type = VolumePostprocessor block = 4 [../] # execute_on = 'final'
 [./volume5] type = VolumePostprocessor block = 5 [../] # execute_on = 'final'
 [./volume6] type = VolumePostprocessor block = 6 [../] # execute_on = 'final'
 [./volume7] type = VolumePostprocessor block = 7 [../] # execute_on = 'final'
 [./volume8] type = VolumePostprocessor block = 8 [../] # execute_on = 'final'

 []
