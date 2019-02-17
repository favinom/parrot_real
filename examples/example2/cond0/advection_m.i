[Problem]
type = ParrotProblem
[]

[Mesh]
  file = ../refineMesh_0_0003_mesh.xda
  #uniform_refine = 1
  #second_order=true
[]

[AuxVariables]
[./P_aux] [../]
[./vel_x]family=MONOMIAL [../]
[./vel_y]family=MONOMIAL [../]
[./vel_z]family=MONOMIAL [../]
[]

[AuxKernels]
[./en]
type = SolutionAux
solution = soln
variable = P_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[vel_x]
type=CondactivityAux
variable=vel_x
comp_i=0
[../]
[vel_y]
type=CondactivityAux
variable=vel_y
comp_i=1
[../]
[vel_z]
type=CondactivityAux
variable=vel_z
comp_i=2
[../]
[]

[Variables]
[./CM] [../]
[]

[MeshModifiers]
 [./createNewSidesetX]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'left'
 boundary_id_new = 10
 block_id = 0
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.1  0.25 0.25'
 [../]

 [./createNewSidesetY]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'bottom'
 boundary_id_new = 11
 block_id = 0
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.25  0.1  0.25' 
[../]

 [./createNewSidesetZ]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'back'
 boundary_id_new = 12
 block_id = 0
 bottom_left = '-0.1  -0.1  -0.1'
 top_right =   '0.25  0.25  0.1'
[../]
[]


[Materials]
[./conductivity2] type = HydraulicConductivity3D
 fn = 9
 pressure=P_aux
 cond0=true
 cond1=false
 phi_m=0.1
 phi_f=0.9
 fx_string = '0.5,0.5,0.5,
              0.749975,0.75,0.749975,
              0.625,0.625,0.625'
 fy_string = '0.5,0.5,0.5,
              0.749975,0.749975,0.75,
              0.625,0.625,0.625'
 fz_string = '0.5,0.5,0.5,
              0.75,0.749975,0.749975,
              0.625,0.625,0.625'
 fa1_string = '0.0,0.0,0.0,
              0.0,0.0,0.0,
              0.0,0.0,0.0'
 fa2_string = '0.0,90.0,0.0,
               0.0,90.0,0.0,
               0.0,90.0,0.0'
 fa3_string = '0.0,0.0,90.0,
               0.0,0.0,90.0,
               0.0,0.0,90.0'
 fd1_string = '1.0,1.0,1.0,
               0.50005,0.50005,0.50005,
               0.2501,0.2501,0.2501'
 fd2_string = '1.0,1.0,1.0,
               0.50005,0.50005,0.50005,
               0.2501,0.2501,0.2501'
 fd3_string = '0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001'
 [../]
[]

[Kernels]
active='convection stab time'

[upwind]
type = AlgebraicDiffusion
variable = CM
[../]

[./convection]
type = Advection
variable = CM
int_by_parts=false
[../]

[./MyDiffusion]
type = MyDiffusion
variable = CM
coef=1.0e-9
[../]

[./stab]
type = AdvectionSUPG
variable = CM
coef=0.8
use_h=true
[../]

[./timestab]
type = TimeAdvectionSUPG
variable = CM
coef=0.5
use_h=true
lumping=false
[../]

[./time]
type = PorosityTimeDerivative
variable = CM
lumping = false
[../]

[./diff]
type = AnisotropicDiffusion
variable = CM
tensor_coeff='1.e-8 0 0 0 1.0e-8 0 0 0 1.0e-8'
[../]
[]

[BCs]
[./u_injection_left] type = DirichletBC boundary = '10 11 12' variable = CM value='1.0' [../]
#[./u_injection_right] type = OutflowBC boundary = 'outflow' variable = CM [../]
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

dt = 0.0025
num_steps=100

[./Quadrature]
order=SIXTH
[../]

[]

[Outputs]
 file_base = advectioOut1_1
exodus = true
csv=true
print_perf_log = true
[]


[UserObjects]
[./soln]
type = SolutionUserObject
mesh = DiffusionOutput1_1.e
timestep = 2
system_variables = pressure
execute_on = 'initial'
[../]
[]

[Postprocessors]

[./intC00]
type = ElementIntegral_Variable_Region
variable = CM region = 0 execute_on = 'timestep_end'
[../]
[./intD00]
type = ElementIntegral_Variable_Region
region = 0 execute_on = 'timestep_end'
[../]

[./intC01]
 type = ElementIntegral_Variable_Region
 variable = CM region = 1 execute_on = 'timestep_end'
 [../]
[./intD01]
 type = ElementIntegral_Variable_Region
 region = 1 execute_on = 'timestep_end'
 [../]

 [./intC02]
 type = ElementIntegral_Variable_Region
 variable = CM region = 2 execute_on = 'timestep_end'
 [../]
[./intD02]
 type = ElementIntegral_Variable_Region
 region = 2 execute_on = 'timestep_end'
 [../]

 [./intC03]
 type = ElementIntegral_Variable_Region
 variable = CM region = 3 execute_on = 'timestep_end'
 [../]
[./intD03]
 type = ElementIntegral_Variable_Region
 region = 3 execute_on = 'timestep_end'
 [../]

 [./intC04]
 type = ElementIntegral_Variable_Region
 variable = CM region = 4 execute_on = 'timestep_end'
 [../]
[./intD04]
 type = ElementIntegral_Variable_Region
 region = 4 execute_on = 'timestep_end'
 [../]

 [./intC05]
 type = ElementIntegral_Variable_Region
 variable = CM region = 5 execute_on = 'timestep_end'
 [../]
[./intD05]
 type = ElementIntegral_Variable_Region
 region = 5 execute_on = 'timestep_end'
 [../]

 [./intC06]
 type = ElementIntegral_Variable_Region
 variable = CM region = 6 execute_on = 'timestep_end'
 [../]
[./intD06]
 type = ElementIntegral_Variable_Region
 region = 6 execute_on = 'timestep_end'
 [../]

 [./intC07]
 type = ElementIntegral_Variable_Region
 variable = CM region = 8 execute_on = 'timestep_end'
 [../]
[./intD07]
 type = ElementIntegral_Variable_Region
 region = 7 execute_on = 'timestep_end'
 [../]

 [./intC08]
 type = ElementIntegral_Variable_Region
 variable = CM region = 8 execute_on = 'timestep_end'
 [../]
[./intD08]
 type = ElementIntegral_Variable_Region
 region = 8 execute_on = 'timestep_end'
 [../]

 [./intC09]
 type = ElementIntegral_Variable_Region
 variable = CM region = 9 execute_on = 'timestep_end'
 [../]
[./intD09]
 type = ElementIntegral_Variable_Region
 region = 9 execute_on = 'timestep_end'
 [../]

 [./intC10]
 type = ElementIntegral_Variable_Region
 variable = CM region = 10 execute_on = 'timestep_end'
 [../]
[./intD10]
 type = ElementIntegral_Variable_Region
 region = 10 execute_on = 'timestep_end'
 [../]

 [./intC11]
 type = ElementIntegral_Variable_Region
 variable = CM region = 11 execute_on = 'timestep_end'
 [../]
[./intD11]
 type = ElementIntegral_Variable_Region
 region = 11 execute_on = 'timestep_end'
 [../]

 [./intC12]
 type = ElementIntegral_Variable_Region
 variable = CM region = 12 execute_on = 'timestep_end'
 [../]
[./intD12]
 type = ElementIntegral_Variable_Region
 region = 12 execute_on = 'timestep_end'
 [../]

 [./intC13]
 type = ElementIntegral_Variable_Region
 variable = CM region = 13 execute_on = 'timestep_end'
 [../]
[./intD13]
 type = ElementIntegral_Variable_Region
 region = 13 execute_on = 'timestep_end'
 [../]

 [./intC14]
 type = ElementIntegral_Variable_Region
 variable = CM region = 14 execute_on = 'timestep_end'
 [../]
[./intD14]
 type = ElementIntegral_Variable_Region
 region = 14 execute_on = 'timestep_end'
 [../]

 [./intC15]
 type = ElementIntegral_Variable_Region
 variable = CM region = 15 execute_on = 'timestep_end'
 [../]
[./intD15]
 type = ElementIntegral_Variable_Region
 region = 15 execute_on = 'timestep_end'
 [../]

 [./intC16]
 type = ElementIntegral_Variable_Region
 variable = CM region = 16 execute_on = 'timestep_end'
 [../]
[./intD16]
 type = ElementIntegral_Variable_Region
 region = 16 execute_on = 'timestep_end'
 [../]

 [./intC17]
 type = ElementIntegral_Variable_Region
 variable = CM region = 17 execute_on = 'timestep_end'
 [../]
[./intD17]
 type = ElementIntegral_Variable_Region
 region = 17 execute_on = 'timestep_end'
 [../]

 [./intC18]
 type = ElementIntegral_Variable_Region
 variable = CM region = 18 execute_on = 'timestep_end'
 [../]
[./intD18]
 type = ElementIntegral_Variable_Region
 region = 18 execute_on = 'timestep_end'
 [../]

 [./intC19]
 type = ElementIntegral_Variable_Region
 variable = CM region = 19 execute_on = 'timestep_end'
 [../]
[./intD19]
 type = ElementIntegral_Variable_Region
 region = 19 execute_on = 'timestep_end'
 [../]

 [./intC20]
 type = ElementIntegral_Variable_Region
 variable = CM region = 20 execute_on = 'timestep_end'
 [../]
[./intD20]
 type = ElementIntegral_Variable_Region
 region = 20 execute_on = 'timestep_end'
 [../]

 [./intC21]
 type = ElementIntegral_Variable_Region
 variable = CM region = 21 execute_on = 'timestep_end'
 [../]
[./intD21]
 type = ElementIntegral_Variable_Region
 region = 21 execute_on = 'timestep_end'
 [../]

 
[]





