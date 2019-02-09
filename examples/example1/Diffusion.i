[Mesh]
  file = Mesh_level3_hex_finer.e
  block_id = '2 4 5 6 7'
  boundary_id = '1 2'
  boundary_name = 'inflow outflow'
  uniform_refine = 0
  #second_order=true
[]


[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[MeshModifiers]
  [./rotate]
    type = Transform
    transform = TRANSLATE
    vector_value = '50 50 50'
  [../]
[]
 
[AuxVariables]
[./k_00]
order = CONSTANT
family = MONOMIAL
[../]
[./k_01]
order = CONSTANT
family = MONOMIAL
[../]
[./k_02]
order = CONSTANT
family = MONOMIAL
[../]
[./k_10]
order = CONSTANT
family = MONOMIAL
[../]
[./k_11]
order = CONSTANT
family = MONOMIAL
[../]
[./k_12]
order = CONSTANT
family = MONOMIAL
[../]
[./k_20]
order = CONSTANT
family = MONOMIAL
[../]
[./k_21]
order = CONSTANT
family = MONOMIAL
[../]
[./k_22]
order = CONSTANT
family = MONOMIAL
[../]
[]

[Kernels]
[./myDiffusion] type = MyDiffusion variable = pressure  [../]
[]

[AuxKernels]
[./auxfiber_00]
fracture=true
type = CondactivityAux
variable = k_00
comp_i=0
comp_j=0
execute_on = timestep_end
[../]

[./auxfiber_01]
fracture=true
type = CondactivityAux
variable = k_01
comp_i=0
comp_j=1
execute_on = timestep_end
[../]

[./auxfiber_02]
fracture=true
type = CondactivityAux
variable = k_02
comp_i=0
comp_j=2
execute_on = timestep_end
[../]

[./auxfiber_10]
fracture=true
type = CondactivityAux
variable = k_10
comp_i=1
comp_j=0
execute_on = timestep_end
[../]

[./auxfiber_11]
fracture=true
type = CondactivityAux
variable = k_11
comp_i=1
comp_j=1
execute_on = timestep_end
[../]

[./auxfiber_12]
fracture=true
type = CondactivityAux
variable = k_12
comp_i=1
comp_j=2
execute_on = timestep_end
[../]

[./auxfiber_20]
fracture=true
type = CondactivityAux
variable = k_20
comp_i=2
comp_j=0
execute_on = timestep_end
[../]

[./auxfiber_21]
fracture=true
type = CondactivityAux
variable = k_21
comp_i=0
comp_j=1
execute_on = timestep_end
[../]

[./auxfiber_22]
fracture=true
type = CondactivityAux
variable = k_22
comp_i=2
comp_j=2
execute_on = timestep_end
[../]
[]


[Materials]
#[./conductivity2] type = HydraulicConductivityFracture block = 2 conductivity = '1e-3 1e-3 20.0' theta = '0.0 30.9638 0.0' [../]
#[./conductivity4] type = HydraulicConductivity block = 4 conductivity = 1e-6 [../]
#[./conductivity5] type = HydraulicConductivity block = 5 conductivity = 1e-6 [../]
#[./conductivity6] type = HydraulicConductivity block = 6 conductivity = 1e-6 [../]
#[./conductivity7] type = HydraulicConductivity block = 7 conductivity = 1e-5 [../]
[./conductivity2] type = HydraulicConductivityFracture block = 2 conductivity = '1e-1 1e-1 10' theta = '0.0 30.9638 0.0' [../]
[./conductivity4] type = HydraulicConductivity block = 4 conductivity = 1e-6 [../]
[./conductivity5] type = HydraulicConductivity block = 5 conductivity = 1e-6 [../]
[./conductivity6] type = HydraulicConductivity block = 6 conductivity = 1e-6 [../]
[./conductivity7] type = HydraulicConductivity block = 7 conductivity = 1e-5 [../]
[]

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization

[BCs]
[./inflowBC]  type = DirichletBC variable = pressure value = 4.0  boundary = inflow  [../]
[./outflowBC] type = DirichletBC variable = pressure value = 1.0  boundary = outflow [../]
#[./inflowBC]  type = PenaltyDirichletBC variable = pressure boundary = inflow  value = 4.0 penalty = 1e10 [../]
#[./outflowBC] type = PenaltyDirichletBC variable = pressure boundary = outflow value = 1.0 penalty = 1e10 [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]

 type=Steady
 solve_type= LINEAR
 line_search = none
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='   preonly   lu  NONZERO               mumps'
# petsc_options_value='   preonly   lu       NONZERO               mumps'

 
# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'
 
# l_tol = 1e-8
 
[./Quadrature]
order=SIXTH
[../]

[]


[Outputs]
 file_base = OutputBenchmark1_hex_fine
 exodus    = true
 print_perf_log=true
[]
