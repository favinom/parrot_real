[Mesh]
  file = MeshBenchmark1.e
  block_id = '2 4 5 6 7'
  boundary_id = '1 2'
  boundary_name = 'inflow outflow'
  uniform_refine = 2
[]

#[MeshModifiers]
#  [./scale]
#    type = Transform
#    transform = SCALE
#    vector_value = '0.8474 0.8474 0.8474'
#  [../]
#[]

[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]
 
[Kernels]
[./myDiffusion] type = MyDiffusion variable = pressure  [../]
[]

 
[BCs]
[./inflowBC]  type = DirichletBC variable = pressure value = 4.0  boundary = inflow  [../]
[./outflowBC] type = DirichletBC variable = pressure value = 1.0  boundary = outflow [../]
#[./inflowBC]  type = PenaltyDirichletBC variable = pressure boundary = inflow  value = 4.0 penalty = 1e10 [../]
#[./outflowBC] type = PenaltyDirichletBC variable = pressure boundary = outflow value = 1.0 penalty = 1e10 [../]
[]

[Materials]
#[./conductivity2] type = HydraulicConductivityFracture block = 2 conductivity = '1e-3 1e-3 20.0' theta = '0.0 30.9638 0.0' [../]
#[./conductivity4] type = HydraulicConductivity block = 4 conductivity = 1e-6 [../]
#[./conductivity5] type = HydraulicConductivity block = 5 conductivity = 1e-6 [../]
#[./conductivity6] type = HydraulicConductivity block = 6 conductivity = 1e-6 [../]
#[./conductivity7] type = HydraulicConductivity block = 7 conductivity = 1e-5 [../]
[./conductivity2] type = HydraulicConductivityFracture block = 2 conductivity = '1e-3 1e-3 20.0' theta = '0.0 30.9638 0.0' [../]
[./conductivity4] type = HydraulicConductivity block = 4 conductivity = 1e-6 [../]
[./conductivity5] type = HydraulicConductivity block = 5 conductivity = 1e-6 [../]
[./conductivity6] type = HydraulicConductivity block = 6 conductivity = 1e-6 [../]
[./conductivity7] type = HydraulicConductivity block = 7 conductivity = 1e-5 [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]

 type=Steady
 solve_type= LINEAR
 line_search = none
#petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
# petsc_options_value='   preonly   lu       NONZERO               mumps'
#petsc_options_value='   preonly   cholesky  NONZERO               mumps'
 
 petsc_options_iname = '-pc_type -pc_hypre_type'
 petsc_options_value = 'hypre boomeramg'
 
 l_tol = 1e-8
 
#[./Quadrature]
#order=TENTH
#[../]

[]


[Outputs]
 file_base = OutputBenchmark1
 exodus    = true
[]
