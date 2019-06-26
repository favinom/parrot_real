[Mesh]
 file = refinedMesh_${origLevel}_000${adapSteps}_mesh.xdr
[]

[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./myDiffusion] type = MyDiffusion variable = pressure coef=1.0 [../]
[]

[Materials]
[./conductivity1] 
 block='1 2 3 4 5 6 7 8'
type =  HydraulicConductivity
 conductivity = 1.0e4
[../]

[./conductivity2] 
 block ='11 12 13'
type =  HydraulicConductivity
conductivity = 1.0
[../]
[]

[Functions]
 [./fun_n]
 type = ParsedFunction
 value = '1.0'
 [../]
[]

[BCs]

[./fluxBC] type = NeumannBC variable = pressure value = 1.2752  boundary = '21' [../] #0_1_1.2650 #0_0 1.390 1.33750 1.252 1.415
[./dirBC]  type = DirichletBC variable = pressure value = 0  boundary = '22'  [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
 
[Postprocessors]
 active='average'
 [./average]
 type = SideAverageValue
 boundary = 21
 variable = pressure
 [../]
[./fluxBoundary]
 type = SideIntegralForFluxPostprocessor
 variable = pressure
 boundary   = '21'
 execute_on = 'timestep_end'
 [../]
 [./integral]
 type = AreaPostprocessor
 boundary = '21'
 [../]
 []
 
[Executioner]

 type=Steady
 solve_type= LINEAR
 line_search = none
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps'
 
# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'

[./Quadrature]
order=SIXTH
[../]

[]


[Outputs]
 file_base      = DiffusionOut_${origLevel}_${adapSteps}
 exodus         = true
 print_perf_log = true
 csv=true
 []

