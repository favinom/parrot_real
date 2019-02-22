[Mesh]
 file = adaptedAMR1_0006_mesh.xda
[]

[MeshModifiers]
 [./createNewSidesetX]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'bottom'
 boundary_id_new = 21
 block_id = 0
 bottom_left = '-0.1 -0.1 0.33333'
 top_right =   '1.1  0.1  0.66667'
 [../]
 
 [./createNewSidesetY]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'top'
 boundary_id_new = 22
 block_id = 0
 bottom_left = '-0.1  2.220  -0.1'
 top_right =   ' 1.1  2.30  0.334'
 [../]
 
 [./createNewSidesetZ]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'top'
 boundary_id_new = 23
 block_id = 0
 bottom_left = '-0.1  2.220  0.6665'
 top_right =   ' 1.1  2.30  1.1'
 [../]
[]

[Variables]
[./pressure] order=FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./myDiffusion] type = MyDiffusion variable = pressure coef=1.0 [../]
[]

 [AuxVariables]
[./pp] order=CONSTANT  family=MONOMIAL [../]
 []
 
 [AuxKernels]
[./ciao] type = MaterialRealTensorValueAux i = 0 j = 0 variable = pp property = conductivityTensor [../]
 []
 
[Materials]
[./conductivity2] type =  HydraulicConductivity3D
 cond0=true
 cond1=false
 phi_m=1.0
 phi_f=1.0
 fn = 8
 fx_string = '0.5  ,0.5  ,0.77,0.83, 0.2,0.2  , 0.5,0.5'
 fy_string = '1.125,0.175,2.05,2.05, 2.05,2.05 , 1.6,1.6'
 fz_string = '0.5  ,0.5  ,0.5 ,0.5, 0.5,0.5 , 0.675,0.31'
 fa1_string = '0,90,90,90,78.6901,-78.6901,0,0'
 fa2_string = '0, 0, 0, 0,0,0,0,0'
 fa3_string = '0,90,90,90,90,90,16.2602,-15.8192'
 fd1_string = '0.9,0.25,0.3,0.3,0.3059,0.3059,0.9,0.9'
 fd2_string = '1.75,0.9,0.4,0.4,0.4,0.4,1.25,1.2472'
 fd3_string = '0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01'
 [../]
[]

[BCs]
[./dirBC22]  type = DirichletBC variable = pressure value = 0  boundary = '22'  [../]
[./dirBC23]  type = DirichletBC variable = pressure value = 0  boundary = '23'  [../]
[./fluxBC] type = NeumannBC   variable = pressure value = 1  boundary = '21'  [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]

 type=Steady
 solve_type= newton
 line_search = none
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps'
 
# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'

[./Quadrature]
order=SEVENTH
[../]

[]


[Outputs]
 file_base      = DiffusionOutput
 exodus         = true
 print_perf_log = true
[]
