[Mesh]

  file = ../final${i}.xda
  uniform_refine = 0
  #second_order=true
[]

[MeshModifiers]
 [./createNewSidesetX]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'right'
 boundary_id_new = 10
 block_id = 0
 bottom_left = '0 0.87501 0.87501'
 top_right =   '1.00001  1.00001 1.00001'
 [../]
 
 [./createNewSidesetY]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'top'
 boundary_id_new = 11
 block_id = 0
 bottom_left = '0.87501 0  0.87501'
 top_right =   '1.00001  1.00001 1.00001'
 [../]
 
 [./createNewSidesetZ]
 type = AddSideSetsFromBoundingBox
 boundary_id_old = 'front'
 boundary_id_new = 12
 block_id = 0
 bottom_left = '0.87501  0.87501 0'
 top_right =   '1.00001  1.00001 1.00001'
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
[./conductivity2] type = HydraulicConductivity3D
 fn = 9
 cond0=true
 cond1=false
 phi_m=1.0
 phi_f=1.0
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

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization

[Functions]
 [./fun_dr]
 type = ParsedFunction
 value = '1*(x>0.875)*(y>0.875)*(z>0.875)'
 [../]
 [./fun_n]
 type = ParsedFunction
 value = '1*(x<0.2500)*(y<0.2500)*(z<0.2500)'
 [../]
[]

[BCs]
[./dirBC]  type = DirichletBC variable = pressure value = 1  boundary = '10 11 12'  [../]
[./fluxBC] type = FunctionNeumannBC variable = pressure function = fun_n  boundary = 'left bottom back' [../]
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
 petsc_options_value='  preonly   lu       NONZERO               superlu_dist'
 
# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'

[./Quadrature]
order=SIXTH
[../]

[]


[Outputs]
 file_base      = DiffusionOutput
 exodus         = true
 print_perf_log = true
[]
