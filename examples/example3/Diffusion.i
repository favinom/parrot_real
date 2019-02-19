[Mesh]
 file = mesh_test_3.e
 #uniform_refine = 0
 #second_order=true
[]

[MeshModifiers]
# [./createNewSidesetX]
# type = AddSideSetsFromBoundingBox
# boundary_id_old = 3
# boundary_id_new = 10
# block_id = 2
# bottom_left = '-10 -10 -10'
# top_right =   '10  10 10'
# [../]
 
# [./createNewSidesetY]
# type = AddSideSetsFromBoundingBox
# boundary_id_old = '2'
# boundary_id_new = 11
# block_id = 2
# bottom_left = '-0.10 2.249 -0.10'
# top_right =   '1.0  2.2510 0.3333'
# [../]
 
# [./createNewSidesetZ]
# type = AddSideSetsFromBoundingBox
# boundary_id_old = '2'
# boundary_id_new = 12
# block_id = 2
# bottom_left = '-0.10 2.249 0.6666'
# top_right =   '1.0   2.251 1.0'
#[../]
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
[./conductivity1] 
block='1' 
type =  HydraulicConductivity
 conductivity = 1e4 
[../]

[./conductivity2] 
block ='2' 
type =  HydraulicConductivity
conductivity = 1.0
[../]
[]
# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization

[Functions]
 [./fun_n]
 type = ParsedFunction
 value = '1*(x<0.2500)*(y<0.2500)*(z<0.2500)'
 [../]
[]

[BCs]
[./dirBC]  type = DirichletBC variable = pressure value = 1  boundary = 11  [../]
[./fluxBC] type = NeumannBC variable = pressure value = '1'  boundary = '10' [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]

 type=Steady
 solve_type= newton
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
 file_base      = DiffusionOutput1
 exodus         = true
 print_perf_log = true
[]
