[Mesh]
 file = test_3_blocchi.e #refineMesh_0_0002_mesh.xdr
 uniform_refine = 0
 #second_order=true
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
 block='1 2 3 4 5 6 7 8'
type =  HydraulicConductivity
 conductivity = 1.e4
[../]

[./conductivity2] 
 block ='11 12 13'
type =  HydraulicConductivity
conductivity = 1.0
[../]
[]
# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization

[Functions]
 [./fun_n]
 type = ParsedFunction
 value = '1*(z>0.33333300)*(z<0.66666600)'
 [../]
[]

[BCs]
[./dirBC]  type = DirichletBC variable = pressure value = 0  boundary = '22'  [../]
[./fluxBC] type = FunctionNeumannBC variable = pressure function = 'fun_n'  boundary = '21' [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
 
[Postprocessors]
 [./average]
 type = SideAverageValue
 boundary = 21
 variable = pressure
 [../]
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
order=SIXTH
[../]

[]


[Outputs]
 file_base      = DiffusionOutput_0
 exodus         = true
 print_perf_log = true
[]
