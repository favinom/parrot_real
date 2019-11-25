#[Mesh]
#type = GeneratedMesh
#dim = 3
#nx = 50
#ny = 50
#nz = 5
#zmin = 0
#zmax = 100
#xmin = 0
#xmax = 1000
#ymin = 0
#ymax = 1000
#elem_type = HEX8
#[]
 

[Mesh]
 file = refineMesh_0_0004_mesh.xdr
#second_order=true
 []
 
 
[MeshModifiers]
active=''
 [./middle_node]
 type = AddExtraNodeset
 new_boundary = 10
 coord = '500 500 100'
 [../]
 [./subdomains]
 type = SubdomainBoundingBox
 bottom_left = '480 480 -1.0'
 block_id = 1
 top_right = '520 520 100'
 [../]
 []

 
 
[Variables]
[./pressure]
 order=FIRST
 family=LAGRANGE
 [../]
[]

[Kernels]
 active='myDiffusion time bf'
[./myDiffusion]
 type = MyDiffusion
 variable = pressure
 coef = 1.0
[../]
 
[./time]
 type = PorosityTimeDerivative
 variable = pressure
 lumping = true
[../]

 [./bf]
 type = BodyForce
 variable = pressure
 function ='0.002*(((x-500)*(x-500)+(y-500)*(y-500))<0.0144)'
 [../]
[]

[Materials]
 active='conductivity1 porosity1'
[./conductivity1]
 type = HydraulicConductivity
 conductivity = 4e-13
[../]
[./porosity1]
 type = Porosity
 phi = 7.0e-11
[../]
[./conductivity2]
 type = HydraulicConductivityFunction
 conductivity = 1.0
 function ='4.4372*(((x-500)*(x-500)+(y-500)*(y-500))<=0.0144)+2e-10*(((x-500)*(x-500)+(y-500)*(y-500))>0.0144)'
 [../]
[./porosity2]
 type = PorosityFunction
 phi = 0.01
 function = 1 #'0.0*(((x-500)*(x-500)+(y-500)*(y-500))<=0.0144)+7.2e-11*(((x-500)*(x-500)+(y-500)*(y-500))>0.0144)'
[../]
 
[]

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
active=''
[./outflowBC]
 type = FunctionNeumannBC
 variable = pressure
 function='7e1*(((x-500)*(x-500)+(y-500)*(y-500))<0.0144)'
 boundary = front  [../]
[]
 
[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]
 
[Executioner]

 type = Transient
 solve_type= LINEAR
 line_search = none
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='  preonly   lu       NONZERO               mumps '
 
# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'
dt = 10.0
num_steps=100
 [./Quadrature]
 order=SIXTH
 [../]
[]

 

[Outputs]
 file_base      = OutputBenchmark1
 exodus         = true
 print_perf_log = true
  xdr = true
[]
