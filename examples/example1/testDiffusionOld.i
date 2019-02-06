[Mesh]
  file = Mesh_level1.e
  block_id = '2 4 5 6 7'
  boundary_id = '1 2'
  boundary_name = 'inflow outflow'
  uniform_refine = 0
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
[]

[Materials]
[./conductivity2] type = HydraulicConductivity block = 2 conductivity = 20.0 [../]
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
 petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
 petsc_options_value='   preonly   lu       NONZERO               mumps'

 [./Quadrature]
 order=SIXTH
 [../]

[]


[Outputs]
 file_base= OutputDiffusion1
 exodus = true
[]
