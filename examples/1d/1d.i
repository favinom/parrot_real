[Problem]
type = ParrotProblem
[]

 [Mesh]
 type = GeneratedMesh
 dim = 1
 nx = 100
 xmin = 0.0
 xmax = 1.0
 []

 [Functions]
 [./parsed_function]
 type = ParsedFunction
 value = '-x'
 [../]
 []
 
[ICs]
 [./u_ic]
 type = FunctionIC
 variable = 'P_aux'
 function = parsed_function
 [../]
 []
 
[AuxVariables]
[./P_aux] [../]
[]

 

[Variables]
[./CM] [../]
[]

[Materials]
[./conductivity4] type = HydraulicConductivity   conductivity = 1   pressure = P_aux [../]


[./porosity2] type = Porosity phi = 1  [../]

[]


[Kernels]

[upwind]
type = AdvectionBubble
variable = CM
[../]
 
[]


[BCs]
[./u_injection_left] type = DirichletBC boundary = left variable = CM value='0.01' [../]
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

dt = 0.1
num_steps=100

[./Quadrature]
order=SIXTH
[../]

[]

[Outputs]
exodus = true
csv=true
print_perf_log = true
[]

