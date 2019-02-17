[Mesh]
file = adapt1.xda
#uniform_refine = 0
#partitioner = parmetis
[]

[Variables]
[./pressure]  order = FIRST  family=LAGRANGE [../]
[]

[Kernels]
[./StressDivergenceParrot_real_x]
 type = Reaction
 variable = pressure
[../]
[]

[Executioner]
 type=Steady
 solve_type=LINEAR
 line_search = 'none'
[]

 
 
[Problem]
 type = FEProblem
 solve = false
[]

[Outputs]
 file_base = adaptedOut
 print_linear_residuals = true
 print_perf_log = true
 exodus = true
[]

