[Mesh]
type = FileMesh
file = mesh_m.e
[]

[Variables]
[./P]
order = FIRST
family = LAGRANGE
[../]
[]


[AuxVariables]
[./p_slave]
[../]
[]


[Kernels]
[./diff_1]
type = Diffusion
variable = P
[../]
[]

[BCs]
# active=''
[./_b_22]
 type = DirichletBC
 variable = P
 boundary = 22
 value=4
[../]
 
[./_b_21]
 type = DirichletBC
 variable = P
 boundary = 21
 value=1
[../]
[]




[Problem]
 type = FEProblem
 solve = false
 kernel_coverage_check = false
[]


[Preconditioning]
 [./SMP]
 type = SMP
 full = true
 ksp_norm = default
 [../]
[]

[Executioner]
type = Transient
dt = 1.0
start_time = 0.0
end_time = 1.0
solve_type ='LINEAR'
[]

[Outputs]
file_base = matrix
exodus = true
execute_on = 'timestep_end'
[]



