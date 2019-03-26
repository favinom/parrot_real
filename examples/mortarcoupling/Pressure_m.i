#[Mesh]
#  type = GeneratedMesh
#  dim = 2
#  xmin =-0.5
#  xmax = 0.5
#  ymin =-0.5
#  ymax = 0.5
#  nx = 40
#  ny = 40
#  elem_type = QUAD9
#[]

[Mesh]
type = FileMesh
file = mesh_m_q3.e
second_order=true
[]


[Variables]
[./P]
order = SECOND
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
 boundary = right
 value=4
[../]
 
[./_b_21]
 type = DirichletBC
 variable = P
 boundary = left
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



