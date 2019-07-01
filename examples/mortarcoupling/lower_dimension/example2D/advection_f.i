[Mesh]
type = FileMesh
file = mesh_fracture_2d_non_conf.e
uniform_refine=0
[]
 
 
[Variables]
[./cm]
order = FIRST
family = LAGRANGE
[../]
[]

[AuxVariables]
[./lambda]
[../]
[./p_aux]
[../]
[]

[Kernels]
[./convection]
type = Advection
variable = cm
p = p_vel
int_by_parts=false
[../]
[]

[Materials]
[./conductivity1]
 type =  HydraulicConductivity
 conductivity = 1e-1
 pressure = p_aux
[../]
[]

[Executioner]
 type = Transient
 solve_type = LINEAR
 dt = 1e7
 num_steps=100.0
 l_tol = 1E-14
 []

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]


[Problem]
type = FEProblem
solve = false
kernel_coverage_check = false
[]

[AuxKernels]
#active=''
[./en]
type = SolutionAux
solution = soln
variable = p_aux
scale_factor = 1.0
execute_on = 'initial'
[../]
[]

[UserObjects]
#active=''
[./soln]
type = SolutionUserObject
mesh = fracture_p.e
timestep = 1
system_variables = p
execute_on = 'initial'
[../]
[]


[Outputs]
exodus = true
execute_on = 'timestep_end'
[]
