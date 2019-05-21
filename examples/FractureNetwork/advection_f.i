[Mesh]
type = FileMesh
file = inclined_line.e
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
[./p_vel]
[../]
[]

[Kernels]
[./convection]
type = Advection_Pressure
variable = CM
conductivity=0.1
p = p_vel
int_by_parts=false
[../]
[]


[Executioner]
type = Transient
solve_type = LINEAR
dt = 0.1
num_steps=1.0
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
variable = p_vel
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
