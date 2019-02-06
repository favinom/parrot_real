
[Mesh]
type = FileMesh
file = Mesh_level3.e
boundary_id = '1 2'
boundary_name = 'inflow outflow'
[]

[Variables]
[./CM]
[../]
[]

[AuxVariables]
[./P_vel]
[../]
[]

[Kernels]
[./convection]
type = AdvectionSUPG
variable = CM
p = P_vel
epsilon=1.0
[../]

[./time]
type = TimeDerivative
variable = CM
[../]

[]


[Materials]
#[./conductivity2] type = HydraulicConductivityFracture block = 2 conductivity = '1e-3 1e-3 20.0' theta = '0.0 30.9638 0.0' [../]
#[./conductivity4] type = HydraulicConductivity block = 4 conductivity = 1e-6 [../]
#[./conductivity5] type = HydraulicConductivity block = 5 conductivity = 1e-6 [../]
#[./conductivity6] type = HydraulicConductivity block = 6 conductivity = 1e-6 [../]
#[./conductivity7] type = HydraulicConductivity block = 7 conductivity = 1e-5 [../]
[./conductivity2] type = HydraulicConductivityFracture block = 2 conductivity = '1e-3 1e-3 20.0' theta = '0.0 30.9638 0.0' [../]
[./conductivity4] type = HydraulicConductivity block = 4 conductivity = 1e-6 [../]
[./conductivity5] type = HydraulicConductivity block = 5 conductivity = 1e-6 [../]
[./conductivity6] type = HydraulicConductivity block = 6 conductivity = 1e-6 [../]
[./conductivity7] type = HydraulicConductivity block = 7 conductivity = 1e-5 [../]
[]


[Executioner]
type = Transient
dt = 1.0e9
start_time = 0.0
end_time = 1.e9
###  nota che il PJFNK fa un secondo step che fa e piu lento
#    solve_type = 'PJFNK'
solve_type = 'NEWTON'
line_search = 'none'
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='   preonly   lu       NONZERO               mumps                '
[]


[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]



[BCs]
[./u_injection_left]
 type = DirichletBC
 boundary = inflow
 variable = CM
 value='0.01'
[../]
[]



[Problem]
type = FEProblem
solve = true
kernel_coverage_check = false
[]


[AuxKernels]
#active=''
[./en]
type = SolutionAux
solution = soln
variable = P_vel
scale_factor = 1.0
execute_on ='timestep_begin'
[../]
[]


[UserObjects]
#active='soln'
[./soln]
type = SolutionUserObject
mesh = OutputBenchmark1.e
timestep = 2
system_variables = pressure
execute_on = 'timestep_begin'
[../]
[]

[Functions]
[./ic_func]
type = ParsedFunction
value = '10*(x<0.9)*(x>0.7)*(y>0.7)*(y<0.9)'
[../]
[]

[Outputs]
exodus = true
execute_on = 'timestep_end'
[]


