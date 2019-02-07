
[Mesh]
type = FileMesh
file = Mesh_level3.e
boundary_id = '1 2'
boundary_name = 'inflow outflow'
[]

[Variables]
[./CM] [../]
[]

[AuxVariables]
[./P_aux]
[../]
[]

[Kernels]

[./convection]
type = Advection # ConservativeAdvection #
variable = CM
[../]

[./stab]
type = AdvectionSUPG # ConservativeAdvection #
variable = CM
coef=0.5
use_h=true
[]

[./time]
type = TimeDerivative
variable = CM
[../]

[]

[Materials]
[./conductivity2] type = HydraulicConductivityFracture block = 2 conductivity = '1e-3 1e-3 20.0' theta = '0.0 30.9638 0.0' pressure=P_aux[../]
[./conductivity4] type = HydraulicConductivity block = 4 conductivity = 1e-6 pressure=P_aux [../]
[./conductivity5] type = HydraulicConductivity block = 5 conductivity = 1e-6 pressure=P_aux [../]
[./conductivity6] type = HydraulicConductivity block = 6 conductivity = 1e-6 pressure=P_aux [../]
[./conductivity7] type = HydraulicConductivity block = 7 conductivity = 1e-5 pressure=P_aux[../]

[]

[Executioner]

type = Transient
solve_type= LINEAR
line_search = none
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='   preonly   lu       NONZERO               mumps'

# petsc_options_iname = '-pc_type -pc_hypre_type'
# petsc_options_value = 'hypre boomeramg'

dt = 0.01
num_steps=1

[./Quadrature]
order=SIXTH
[../]

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
variable = P_aux
scale_factor = 1.0
execute_on ='initial'
[../]
[]


[UserObjects]
[./soln]
type = SolutionUserObject
mesh = OutputBenchmark1.e
timestep = 2
system_variables = pressure
execute_on = 'initial'
[../]
[]

[Postprocessors]

[./flux_left]
type = SideFlux
variable = CM
boundary = outflow
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
print_perf_log = true
[]



