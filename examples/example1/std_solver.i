# [Problem] type = ParrotProblem []


[Executioner]

type = Transient
solve_type= LINEAR
line_search = none
petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
petsc_options_value='   preonly   lu      NONZERO               mumps                       '
 
dt = 1e7
num_steps=100

[./Quadrature]
order=SIXTH
[../]

[]
