#ifndef SystemInitialize_H
#define SystemInitialize_H

#include "MultiAppTransfer.h"
#include <utopia_fe.hpp>
#include <utopia_LibMeshBackend.hpp>
#include "utopia.hpp"
#include <memory>
#include "UserObjectInterface.h"
#include "StoreTransferOperators.h"
#include "utopia_fe_base.hpp"

class SystemInitialize
;

/*
 *                                  _.._
 *                                .' .-'`
 *                               /  /
 *                               |  |
 *                               \  '.___.;
 *                                '._  _.'
 *
 *                 ) \     /_(
 *                  )_`-)-_(_
 *                   `,' `__,)
 *                  _/   ((
 *         ______,-'    )
 *        (            ,
 *         \  )    (   |
 *        /| /`-----` /|
 *        \ \        / |
 *        |\|\      /| |\
 */

template <>
InputParameters validParams<SystemInitialize>();

/**
 * Project values from one domain to another
 */
class SystemInitialize : public GeneralUserObject

{
public:
    SystemInitialize(const InputParameters & parameters);
    
    virtual void initialize() override;
    
    virtual void execute() override;
    
    virtual void finalize() override {}

    VariableName _m_var_name;

    void stabilize_A_matrix(FEProblemBase & _problem, utopia::USparseMatrix &A_0, utopia::USparseMatrix &S_matrix);

    void constraint_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints);

    void unconstraint_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints);
    
    void CopyMatrixSolution(utopia::UVector _sol_m);

    void boundary_constraint_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);

    void unconstraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);

    void assemble_mass_matrix(FEProblemBase & _problem, utopia::USparseMatrix &mass_matrix, utopia::USparseMatrix &lumped_mass_matrix);

    void assemble_rhs(utopia::USparseMatrix &A_matrx, utopia::UVector &rhs);

    
};

#endif /* FractureAppConforming */
