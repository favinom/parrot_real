#ifndef FractureAppConforming_H
#define FractureAppConforming_H

#include "MultiAppTransfer.h"
#include <utopia_fe.hpp>
#include <utopia_LibMeshBackend.hpp>
#include "utopia.hpp"
#include <memory>
#include "UserObjectInterface.h"
#include "StoreTransferOperators.h"
#include "utopia_fe_base.hpp"

class FractureAppConforming;

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
InputParameters validParams<FractureAppConforming>();

/**
 * Project values from one domain to another
 */
class FractureAppConforming : public GeneralUserObject

{
public:
    FractureAppConforming(const InputParameters & parameters);
    
    virtual void initialize() override;
    
    virtual void execute() override;
    
    virtual void finalize() override {}
    
    typedef utopia::DSMatrixd SparseMatT;
    
    typedef utopia::DVectord VecT;
    
    void stabilize_A_matrix(FEProblemBase & _problem, utopia::USparseMatrix &A_matrix, utopia::USparseMatrix &S_matrix);
    
    void assemble_poro_mass_matrix(FEProblemBase & problem, utopia::USparseMatrix &Mass_p, utopia::USparseMatrix &L_Mass_p);

        
    VariableName _m_var_name;

    std::vector<int> _vector_p;

    std::vector<Real> _vector_value;

    
    void boundary_constraint_vec(utopia::UVector &boundary, utopia::UVector &vec, bool flag);
    
    void constraint_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints);

    void unstabilize_coeffiecient(FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs);
    
    bool  _boundary, _constraint_m;
    
    utopia::USparseMatrix A_m_t;
    
    utopia::UVector x_m, x_f, lagr;
    
    utopia::UVector c_m, mass_c_m_old_dot;
    
    utopia::UVector rhs_f_t, rhs_f_c;
    
    utopia::UVector rhs_m_t, rhs_m_c, int_dof;
    
    std::shared_ptr<utopia::USparseMatrix> _A_t_store = NULL;
    
    utopia::USparseMatrix mass, mass_p;
    
    utopia::USparseMatrix mass_lumped, mass_lumped_p;
    
    utopia::UVector c_m_old;
    
    utopia::UVector diag_elem, diag_elem_p, c_dot;
    
    void zero_constraint_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);
    
    void one_constraint_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);
    
    std::string _userobject_name_1 = "Matrix_System";
    
    void CopyMatrixSolution(utopia::UVector _sol_m);

    void send_list(FEProblemBase & _problem, utopia::USparseMatrix &_M, std::vector<int> &a);

    void assemble_mass_matrix(FEProblemBase & _problem, utopia::USparseMatrix &mass_matrix, utopia::USparseMatrix &lumped_mass_matrix);

    void unconstraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);

    Real ComputeMaterialProprties(const Elem *elem);

    void stabilize_coeffiecient(utopia::UVector &rhs_m_c, FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs);



    
};

#endif /* FractureAppConforming */
