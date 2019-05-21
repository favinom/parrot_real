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

protected:

    // std::string _name_mesh;

    // std::vector<BoundaryID> _slave_id, _master_id;

    const MaterialProperty<Real> &_poro;

    const int _slave_id;

    VariableName _f_var_name, _m_var_name;
    
    //AuxVariableName _l_var_name;
  
    //VariableName _from_var_name;
    
    bool ok = false;

    bool _impact;

    // userobject to store our operators
    const StoreTransferOperators & _operator_storage;

    std::string _operator_type;
    
    virtual void CopyMatrixSolution(utopia::UVector to_sol);
    
    virtual void CopyFractureSolution(utopia::UVector to_sol);
    
    bool solve();
    
    bool solve_cg_dual();
    
    bool solve_monolithic();\
    
    void solve_stabilize_monolithic();

    void assemble_mass_matrix(double _porosity, FEProblemBase & problem, utopia::USparseMatrix &Mass, utopia::USparseMatrix &L_Mass);

    void assemble_poro_mass_matrix(double _porosity, FEProblemBase & problem, utopia::USparseMatrix &Mass_p, utopia::USparseMatrix &L_Mass_p);
    
    void constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool flag);

    void constraint_concentration_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints);
    
    bool _solve_cg, _solve_mg, _pressure, _transport, _boundary, _constraint_m, _constraint_f;

    double _porosity_m, _porosity_f;

    std::string _multiapp_name;
    
    utopia::USparseMatrix D, B, D_t, B_t, A_m_t,A_f_t; //, A_t2;

    utopia::USparseMatrix D_temp, B_temp;
    
    utopia::UVector x_m, x_f, lagr;
    
    utopia::UVector c_m, c_f, mass_c_m_old_dot, mass_c_f_old_dot, lagr_t, dlagr_t, x_tot;
    
    utopia::UVector rhs_f_t, rhs_f_c;
    
    utopia::UVector rhs_m_t, rhs_m_c, int_dof;

    std::shared_ptr<utopia::USparseMatrix> _A_t_store = NULL;

    utopia::USparseMatrix mass_m, mass_mp, mass_lumrhs_f_tdp;

    utopia::USparseMatrix mass_lumpedp, mass_lumped;

    utopia::UVector c_m_old, c_f_old;

    std::shared_ptr<SparseMatT> _P= NULL;  

    utopia::UVector diag_elem, diag_elem_p, c_dot;
    
    void bundary_volume_permulation_matrix(EquationSystems &_b_es, unsigned int _b_var_num, MooseVariable & _v_var, std::shared_ptr<SparseMatT> P_matrix);

    void stabilize_A_matrix(FEProblemBase & _problem, utopia::USparseMatrix &S_matrix);

    void stabilize_coeffiecient(FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs);
    
    void unstabilize_coeffiecient(FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs);

    void send_list(FEProblemBase & _problem, utopia::USparseMatrix &_M, std::vector<int> &a);

    void unconstraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);

    void one_constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);

};

#endif /* FractureAppConforming */
