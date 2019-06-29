#ifndef FractureAppNConforming_H
#define FractureAppNConforming_H

#include "MultiAppTransfer.h"
#include <utopia_fe.hpp>
#include <utopia_LibMeshBackend.hpp>
#include "utopia.hpp"
#include <memory>
#include "UserObjectInterface.h"
#include "StoreTransferOperators.h"
#include "utopia_fe_base.hpp"

class FractureAppNConforming;

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
InputParameters validParams<FractureAppNConforming>();

/**
 * Project values from one domain to another
 */
class FractureAppNConforming : public GeneralUserObject

{
public:
    FractureAppNConforming(const InputParameters & parameters);

    virtual void initialize() override;
    
    virtual void execute() override;
    
    virtual void finalize() override {}
  
    typedef utopia::DSMatrixd SparseMatT;
  
    typedef utopia::DVectord VecT;

protected:



    VariableName _f_var_name, _m_var_name;
    
    bool ok = false;

    bool _impact;

    // userobject to store our operators
    const StoreTransferOperators & _operator_storage;

    std::string _operator_type;
    
    virtual void CopyMatrixSolution(utopia::UVector to_sol);
    
    virtual void CopyFractureSolution(utopia::UVector to_sol);
    
    bool solve();
    
    bool solve_cg_dual();
    
    bool solve_monolithic();

    void solve_transport_monolithic();
    
    void solve_transport_stabilize();

    void assemble_mass_matrix(double _porosity, FEProblemBase & problem, utopia::USparseMatrix &Mass, utopia::USparseMatrix &L_Mass);

    void assemble_poro_mass_matrix(double _porosity, FEProblemBase & problem, utopia::USparseMatrix &Mass_p, utopia::USparseMatrix &L_Mass_p);
    
    void constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool flag);

    void constraint_concentration_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints);
    
    bool _solve_cg, _solve_mg, _pressure, _transport, _boundary, _biorth, _stabilize, _constraint_m, _constraint_f;

    std::vector<int> _vector_p;

    std::vector<Real> _vector_value;

    double _porosity_m, _porosity_f;

    std::string _multiapp_name;
    
    utopia::USparseMatrix D, B, D_t, B_t, A_m_t,A_f_t; //, A_t2;

    utopia::USparseMatrix D_temp, B_temp;
    
    utopia::UVector x_m, x_f, lagr;
    
    utopia::UVector c_m, c_f, mass_c_m_old_dot, mass_c_f_old_dot, lagr_t, dlagr_t, x_tot;
    
    utopia::UVector rhs_f_t, rhs_f_c;
    
    utopia::UVector rhs_m_t, rhs_m_c, int_dof;

    std::shared_ptr<utopia::USparseMatrix> _A_t_store = NULL;

    std::shared_ptr<utopia::USparseMatrix> _S_m_store = NULL;

    std::shared_ptr<utopia::USparseMatrix> _S_f_store = NULL;

    utopia::USparseMatrix mass_m, mass_mp, mass_f, mass_fp;

    utopia::USparseMatrix mass_lumped_mp, mass_lumped_fp, mass_lumped_m, mass_lumped_f;

    utopia::UVector c_m_old, c_f_old;

    std::shared_ptr<SparseMatT> _P= NULL;  

    utopia::UVector diag_elem, diag_elem_p, c_dot_m, c_dot_f;
    
    void bundary_volume_permulation_matrix(EquationSystems &_b_es, unsigned int _b_var_num, MooseVariable & _v_var, std::shared_ptr<SparseMatT> P_matrix);

    void stabilize_A_matrix(FEProblemBase & _problem, utopia::USparseMatrix &S_matrix);

    void stabilize_coeffiecient(FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs);
    
    void unstabilize_coeffiecient(FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs);

    void send_list(FEProblemBase & _problem, utopia::USparseMatrix &_M, std::vector<int> &a);

    void unconstraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);

    void one_constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints);

    void solve_cg_dual_stabilize();

    void compute_transport_monolithic_static_condenstation();

    void stabilize_A_matrix_cons(utopia::USparseMatrix &A_0, utopia::USparseMatrix &S_matrix);

    Real ComputeMaterialProprties(const Elem *elem);

    std::string _userobject_name_1 = "Matrix_m";

    std::string _userobject_name_2 = "Matrix_f";

    std::string _userobject_name_T = "Matrix_t";

    std::string _userobject_name_A_m = "Matrix_Am";

    std::string _userobject_name_A_f = "Matrix_Af";

    std::string _userobject_name_B = "Matrix_B";

    std::string _userobject_name_D = "Matrix_D";

    
    std::unique_ptr<ExodusII_IO> _ex_writer;

};

#endif /* FractureAppConforming */
