#ifndef FractureNetworkUserObject_H
#define FractureNetworkUserObject_H

#include "MultiAppTransfer.h"
#include <utopia_fe.hpp>
#include <utopia_LibMeshBackend.hpp>
#include "utopia.hpp"
#include <memory>
#include "UserObjectInterface.h"
#include "StoreTransferOperators.h"
#include "utopia_fe_base.hpp"

class FractureNetworkUserObject;

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
InputParameters validParams<FractureNetworkUserObject>();

/**
 * Project values from one domain to another
 */
class FractureNetworkUserObject : public GeneralUserObject

{
public:
  FractureNetworkUserObject(const InputParameters & parameters);

    virtual void initialize() override;
    
    virtual void execute() override;
    
    virtual void finalize() override {}
  
    typedef utopia::DSMatrixd SparseMatT;
  
    typedef utopia::DVectord VecT;

protected:

    // std::string _name_mesh;

    // std::vector<BoundaryID> _slave_id, _master_id;
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
    
    bool solve_monolithic();
    
   // void setup_transport_monolithic();
    
    void solve_transport_monolithic();
    
    void assemble_mass_matrix(double _porosity, FEProblemBase & problem, utopia::USparseMatrix &Mass);
    
    void constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool flag);

    void constraint_concentration_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints);
    
    bool _solve_cg, _pressure, _transport, _boundary, _constraint_m, _constraint_f;

    double _porosity_m, _porosity_f;

    std::string _multiapp_name;
    
    utopia::USparseMatrix D, B, D_t, B_t, A_m_t, A_f_t; //, A_t2;

    utopia::USparseMatrix D_temp, B_temp;
    
    utopia::UVector x_m, x_f, lagr;
    
    utopia::UVector c_m, c_f, mass_c_m_old_dot, mass_c_f_old_dot, lagr_t, dlagr_t, x_tot;
    
    utopia::UVector rhs_f_t, rhs_f_c;
    
    utopia::UVector rhs_m_t, rhs_m_c;

    std::shared_ptr<utopia::USparseMatrix> _A_t_store = NULL;

    utopia::USparseMatrix mass_m;

    utopia::USparseMatrix mass_f;

    utopia::UVector c_m_old, c_f_old;

    std::shared_ptr<SparseMatT> _P= NULL;  
    
    void bundary_volume_permulation_matrix(EquationSystems &_b_es, unsigned int _b_var_num, MooseVariable & _v_var, std::shared_ptr<SparseMatT> P_matrix);
    //utopia::USparseMatrix A_m_tot, A_f_tot;

    
};

#endif /* L2ProjectionLibMeshTransferOld_H */
