#ifndef Test_Assembly_H
#define Test_Assembly_H

#include "MultiAppTransfer.h"
#include <utopia_fe.hpp>
#include <utopia_LibMeshBackend.hpp>
#include "utopia.hpp"
#include <memory>
#include "UserObjectInterface.h"
#include "StoreTransferOperators.h"
#include "utopia_fe_base.hpp"

class Test_Assembly;

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
InputParameters validParams<Test_Assembly>();

/**
 * Project values from one domain to another
 */
class Test_Assembly : public GeneralUserObject

{
public:
    Test_Assembly(const InputParameters & parameters);

    virtual void initialize() override;
    
    virtual void execute() override;
    
    virtual void finalize() override {}
  
    typedef utopia::DSMatrixd SparseMatT;
  
    typedef utopia::DVectord VecT;

protected:

   

    void assemble_mass_matrix(FEProblemBase & problem, utopia::USparseMatrix &Mass, utopia::USparseMatrix &L_Mass);

    utopia::USparseMatrix mass_m, mass_lumped;

   



};

#endif /* FractureAppConforming */
