/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
/****************************************************************/
/* Barna Becsek, 2017, Uni Bern                                 */
/****************************************************************/

#ifndef STORETRANSFEROPERATORS_H
#define STORETRANSFEROPERATORS_H

// MOOSE includes
#include "GeneralUserObject.h"

// utopia includes
#include "utopia.hpp"
#include <memory>
// Forward Declarations
class StoreTransferOperators;

template<>
InputParameters validParams<StoreTransferOperators>();

class StoreTransferOperators :
  public GeneralUserObject
{
public:
  typedef utopia::DSMatrixd SparseMatT;
  typedef utopia::DVectord VecT;

  // constructor
  StoreTransferOperators(const InputParameters & parameters);

  // returns pointer to stored transfer operators
  std::shared_ptr<SparseMatT> &
  getMatrix()
  {
    return _B;
  };
  std::shared_ptr<SparseMatT> &
  getTransferOperatorReverse()
  {
    return _B_reverse;
  };
    
  std::shared_ptr<VecT> &
  getRHS()
  {
        return _rhs;
  };
    
  std::shared_ptr<SparseMatT> &
  setMatrix()
  {
    return _B;
  };
  std::shared_ptr<SparseMatT> &
  setTransferOperatorReverse()
  {
    return _B_reverse;
  };
    
  std::shared_ptr<VecT> &
  setRHS()
  {
        return _rhs;
  };

  std::shared_ptr<VecT> &
  setOldSol()
  {
        return _sol;
  };

    std::shared_ptr<VecT> &
  getOldSol()
  {
        return _sol;
  };
    
  std::shared_ptr<void> & getVoidPointer()
    {
        return void_ptr;
    }
    
    std::shared_ptr<void> void_ptr;

  void execute(){};
  void initialize(){};
  void finalize(){};

protected:
  std::shared_ptr<SparseMatT> _B;
  std::shared_ptr<SparseMatT> _B_reverse;
  std::shared_ptr<VecT> _rhs, _sol;
};

#endif
