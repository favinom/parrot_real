//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef POROSITYTIMEDERIVATIVE_H
#define POROSITYTIMEDERIVATIVE_H

#include "TimeDerivative.h"

// Forward Declarations
class PorosityTimeDerivative;

/**
 * validParams returns the parameters that this Kernel accepts / needs
 * The actual body of the function MUST be in the .C file.
 */
template <>
InputParameters validParams<PorosityTimeDerivative>();

class PorosityTimeDerivative : public TimeDerivative
{
public:
  PorosityTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
//    virtual void computeJacobian();
    virtual void computeResidual();
  /**
   * This MooseArray will hold the reference we need to our
   * material property from the Material class
   */
  const MaterialProperty<Real> &_poro;
    
    const VariableValue & _u_dot_nodal;
    
};
#endif // COEFTIMEDERIVATIVE_H
