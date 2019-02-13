//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ALGEBRA_G
#define ALGEBRA_G

#include "Kernel.h"

// Forward Declaration
class AlgebraicDiffusion;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
template <>
InputParameters validParams<AlgebraicDiffusion>();

class AlgebraicDiffusion : public Kernel
{
public:
    AlgebraicDiffusion(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;
    
    virtual Real computeQpJacobian() override {return 0.0;}

    virtual void computeResidual() override;

    virtual void computeJacobian() override;
    
    const VariableValue & _u_old;
    
    const MaterialProperty<RealVectorValue> &_U;
    
    const MaterialProperty<Real> &_poro;
    

};

#endif
