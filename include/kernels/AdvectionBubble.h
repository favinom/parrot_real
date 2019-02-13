//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADVECTIONBUBBLE_H
#define ADVECTIONBUBBLE_H

#include "Kernel.h"

// Forward Declaration
class AdvectionBubble;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
template <>
InputParameters validParams<AdvectionBubble>();

class AdvectionBubble : public Kernel
{
public:
    AdvectionBubble(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;
    
    virtual Real computeQpJacobian() override {return 0.0;}

    virtual void computeResidual() override;

    virtual void computeJacobian() override;
    
    const MaterialProperty<RealVectorValue> &_U;
    
    const MaterialProperty<Real> &_poro;
    

};

#endif
