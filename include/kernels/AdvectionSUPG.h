//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADVECTIONSUPG_H
#define ADVECTIONSUPG_H

#include "Kernel.h"

// Forward Declaration
class AdvectionSUPG;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
template <>
InputParameters validParams<AdvectionSUPG>();

class AdvectionSUPG : public Kernel
{
public:
    AdvectionSUPG(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;

    virtual Real computeQpJacobian() override;
    
    virtual void computeJacobian() override;

    Real _coef;
    
    // const VariableGradient &_gradP;
    
    const MaterialProperty<RealVectorValue> &_U;

    bool _use_h;
    // const MaterialProperty<Real> &_Hsupg;


};

#endif // CONSERVATIVEADVECTION_H
