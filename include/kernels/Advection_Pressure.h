//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADVECTION_PRESSURE_HPP
#define ADVECTION_PRESSURE_HPP

#include "Kernel.h"

// Forward Declaration
class Advection_Pressure;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
template <>
InputParameters validParams<Advection_Pressure>();

class Advection_Pressure : public Kernel
{
public:
    Advection_Pressure(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;
    
    virtual Real computeQpJacobian() override;

    //virtual void computeJacobian() override;

    RealVectorValue negSpeedQp() const;

    const MaterialProperty<Real> &_poro;

    const VariableGradient &_vel;

    Real _k;
    
    

};

#endif
