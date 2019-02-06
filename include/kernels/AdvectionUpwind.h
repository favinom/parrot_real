//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADVECTIONUPWIND_H
#define ADVECTIONUPWIND_H

#include "Kernel.h"

class AdvectionUpwind;

template <>
InputParameters validParams<AdvectionUpwind>();

class AdvectionUpwind : public Kernel
{
public:
    AdvectionUpwind(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual void computeResidual() override;
    virtual void computeJacobian() override;
    

    
    enum class JacRes
    {
        CALCULATE_RESIDUAL = 0,
        CALCULATE_JACOBIAN = 1
    };
    
    const enum class UpwindingType { none, full} _upwinding;
    
    const VariableGradient &_gradP;
    
    const MaterialProperty<RealTensorValue> & _K;
    
    const VariableValue & _u_nodal;
    
    std::vector<bool> _upwind_node;
    
    std::vector<Real> _dtotal_mass_out;
    
    void fullUpwind(JacRes res_or_jac);
};

#endif // CONSERVATIVEADVECTION_H
