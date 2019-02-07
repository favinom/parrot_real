//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Advection.h"
#include "Function.h"

registerMooseObject("parrot_realApp", Advection);

template <>
InputParameters
validParams<Advection>()
{
    InputParameters params = validParams<Kernel>();
    params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                               "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
    return params;
}

Advection::Advection(const InputParameters & parameters)
: Kernel(parameters),
_U(getMaterialProperty<RealVectorValue>("VelocityVector"))
{}

Real
Advection::computeQpResidual()
{

    
    return _u[_qp]*(_U[_qp] * _grad_test[_i][_qp]);
}


Real
Advection::computeQpJacobian()
{
    
    return _phi[_j][_qp]*(_U[_qp] * _grad_test[_i][_qp]);
}




