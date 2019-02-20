//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SideIntegralFunctionVariable.h"

registerMooseObject("parrot_realApp", SideIntegralFunctionVariable);

template <>
InputParameters
validParams<SideIntegralFunctionVariable>()
{
  InputParameters params = validParams<SideIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable",
                               "The name of the variable that this boundary condition applies to");
  return params;
}

SideIntegralFunctionVariable::SideIntegralFunctionVariable
(const InputParameters & parameters) :
SideIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
_u(coupledValue("variable")),
_grad_u(coupledGradient("variable"))
{
  addMooseVariableDependency(mooseVariable());
}

Real
SideIntegralFunctionVariable::computeQpIntegral()
{
    RealVectorValue punto=_q_point[_qp];
    if (1.0/3.0 <= punto(2) && punto(2) <= 2.0/3.0)
        return _u[_qp]/3.0;
    else
        return 0.0;
}
