//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegral_c1_MatProp.h"

registerMooseObject("parrot_realApp", ElementIntegral_c1_MatProp);

template <>
InputParameters
validParams<ElementIntegral_c1_MatProp>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  return params;
}

ElementIntegral_c1_MatProp::ElementIntegral_c1_MatProp(const InputParameters & parameters) :
ElementIntegralPostprocessor(parameters),
_u(coupledValue("variable")),
_level_set_1(getMaterialProperty<Real>("level_set_1"))
{
}

Real
ElementIntegral_c1_MatProp::computeQpIntegral()
{
  return _u[_qp]*_level_set_1[_qp];
}
