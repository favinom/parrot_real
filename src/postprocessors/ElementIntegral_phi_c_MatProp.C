//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegral_phi_c_MatProp.h"

registerMooseObject("parrot_realApp", ElementIntegral_phi_c_MatProp);

template <>
InputParameters
validParams<ElementIntegral_phi_c_MatProp>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addRequiredParam<MaterialPropertyName>("mat_prop", "The name of the material property");
  return params;
}

ElementIntegral_phi_c_MatProp::ElementIntegral_phi_c_MatProp(const InputParameters & parameters) :
ElementIntegralPostprocessor(parameters),
_u(coupledValue("variable")),
_porosity(getMaterialProperty<Real>("Porosity")),
_scalar(getMaterialProperty<Real>("mat_prop"))
{
}

Real
ElementIntegral_phi_c_MatProp::computeQpIntegral()
{
  return _u[_qp]*_scalar[_qp]*_porosity[_qp];
}
