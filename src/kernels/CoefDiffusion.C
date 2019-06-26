//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoefDiffusion.h"

registerMooseObject("parrot_realApp", CoefDiffusion);

template <>
InputParameters
validParams<CoefDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addCustomTypeParam("coef", 0.0, "CoefficientType", "The coefficient of diffusion");
  params.addParam<Real>("non_controllable", "A parameter we cannot control.");

  params.declareControllable("coef");

  return params;
}

CoefDiffusion::CoefDiffusion(const InputParameters & parameters)
  : Kernel(parameters), _coef(getParam<Real>("coef")), _poro(getMaterialProperty<Real>("Porosity"))
{
}

Real
CoefDiffusion::computeQpResidual()
{
     Real stab = _coef * _poro[_qp] * ( _current_elem->hmin() * _current_elem->hmin())/(6.0 *_dt);
  return stab * _grad_test[_i][_qp] * _grad_u[_qp];
}

Real
CoefDiffusion::computeQpJacobian()
{
   Real stab = _coef * _poro[_qp] * ( _current_elem->hmin() * _current_elem->hmin())/(6.0 * _dt);
  return stab * _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}
