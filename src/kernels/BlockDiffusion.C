//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlockDiffusion.h"

registerMooseObject("parrot_realApp", BlockDiffusion);

template <>
InputParameters
validParams<BlockDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak " "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addRequiredParam<Real>("coef","coef");
  return params;
}

BlockDiffusion::BlockDiffusion(const InputParameters & parameters) :
Kernel(parameters),
//_K(getMaterialProperty<RealTensorValue>("conductivityTensor")),
_coef(getParam<Real>("coef"))
{}

Real
BlockDiffusion::computeQpResidual()
{
  return _coef * _grad_u[_qp] * (_grad_test[_i][_qp]);
}

Real
BlockDiffusion::computeQpJacobian()
{
  return _coef * _grad_phi[_j][_qp] * (_grad_test[_i][_qp]);
}
