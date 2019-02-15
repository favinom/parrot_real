//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "myDiffusion.h"

registerMooseObject("parrot_realApp", MyDiffusion);

template <>
InputParameters
validParams<MyDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak " "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addRequiredParam<Real>("coef","coef");
  return params;
}

MyDiffusion::MyDiffusion(const InputParameters & parameters) :
Kernel(parameters),
_K(getMaterialProperty<RealTensorValue>("conductivityTensor")),
_coef(getParam<Real>("coef"))
{}

Real
MyDiffusion::computeQpResidual()
{
  return _coef * _grad_u[_qp] * (_K[_qp] *_grad_test[_i][_qp]);
}

Real
MyDiffusion::computeQpJacobian()
{
  return _coef * _grad_phi[_j][_qp] * (_K[_qp]* _grad_test[_i][_qp]);
}
