//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorosityTimeDerivative.h"

registerMooseObject("parrot_realApp", PorosityTimeDerivative);

template <>
InputParameters
validParams<PorosityTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}

PorosityTimeDerivative::PorosityTimeDerivative(const InputParameters & parameters) :
TimeDerivative(parameters),
_phi(getMaterialProperty<Real>("Porosity"))
{}

Real
PorosityTimeDerivative::computeQpResidual()
{
  // We're reusing the TimeDerivative Kernel's residual
  // so that we don't have to recode that.
  return _phi[_qp] * TimeDerivative::computeQpResidual();
}

Real
PorosityTimeDerivative::computeQpJacobian()
{
  return _phi[_qp] * TimeDerivative::computeQpJacobian();
}
