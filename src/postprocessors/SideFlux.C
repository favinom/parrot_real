//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SideFlux.h"
//#include "HydraulicConductivity.h"

registerMooseObject("parrot_realApp", SideFlux);

template <>
InputParameters
validParams<SideFlux>()
{
  InputParameters params = validParams<SideIntegralVariablePostprocessor>();
  params.addClassDescription("Computes the integral of the flux over the specified boundary");
  return params;
}

SideFlux::SideFlux(const InputParameters & parameters)
  : SideIntegralVariablePostprocessor(parameters),
   _K(getMaterialProperty<RealTensorValue>("conductivityTensor"))
{
}

Real
SideFlux::computeQpIntegral()
{
  return 1.0 * _K[_qp] * _grad_u[_qp] * _normals[_qp];
}
