//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Void.h"
#include "NonlinearSystem.h"

registerMooseObject("MooseApp", Void);

template <>
InputParameters
validParams<Void>()
{
  InputParameters params = validParams<TimeIntegrator>();

  return params;
}

Void::Void(const InputParameters & parameters) : TimeIntegrator(parameters) {}

Void::~Void() {}

void
Void::computeTimeDerivatives()
{

//   std::cout<<"I am here A"<<std::endl;
}

void
Void::postResidual(NumericVector<Number> & residual)
{

}
