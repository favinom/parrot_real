//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SIDEFLUX_H
#define SIDEFLUX_H

//#include "HydraulicConductivity.h"

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
class SideFlux;

template <>
InputParameters validParams<SideFlux>();

/**
 * This postprocessor computes a side integral of the mass flux.
 */
class SideFlux: public SideIntegralVariablePostprocessor
{
public:
  SideFlux(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const MaterialProperty<RealTensorValue> &_K;

};

#endif // SIDEFLUXINTEGRAL_H
