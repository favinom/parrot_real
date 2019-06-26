/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef SIDEINTEGRALFORFLUXPOSTPROCESSOR_H
#define SIDEINTEGRALFORFLUXPOSTPROCESSOR_H

#include "SideIntegralVariablePostprocessor.h"
//#include "MooseVariableInterface.h"

// Forward Declarations
class SideIntegralForFluxPostprocessor;


template <>
InputParameters validParams<SideIntegralForFluxPostprocessor>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class SideIntegralForFluxPostprocessor : public SideIntegralVariablePostprocessor
{
public:
  SideIntegralForFluxPostprocessor(const InputParameters &parameters);

protected:
  virtual Real computeQpIntegral();



};

#endif
