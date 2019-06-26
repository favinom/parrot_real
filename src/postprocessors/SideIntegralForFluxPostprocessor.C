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

#include "SideIntegralForFluxPostprocessor.h"
#include "MooseMesh.h"


registerMooseObject("parrot_realApp",SideIntegralForFluxPostprocessor);

template <>
InputParameters
validParams<SideIntegralForFluxPostprocessor>()
{
  InputParameters params = validParams<SideIntegralVariablePostprocessor>();
  return params;
}

SideIntegralForFluxPostprocessor::SideIntegralForFluxPostprocessor(
    const InputParameters &parameters)
    : SideIntegralVariablePostprocessor(parameters)
{

}

Real
SideIntegralForFluxPostprocessor::computeQpIntegral()
{

  Real J0 = -1.0 * _grad_u[_qp]  * _normals[_qp];


  return J0;
}
