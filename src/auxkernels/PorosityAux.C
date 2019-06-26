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

#include "PorosityAux.h"
registerMooseObject("parrot_realApp", PorosityAux);


template <>
InputParameters
validParams<PorosityAux>()
{
  // inherit the parameters of AuxKernel:
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned>("comp_i", "component");


  
  
  return params;
}

PorosityAux::PorosityAux(const InputParameters &parameters): 
   AuxKernel(parameters),
   _poro(getMaterialProperty<Real>("Porosity"))
{
}

Real
PorosityAux::computeValue()
{



    return _poro[_qp];
}
