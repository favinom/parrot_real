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

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*       MECH - ICS Mechanical simulation framework             */
/*                Prepared by Maria Nestola,                     */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/* Auxiliary Kernel to visualize the fibers.                    */
/****************************************************************/


#ifndef PorosityAux_H
#define PorosityAux_H

#include "AuxKernel.h"

class PorosityAux;

template <>
InputParameters validParams<PorosityAux>();

class PorosityAux : public AuxKernel
{
public:
  PorosityAux(const InputParameters &parameters);

protected:
  virtual Real computeValue();


  // unsigned int _comp_j;
  // bool _fracture;

  MaterialProperty<Real> const &_poro;
};
#endif
