//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MYDIFFUSION_H
#define MYDIFFUSION_H

#include "Kernel.h"

class MyDiffusion;

template <>
InputParameters validParams<MyDiffusion>();

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class MyDiffusion : public Kernel
{
public:
  MyDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
    
  const MaterialProperty<RealTensorValue> &_K;

  Real _coef;    
};

#endif /* MYDIFFUSION_H */
