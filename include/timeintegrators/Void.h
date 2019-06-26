//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef VOID_H
#define VOID_H

#include "TimeIntegrator.h"

class Void;

template <>
InputParameters validParams<Void>();

/**
 * Implicit Euler's method
 */
class Void : public TimeIntegrator
{
public:
  Void(const InputParameters & parameters);
  virtual ~Void();

  virtual int order() override { return 1; }
  virtual void computeTimeDerivatives() override;
  virtual void postResidual(NumericVector<Number> & residual) override;

protected:
};

#endif /* IMPLICITEULER_H */
