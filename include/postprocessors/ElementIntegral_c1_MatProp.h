//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ElementIntegral_c1_MatProp_H
#define ElementIntegral_c1_MatProp_H

#include "ElementIntegralPostprocessor.h"

class ElementIntegral_c1_MatProp;

template <>
InputParameters validParams<ElementIntegral_c1_MatProp>();

class ElementIntegral_c1_MatProp : public ElementIntegralPostprocessor
{
public:
  ElementIntegral_c1_MatProp(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const VariableValue & _u;
  
  const MaterialProperty<Real> & _level_set_1;
    
};

#endif
