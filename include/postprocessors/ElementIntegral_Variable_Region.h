//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ElementIntegral_variable_region_H
#define ElementIntegral_variable_region_H

#include "ElementIntegralPostprocessor.h"

class ElementIntegral_Variable_Region;

template <>
InputParameters validParams<ElementIntegral_Variable_Region>();

class ElementIntegral_Variable_Region : public ElementIntegralPostprocessor
{
public:
  ElementIntegral_Variable_Region(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const VariableValue & _u;
  
    bool _hasVariable;
    
  const MaterialProperty<int> & _regionID;
  int _region;
    
};

#endif
