//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ElementIntegral_phi_c_MatProp_H
#define ElementIntegral_phi_c_MatProp_H

#include "ElementIntegralPostprocessor.h"

class ElementIntegral_phi_c_MatProp;

template <>
InputParameters validParams<ElementIntegral_phi_c_MatProp>();

class ElementIntegral_phi_c_MatProp : public ElementIntegralPostprocessor
{
public:
  ElementIntegral_phi_c_MatProp(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  const VariableValue & _u;
  
  const MaterialProperty<Real> & _porosity;
    
  const MaterialProperty<Real> & _scalar;
};

#endif
