#ifndef HYDRAULICCONDUCTIVITY_H
#define HYDRAULICCONDUCTIVITY_H

#include "Material.h"

class HydraulicConductivity;

template <>
InputParameters validParams<HydraulicConductivity>();

/**
 * Material for supplying a conductivity tensor.
 * \f$G=\sum_{i=f,s,n}\sigma_i(\vec{e}_i\otimes\vec{e}_i)\f$
 * Conductivities in fibre, sheet and sheet-normal directions are read from the
 * input file.
 * For sensible values, see e.g. \ref Potse2006 "Potse, 2006, Table I".
 */
class HydraulicConductivity : public Material
{
public:
  HydraulicConductivity(const InputParameters &parameters);

protected:
  virtual void computeQpProperties();

protected:



  //const MaterialProperty<RealTensorValue> &_fXf, &_sXs, &_nXn;


  //std::vector<Real> _intraConductivities;
  //std::vector<Real> _extraConductivities;
  //std::vector<Real> _monoConductivities;
    
    Real _cond;

  //MaterialProperty<RealTensorValue> &_Gmono;
  MaterialProperty<RealTensorValue> &_K;
  const VariableGradient &_gradP;
  MaterialProperty<RealVectorValue> &_U;
  // MaterialProperty<Real> &_Hsupg;

};

#endif
