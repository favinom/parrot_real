#ifndef POROSITY_H
#define POROSITY_H

#include "Material.h"

class Porosity;

template <>
InputParameters validParams<Porosity>();

/**
 * Material for supplying a conductivity tensor.
 * \f$G=\sum_{i=f,s,n}\sigma_i(\vec{e}_i\otimes\vec{e}_i)\f$
 * Conductivities in fibre, sheet and sheet-normal directions are read from the
 * input file.
 * For sensible values, see e.g. \ref Potse2006 "Potse, 2006, Table I".
 */
class Porosity : public Material
{
public:
  Porosity(const InputParameters &parameters);

protected:
  virtual void computeQpProperties();
    
    protected:
    
  Real _phiInput;
  MaterialProperty<Real> &_poro;

};

#endif
