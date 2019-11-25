#ifndef POROSITYFUNCTION_H
#define POROSITYFUNCTION_H

#include "Material.h"

class PorosityFunction;

template <>
InputParameters validParams<PorosityFunction>();

/**
 * Material for supplying a conductivity tensor.
 * \f$G=\sum_{i=f,s,n}\sigma_i(\vec{e}_i\otimes\vec{e}_i)\f$
 * Conductivities in fibre, sheet and sheet-normal directions are read from the
 * input file.
 * For sensible values, see e.g. \ref Potse2006 "Potse, 2006, Table I".
 */
class PorosityFunction : public Material
{
public:
  PorosityFunction(const InputParameters &parameters);

protected:
  virtual void computeQpProperties();
    
    protected:
    
  Real _phiInput;
  MaterialProperty<Real> &_phi;
  Function &_function;

};

#endif
