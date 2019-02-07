#ifndef HYDRAULICCONDUCTIVITYFRACTURE_H
#define HYDRAULICCONDUCTIVITYFRACTURE_H

#include "Material.h"

class HydraulicConductivityFracture;

template <>
InputParameters validParams<HydraulicConductivityFracture>();

/**
 * Material for supplying a conductivity tensor.
 * \f$G=\sum_{i=f,s,n}\sigma_i(\vec{e}_i\otimes\vec{e}_i)\f$
 * Conductivities in fibre, sheet and sheet-normal directions are read from the
 * input file.
 * For sensible values, see e.g. \ref Potse2006 "Potse, 2006, Table I".
 */
class HydraulicConductivityFracture : public Material
{
public:
    HydraulicConductivityFracture(const InputParameters &parameters);
    
protected:
    virtual void computeQpProperties();
    
protected:
    
    
    
    //const MaterialProperty<RealTensorValue> &_fXf, &_sXs, &_nXn;
    
    std::vector<Real> _cond;
    std::vector<Real> _theta;
    
    MaterialProperty<RealTensorValue> &_K;
    const VariableGradient &_gradP;
    MaterialProperty<RealVectorValue> &_U;
    
};

#endif
