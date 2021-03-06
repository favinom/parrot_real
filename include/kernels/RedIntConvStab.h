//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef RedIntConvStab_H
#define RedIntConvStab_H

#include "Kernel.h"

// Forward Declaration
class RedIntConvStab;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
template <>
InputParameters validParams<RedIntConvStab>();

class RedIntConvStab : public Kernel
{
public:
    RedIntConvStab(const InputParameters & parameters);
    
    void myAssembleJacobian();
    
    void myAssembleConvection();
    
    void myAssembleArtificialDiffusion();
    
    void verifyMatrix(DenseMatrix<Number> & in);
    
protected:
    virtual Real computeQpResidual() override {return 0.0;};
    
    virtual Real computeQpJacobian() override {return 0.0;}

    virtual void computeResidual() override;

    virtual void computeJacobian() override;
    
    const VariableValue & _u_old;
    
    const MaterialProperty<RealVectorValue> &_U;
    
    const VariableValue & _u_nodal;
    
    RealVectorValue _vel;
    
    DenseMatrix<Number> _convX;
    DenseMatrix<Number> _convY;
    DenseMatrix<Number> _convZ;
    
    DenseMatrix<Number> _artifDiffX;
    DenseMatrix<Number> _artifDiffY;
    DenseMatrix<Number> _artifDiffZ;
};

#endif
