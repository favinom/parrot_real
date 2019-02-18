//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorosityTimeDerivative.h"

registerMooseObject("parrot_realApp", PorosityTimeDerivative);

template <>
InputParameters
validParams<PorosityTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}

PorosityTimeDerivative::PorosityTimeDerivative(const InputParameters & parameters) :
TimeDerivative(parameters),
_poro(getMaterialProperty<Real>("Porosity")),
_u_dot_nodal( _var.dofValuesDot() )
{}

Real
PorosityTimeDerivative::computeQpResidual()
{
  // We're reusing the TimeDerivative Kernel's residual
  // so that we don't have to recode that.
  return _poro[_qp] * _u_dot[_qp] * _test[_i][_qp];
}

Real
PorosityTimeDerivative::computeQpJacobian()
{
  return _poro[_qp] * _du_dot_du[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}

void
PorosityTimeDerivative::computeResidual()
{
    if (_lumping)
    {
        // These have to be called from mother class
        prepareVectorTag(_assembly, _var.number());
        precalculateResidual();
        
        Real volume=0.0;
        Real meanPoro=0.0;
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        {
            volume+=_JxW[_qp];
            meanPoro+=_poro[_qp];
        }
        
        meanPoro/=_qrule->n_points();
        Real entry=meanPoro*volume/_test.size();
        
        prepareMatrixTag(_assembly, _var.number(), _var.number());
        
        precalculateJacobian();
        for (_i = 0; _i < _test.size(); _i++)
            _local_re(_i) += entry*_u_dot_nodal[_i];
        
        accumulateTaggedLocalResidual();
    }
    else
        TimeKernel::computeResidual();
}
