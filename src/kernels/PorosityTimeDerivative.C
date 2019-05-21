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
  return _poro[_qp] * TimeDerivative::computeQpResidual();
}

Real
PorosityTimeDerivative::computeQpJacobian()
{
  return _poro[_qp] * TimeDerivative::computeQpJacobian();
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

void
PorosityTimeDerivative::computeJacobian()
{
    Real volume=0.0;
    Real meanPoro=0.0;
    
    // We need these just to verify
    Real minPoro=1000000.0;
    Real maxPoro=-1000000.0;
    Real volM=(*_current_elem).volume();
    DenseMatrix<Real> testKE;
    testKE.resize(_test.size(),_test.size());
    testKE.zero();
    
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    {
        volume+=_JxW[_qp];
        meanPoro+=_poro[_qp];
        
        // We need these just to verify
        minPoro=std::min(minPoro,_poro[_qp]);
        maxPoro=std::max(maxPoro,_poro[_qp]);
    }

     // We need these just to verify
    if ( fabs(volM-volume)>1e-10 )
    {
        std::cout<<"The volume is different from the sum of nodes\n";
        exit(1);
    }

    if ( fabs(maxPoro-minPoro)>1e-10 )
    {
        std::cout<<"Poro is not constant over the element\n";
        exit(1);
    }
    // up to here
    
    meanPoro/=_qrule->n_points();
    Real entry=meanPoro/_dt*volume/_test.size();

    if (_lumping)
    {
        prepareMatrixTag(_assembly, _var.number(), _var.number());

        precalculateJacobian();
        for (_i = 0; _i < _test.size(); _i++)
            _local_ke(_i, _i) += entry;
        
        for (_i = 0; _i < _test.size(); _i++)
            for (_j = 0; _j < _phi.size(); _j++)
                for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                    testKE(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
        
        testKE-=_local_ke;
        
        if ( std::fabs(testKE.l1_norm() )>1e-10 )
        {
            std::cout<<"There is a difference between algebraic lumping and trapezoidal quadrature\n";
            std::cout<<"Most likely the map between element "<<_current_elem[0].id()<< " and the reference element\n";
            std::cout<<"is non-linear\n";
            std::cout<<"Once we used to exit(1) at this point, nowadays we continue\n";
            //exit(1);
        }

        accumulateTaggedLocalMatrix();
    }
    else
        TimeKernel::computeJacobian();
}
