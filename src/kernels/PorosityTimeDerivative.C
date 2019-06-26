//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorosityTimeDerivative.h"

#include "libmesh/quadrature_trap.h"

#define TOLL 1e-8
#define DIM 3

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
        
        
        myComputeLumpedJacobian();
        
        
        prepareVectorTag(_assembly, _var.number());
        
        precalculateResidual();
        for (_i = 0; _i < _test.size(); _i++)
            _local_re(_i) += _local_ke(_i,_i)*_dt *_u_dot_nodal[_i];
        
        accumulateTaggedLocalResidual();
        
        
        
    }
    else
        TimeKernel::computeResidual();
}

void
PorosityTimeDerivative::computeJacobian()
{
    if (_lumping)
    {
        prepareMatrixTag(_assembly, _var.number(), _var.number());
        precalculateJacobian();
        myComputeLumpedJacobian();
        accumulateTaggedLocalMatrix();
        
    }
    else
        TimeKernel::computeJacobian();
}


void PorosityTimeDerivative::myComputeLumpedJacobian()
{
    Real meanPoro= 0.0;
    Real minPoro = 1.0e15;
    Real maxPoro =-1.0e15;
    
    Real volume0=0.0;
    Real volume1=0.0;
    Real volume2=0.0;
    
    Real volM=(*_current_elem).volume();
    
    
    // Here we compute the sum of the weigths with standard quadrature rule
    // and the max, min, and mean poro
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    {
        meanPoro+=_poro[_qp];
        minPoro=std::min(minPoro,_poro[_qp]);
        maxPoro=std::max(maxPoro,_poro[_qp]);
    }
    
    if ( std::fabs(maxPoro-minPoro)>TOLL && _fe_problem.currentlyComputingJacobian() )
    {
        std::cout<<"Porosity is not constant over the element\n";
        exit(1);
    }
    
    meanPoro/=_qrule->n_points();
    
    if ( std::fabs(meanPoro-minPoro)>TOLL && _fe_problem.currentlyComputingJacobian() )
    {
        std::cout<<"meanPoro is wrong\n";
        exit(1);
    }
    if ( std::fabs(meanPoro-maxPoro)>TOLL && _fe_problem.currentlyComputingJacobian() )
    {
        std::cout<<"meanPoro is wrong\n";
        exit(1);
    }
    
    _Ke0.resize(_test.size(),_test.size());
    _Ke1.resize(_test.size(),_test.size());
    _Ke2.resize(_test.size(),_test.size());
    _Ke0.zero();
    _Ke1.zero();
    _Ke2.zero();
    
    // In _Ke0 we assemble the algebraic lumping,
    // we sum the extra-diagonal entries to the diagonal one
    
    for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _test.size(); _j++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
                _Ke0(_i, _i) += _JxW[_qp] * _coord[_qp] * _poro[_qp]/_dt* _test[_j][_qp]* _test[_i][_qp];
    
    for (int _qp=0; _qp<_qrule->n_points(); _qp++)
        volume0+=_JxW[_qp];
    
    // In _Ke1 we use the trapezoidal rule from libmesh
    
    int dim=DIM;
    UniquePtr<FEBase> fe1 (FEBase::build(dim, FIRST));
    QTrap qrule1 (dim);
    fe1->attach_quadrature_rule (&qrule1);
    fe1->reinit(_current_elem);
    const std::vector<Real>               & JxW1 = fe1->get_JxW();
    const std::vector<std::vector<Real> > & phi1 = fe1->get_phi();
    
    for (_i = 0; _i < phi1.size(); _i++)
        for (_j = 0; _j < phi1.size(); _j++)
            for (_qp = 0; _qp < qrule1.n_points(); _qp++)
            {
                Real test = phi1[_i][_qp];
                Real fhi =  phi1[_j][_qp];
                Real w    = JxW1[_qp];
                
                _Ke1(_i, _j)+= meanPoro/_dt * w *  fhi * test;
            }
    
    for (int _qp=0; _qp<qrule1.n_points(); _qp++)
        volume1+=JxW1[_qp];
    
    // In _Ke2 we use our trapezoidal (in the current configuration)
    
    std::vector<Number> myW(_test.size(),0.0);
    
    for (_i=0; _i<_test.size(); ++_i)
    {
        for (_qp=0; _qp<_qrule->n_points(); _qp++)
        {
            myW.at(_i)+=_JxW[_qp]*_test[_i][_qp];
        }
        
        _Ke2(_i, _i)=myW.at(_i)*meanPoro/_dt;
    }
    
    // Here we compute the sum of the weigths with the trapezoidal quadrature rule
    for (_qp = 0; _qp < myW.size(); _qp++)
    {
        volume2+=myW.at(_qp);
    }
    
    // We need these just to verify
    if ( std::fabs(volM-volume0)>TOLL && _fe_problem.currentlyComputingJacobian() )
    {
        std::cout<<"volume0="<<volume0<<std::endl;
        std::cout<<"volM="<<volM<<std::endl;
        std::cout<<"difference="<<std::fabs(volM-volume0)<<std::endl;
        std::cout<<"The volume0 is different from volM\n";
    }
    
    // We need these just to verify
    if ( std::fabs(volM-volume1)>TOLL && _fe_problem.currentlyComputingJacobian() )
    {
        std::cout<<"The volume1 is different from volM\n";
        exit(1);
    }
    
    // We need these just to verify
    if ( std::fabs(volM-volume2)>TOLL && _fe_problem.currentlyComputingJacobian() )
    {
        std::cout<<"volume2="<<volume2<<std::endl;
        std::cout<<"volM="<<volM<<std::endl;
        std::cout<<"difference="<<std::fabs(volM-volume2)<<std::endl;
        std::cout<<"The volume2 is different from volM\n";
    }
    
    DenseMatrix<Number> diff02(_Ke0);
    diff02-=_Ke2;

    if ( std::fabs(diff02.l1_norm() )>TOLL && _fe_problem.currentlyComputingJacobian() )
    {
        std::cout<<"Algebraic lumping and mass lumping do not coincide, exiting...\n";
        exit(1);
    }
    
    DenseMatrix<Number> diff01(_Ke0);
    diff01-=_Ke1;
    
    if ( std::fabs(diff01.l1_norm() )>TOLL && _fe_problem.currentlyComputingJacobian() )
    {
        std::cout<<"Element "<<_current_elem[0].id()<<std::endl;
        std::cout<<"Difference between mass lumping in the reference and current configuration \n";
        std::cout<<"Most likely the map between the element and the reference element is non-linear.\n";
        std::cout<<"Once at this point we used to exit(1), nowadays we continue\n\n\n\n";
        // exit(1);
        
    }
    
    _local_ke=_Ke0;
    
}

