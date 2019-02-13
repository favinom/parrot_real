//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AlgebraicDiffusion.h"

registerMooseObject("parrot_realApp", AlgebraicDiffusion);

template <>
InputParameters
validParams<AlgebraicDiffusion>()
{
    InputParameters params = validParams<Kernel>();
    params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                               "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
    //params.addRequiredParam<bool>("int_by_parts", "true if you want to integrate by parts");
    return params;
}

AlgebraicDiffusion::AlgebraicDiffusion(const InputParameters & parameters) :
Kernel(parameters),
 _u_old(_var.slnOld()),
_U(getMaterialProperty<RealVectorValue>("VelocityVector")),
_u_nodal(_var.dofValues()),
_poro(getMaterialProperty<Real>("Porosity"))
//_int_by_parts(getParam<bool>("int_by_parts"))
{}

Real
AlgebraicDiffusion::computeQpResidual()
{

        return 1.0 * _grad_u[_qp] * ( _U[_qp] * _test[_i][_qp] );
    
}


void
AlgebraicDiffusion::computeResidual()
{
    prepareVectorTag(_assembly, _var.number());
    
    precalculateResidual();
    
    for (_i = 0; _i < _test.size(); _i++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        {
            Real test;

                test=_test[_i][_qp];
            
            _local_re(_i)+=_JxW[_qp] * _coord[_qp] *
            (_poro[_qp]*(_u[_qp]-_u_old[_qp])/_dt+_U[_qp]*_grad_u[_qp])*test;
        }
    

    
    accumulateTaggedLocalResidual();
    
}

void
AlgebraicDiffusion::computeJacobian()
{
    DenseMatrix<Number> artifDiff;
    artifDiff.resize(_test.size(), _phi.size());
    artifDiff.zero();
    
    prepareMatrixTag(_assembly, _var.number(), _var.number());
    
    precalculateJacobian();
    
    for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _phi.size(); _j++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            {
                Real phi, test;
                RealVectorValue gradPhi;
                    test=_test[_i][_qp];
                    phi=_phi[_j][_qp];
                    gradPhi=_grad_phi[_j][_qp];
                
                
                _local_ke(_i, _j)+=_JxW[_qp] * _coord[_qp] *
                (_poro[_qp]/_dt*test * phi+
                (gradPhi * _U[_qp]) * test) ;
                
            }
    
    for (_i = 0; _i < _test.size(); _i++)
    {
        Real sum =0.0;
        for (_j = 0; _j < _phi.size(); _j++)
        {
            if (_i!=_j)
            {
                if (_local_ke(_i,_j)>0.0)
                {
                    sum+=_local_ke(_i,_j);
                    artifDiff(_i,_j)=-_local_ke(_i,_j);
                }
            }
        }
        artifDiff(_i,_i)+=sum;
    }
//    std::cout<<_local_ke<<std::endl;
    _local_ke+=artifDiff;
//    std::cout<<_local_ke<<std::endl;
    
    accumulateTaggedLocalMatrix();
    
    prepareVectorTag(_assembly, _var.number());
    DenseVector<Number> artifResidual;
    artifResidual.resize(_test.size());
    artifResidual.zero();
    
    for (_i=0; _i<_test.size(); ++_i)
    {
        artifResidual(_i)=0;
        for (_j=0; _j<_phi.size(); ++_j)
            artifResidual(_i)+=artifDiff(_i,_j)*_u_nodal[_j];
        _local_re(_i)+=artifResidual(_i);
    }
    
    accumulateTaggedLocalResidual();
    
}
