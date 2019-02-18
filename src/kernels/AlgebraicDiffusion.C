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
            (_U[_qp]*_grad_u[_qp])*test;
        }
    
    // PARTE DEL RESIDUO LEGATA A DIFFUSIONE ARTIFICIALE
    
    DenseVector<Number> artifResidual;
    artifResidual.resize(_test.size());
    artifResidual.zero();
    
    DenseMatrix<Number> myOperator;
    myOperator.resize(_test.size(), _test.size());
    myOperator.zero();

    DenseMatrix<Number> artifDiff;
    artifDiff.resize(_test.size(), _test.size());
    artifDiff.zero();

    myAssembleJacobian(myOperator);
    //std::cout<<myOperator<<std::endl;
    myComputeArtificialDiffusion(myOperator,artifDiff);

    for (_i=0; _i<_test.size(); ++_i)
    {
        artifResidual(_i)=0;
        for (_j=0; _j<_phi.size(); ++_j)
            artifResidual(_i)+=artifDiff(_i,_j)*_u_nodal[_j];
        _local_re(_i)+=artifResidual(_i);
    }

    accumulateTaggedLocalResidual();
    
}

void
AlgebraicDiffusion::computeJacobian()
{
    DenseMatrix<Number> myOperator;
    myOperator.resize(_test.size(), _phi.size());
    myOperator.zero();

    DenseMatrix<Number> artifDiff;
    artifDiff.resize(_test.size(), _phi.size());
    artifDiff.zero();
    
    myAssembleJacobian(myOperator);
    myComputeArtificialDiffusion(myOperator,artifDiff);

    prepareMatrixTag(_assembly, _var.number(), _var.number());
    precalculateJacobian();
    
    _local_ke+=myOperator;
    _local_ke+=artifDiff;
    
    accumulateTaggedLocalMatrix();
    
}

void AlgebraicDiffusion::myAssembleJacobian(DenseMatrix<Number> & in)
{
    for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _test.size(); _j++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            {
                Real phi, test;
                RealVectorValue gradPhi;
                test=_test[_i][_qp];
                phi=_test[_j][_qp];
                gradPhi=_grad_test[_j][_qp];
                
                
                in(_i, _j)+=_JxW[_qp] * _coord[_qp] *
                ((gradPhi * _U[_qp]) * test) ;
                
            }

}

void AlgebraicDiffusion::myComputeArtificialDiffusion(DenseMatrix<Number> const & op, DenseMatrix<Number> & diff)
{
    for (_i = 0; _i < _test.size(); _i++)
    {
        Real sum =0.0;
        for (_j = 0; _j < _test.size(); _j++)
        {
            diff(_i,_j)=0.0;
            if (_i!=_j)
            {
                if (op(_i,_j)>0)
                {
                    //double max_1=std::max(op(_i,_j),op(_j,_i));
                    //sum+=-1.0*std::max(0.0,max_1);
                    sum+=op(_i,_j);
                    diff(_i,_j)=-1*op(_i,_j);
                }
            }
        }
        diff(_i,_i)+=sum;
    }
}
