//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AdvectionBubble.h"
#include "Function.h"

registerMooseObject("parrot_realApp", AdvectionBubble);

template <>
InputParameters
validParams<AdvectionBubble>()
{
    InputParameters params = validParams<Kernel>();
    params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                               "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
    //params.addRequiredParam<bool>("int_by_parts", "true if you want to integrate by parts");
    return params;
}

AdvectionBubble::AdvectionBubble(const InputParameters & parameters) :
Kernel(parameters),
 _u_old(_var.slnOld()),
_U(getMaterialProperty<RealVectorValue>("VelocityVector")),
_poro(getMaterialProperty<Real>("Porosity")),
_u_dot(_var.uDot())
//_int_by_parts(getParam<bool>("int_by_parts"))
{}

Real
AdvectionBubble::computeQpResidual()
{

        return 1.0 * _grad_u[_qp] * ( _U[_qp] * _test[_i][_qp] );
    
}


void
AdvectionBubble::computeResidual()
{
    DenseVector<Number> my_re;
    my_re.resize(_test.size()+1);
    my_re.zero();
    
    DenseVector<Number> my_re_b;
    my_re_b.resize(_test.size()+1);
    my_re.zero();
    
    prepareVectorTag(_assembly, _var.number());
    
    precalculateResidual();
    
    for (_i = 0; _i < _test.size()+1; _i++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        {
            Real test;
            if (_i!=_test.size())
                test=_test[_i][_qp];
            else
            {
                test=_test[0][_qp]*_test[1][_qp];
            }
            
            my_re(_i)+=_JxW[_qp] * _coord[_qp] *
            (_poro[_qp]*(_u[_qp]-_u_old[_qp])/_dt+_U[_qp]*_grad_u[_qp])*test;
            //my_re(_i)+=_JxW[_qp] * _coord[_qp] * (_poro[_qp]*_u_old[_qp]/_dt)*test;
        }

    std::vector<Real> bubble;
    std::vector<RealVectorValue> bubbleGrad;
    bubble.resize(_qrule->n_points());
    bubbleGrad.resize(_qrule->n_points());
    
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    {
        bubble[_qp]=_test[0][_qp]*_test[1][_qp];
        bubbleGrad[_qp]=_grad_test[0][_qp]*_test[1][_qp]+
        _test[0][_qp]*_grad_test[1][_qp];
    }
    
    for (_i = 0; _i < _test.size()+1; _i++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        {
            Real test;
            if (_i!=_test.size())
                test=_test[_i][_qp];
            else
            {
                test=_test[0][_qp]*_test[1][_qp];
            }
            
            
            
            
            my_re_b(_i)+=_JxW[_qp] * _coord[_qp] * (_poro[_qp]/_dt*bubble[_qp]+_U[_qp]*bubbleGrad[_qp])*test;
            //my_re_b(_i)+=_JxW[_qp] * _coord[_qp] * (_poro[_qp]/_dt*bubble[_qp])*test;
        }
    
    for (_i = 0; _i < _test.size(); _i++)
        _local_re(_i) = my_re(_i)-my_re_b(_i)/my_re_b(_test.size())*my_re(_test.size());

    
    
//    for (_i = 0; _i < _test.size(); _i++)
//        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//            _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();
    
    accumulateTaggedLocalResidual();
    
}

void
AdvectionBubble::computeJacobian()
{
    prepareMatrixTag(_assembly, _var.number(), _var.number());
    
    DenseMatrix<Number> my_ke;
    my_ke.resize(_test.size()+1, _phi.size()+1);
    
    precalculateJacobian();
    
    for (_i = 0; _i < _test.size()+1; _i++)
        for (_j = 0; _j < _phi.size()+1; _j++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            {
                Real phi, test;
                RealVectorValue gradPhi;
                if (_i!=_test.size())
                    test=_test[_i][_qp];
                else
                {
                    test=_test[0][_qp]*_test[1][_qp];
                    
                    //                for (int f=0; f<_test.size(); ++f)
                    //                    test=test*_test[f][_qp];
                }
                
                if (_j!=_phi.size())
                {
                    phi=_phi[_j][_qp];
                    gradPhi=_grad_phi[_j][_qp];
                }
                else
                {
                    phi=_phi[0][_qp]*_phi[1][_qp];
                    gradPhi=_grad_phi[0][_qp]*_phi[1][_qp]+_phi[0][_qp]*_grad_phi[1][_qp];
                    //                    for (int f=0; f<_test.size(); ++f)
                    //                        phi=phi*_phi[f][_qp];
                }
                
                
                
                my_ke(_i, _j)+=_JxW[_qp] * _coord[_qp] *
                (_poro[_qp]/_dt*test * phi+
                (gradPhi * _U[_qp]) * test) ;
                
                
                //if (_i!=_test.size() && _j!=_phi.size())
                //    _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * _test[_i][_qp] * _phi[_j][_qp];
            }
    
    for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _phi.size(); _j++)
            _local_ke(_i, _j)=my_ke(_i, _j)
            -my_ke(_test.size(), _j )*my_ke( _i , _phi.size())/my_ke(_test.size(), _phi.size());
    
    accumulateTaggedLocalMatrix();
    
}

