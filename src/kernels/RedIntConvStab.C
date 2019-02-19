//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RedIntConvStab.h"

#include "libmesh/quadrature_trap.h"

registerMooseObject("parrot_realApp", RedIntConvStab);

template <>
InputParameters
validParams<RedIntConvStab>()
{
    InputParameters params = validParams<Kernel>();
    params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                               "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
    //params.addRequiredParam<bool>("int_by_parts", "true if you want to integrate by parts");
    return params;
}

RedIntConvStab::RedIntConvStab(const InputParameters & parameters) :
Kernel(parameters),
_u_old(_var.slnOld()),
_U(getMaterialProperty<RealVectorValue>("VelocityVector")),
_u_nodal(_var.dofValues())
{}

void
RedIntConvStab::computeResidual()
{
    _vel=RealVectorValue(0.0,0.0,0.0);
    
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _vel+=_U[_qp];
    
    _vel/=_qrule->n_points();
    
    DenseMatrix<Number> convX,convY,convZ;
    convX.resize(_test.size(), _test.size());
    convY.resize(_test.size(), _test.size());
    convZ.resize(_test.size(), _test.size());
    convX.zero();
    convY.zero();
    convZ.zero();
    
    myAssembleJacobian(convX,convY,convZ);
    
    DenseMatrix<Number> artifDiffX;
    DenseMatrix<Number> artifDiffY;
    DenseMatrix<Number> artifDiffZ;
    artifDiffX.resize(_test.size(), _test.size());
    artifDiffY.resize(_test.size(), _test.size());
    artifDiffZ.resize(_test.size(), _test.size());
    artifDiffX.zero();
    artifDiffY.zero();
    artifDiffZ.zero();
    
    myComputeArtificialDiffusion(convX,artifDiffX);
    myComputeArtificialDiffusion(convY,artifDiffY);
    myComputeArtificialDiffusion(convZ,artifDiffZ);
    
    prepareMatrixTag(_assembly, _var.number(), _var.number());
    precalculateJacobian();
    
    _local_ke.resize(_test.size(), _test.size());
    _local_ke.zero();
    _local_ke+=convX;
    _local_ke+=convY;
    _local_ke+=convZ;
    _local_ke+=artifDiffX;
    _local_ke+=artifDiffY;
    _local_ke+=artifDiffZ;
    
    prepareVectorTag(_assembly, _var.number());
    precalculateResidual();
    
    for (_i=0; _i<_test.size(); ++_i)
    {
        for (_j=0; _j<_phi.size(); ++_j)
            _local_re(_i)+=_local_ke(_i,_j)*_u_nodal[_j];
    }

    accumulateTaggedLocalResidual();
    
}

void
RedIntConvStab::computeJacobian()
{
    _vel=RealVectorValue(0.0,0.0,0.0);
    
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _vel+=_U[_qp];
    
    _vel/=_qrule->n_points();

    DenseMatrix<Number> convX,convY,convZ;
    convX.resize(_test.size(), _test.size());
    convY.resize(_test.size(), _test.size());
    convZ.resize(_test.size(), _test.size());
    convX.zero();
    convY.zero();
    convZ.zero();
    
    myAssembleJacobian(convX,convY,convZ);
    
    DenseMatrix<Number> artifDiffX;
    DenseMatrix<Number> artifDiffY;
    DenseMatrix<Number> artifDiffZ;
    artifDiffX.resize(_test.size(), _test.size());
    artifDiffY.resize(_test.size(), _test.size());
    artifDiffZ.resize(_test.size(), _test.size());
    artifDiffX.zero();
    artifDiffY.zero();
    artifDiffZ.zero();

    myComputeArtificialDiffusion(convX,artifDiffX);
    convX+=artifDiffX;
    myComputeArtificialDiffusion(convY,artifDiffY);
    convY+=artifDiffY;
    myComputeArtificialDiffusion(convZ,artifDiffZ);
    convZ+=artifDiffZ;
    
    prepareMatrixTag(_assembly, _var.number(), _var.number());
    precalculateJacobian();

    _local_ke+=convX;
    _local_ke+=convY;
    _local_ke+=convZ;
//    std::cout<<convX<<std::endl;
//        std::cout<<convY<<std::endl;
//        std::cout<<convZ<<std::endl;
//    std::cout<<_local_ke<<std::endl;
//    exit(1);
    
    accumulateTaggedLocalMatrix();
    
}

void RedIntConvStab::myAssembleJacobian(
                                        DenseMatrix<Number> & inX,
                                        DenseMatrix<Number> & inY,
                                        DenseMatrix<Number> & inZ)
{
    UniquePtr<FEBase> fe (FEBase::build(3, FIRST));
    QTrap qrule (3);
    fe->attach_quadrature_rule (&qrule);
    fe->reinit(_current_elem);
    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    const std::vector<std::vector<RealVectorValue> > & dphi = fe->get_dphi();
    
//    DenseMatrix<Number> Me;
//    Me.resize(phi.size(),phi.size());
//    Me.zero();
    
    Real test;
    RealVectorValue gradPhi;
    Real weight;
    for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _test.size(); _j++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            {
                //test=_test[_i][_qp];
                //gradPhi=_grad_test[_j][_qp];
                //weight=_JxW[_qp];
                test=phi[_i][_qp];
                gradPhi=dphi[_j][_qp];
                weight=JxW[_qp];

                inX(_i, _j)+= weight *  gradPhi(0) * test * _vel(0);
                inY(_i, _j)+= weight *  gradPhi(1) * test * _vel(1);
                inZ(_i, _j)+= weight *  gradPhi(2) * test * _vel(2);
                //_coord[_qp] *
            }
}


void RedIntConvStab::myComputeArtificialDiffusion(DenseMatrix<Number> const & op, DenseMatrix<Number> & diff)
{
    for (_i = 0; _i < _test.size(); _i++)
    {
        Real sum = 0.0;
        for (_j = 0; _j < _test.size(); _j++)
        {
            diff(_i,_j)=0.0;
            if (_i!=_j)
            {
                if (op(_i,_j)>0)
                {
                    sum+=op(_i,_j);
                    diff(_i,_j)=-1.0*op(_i,_j);
                }
            }
        }
        diff(_i,_i)+=sum;
    }
}
