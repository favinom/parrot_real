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
    
    _local_ke.resize(_test.size(), _test.size());
    _local_ke.zero();
    myAssembleJacobian();
    
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
    
    prepareMatrixTag(_assembly, _var.number(), _var.number());
    precalculateJacobian();

    myAssembleJacobian();
    
    accumulateTaggedLocalMatrix();
    
}

void RedIntConvStab::myAssembleConvection()
{
    _convX.resize(_test.size(), _test.size());
    _convY.resize(_test.size(), _test.size());
    _convZ.resize(_test.size(), _test.size());
    
    _convX.zero();
    _convY.zero();
    _convZ.zero();
    
    UniquePtr<FEBase> fe (FEBase::build(3, FIRST));
    QTrap qrule (3);
    fe->attach_quadrature_rule (&qrule);
    fe->reinit(_current_elem);
    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    const std::vector<std::vector<RealVectorValue> > & dphi = fe->get_dphi();
    
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

                _convX(_i, _j)+= weight *  gradPhi(0) * test * _vel(0);
                _convY(_i, _j)+= weight *  gradPhi(1) * test * _vel(1);
                _convZ(_i, _j)+= weight *  gradPhi(2) * test * _vel(2);
            }
}


void RedIntConvStab::myAssembleArtificialDiffusion()
{
    UniquePtr<FEBase> fe (FEBase::build(3, FIRST));
    QTrap qrule (3);
    fe->attach_quadrature_rule (&qrule);
    fe->reinit(_current_elem);
    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<RealVectorValue> > & dphi = fe->get_dphi();
    
    
    const std::vector<std::vector<Real> > & phi = fe->get_phi();
    
    DenseMatrix<Number> A(phi.size(),phi.size());
    
    for (_i = 0; _i < phi.size(); _i++)
        for (_j = 0; _j < phi.size(); _j++)
            for (_qp = 0; _qp < qrule.n_points(); _qp++)
            {
                Real test = phi[_i][_qp];
                Real fhi =  phi[_j][_qp];
                Real w    = JxW[_qp];
                
                A(_i, _j)+= w *  fhi * test;
                
            }
    std::cout<<A<<std::endl;

    
    
    _artifDiffX.resize(dphi.size(), dphi.size());
    _artifDiffY.resize(dphi.size(), dphi.size());
    _artifDiffZ.resize(dphi.size(), dphi.size());
    
    _artifDiffX.zero();
    _artifDiffY.zero();
    _artifDiffZ.zero();
    
    RealVectorValue gradPhi,gradTest;
    Real weight;
    for (_i = 0; _i < dphi.size(); _i++)
        for (_j = 0; _j < dphi.size(); _j++)
            for (_qp = 0; _qp < qrule.n_points(); _qp++)
            {
                gradTest= dphi[_i][_qp];
                gradPhi = dphi[_j][_qp];
                weight  = JxW[_qp];
                
                _artifDiffX(_i, _j)+= weight *  gradPhi(0) * gradTest(0);
                _artifDiffY(_i, _j)+= weight *  gradPhi(1) * gradTest(1);
                _artifDiffZ(_i, _j)+= weight *  gradPhi(2) * gradTest(2);
                
            }

}

void RedIntConvStab::verifyMatrix(DenseMatrix<Number> & in)
{
    for (int i=0; i<in.m(); ++i)
        for (int j=0; j<in.n(); ++j)
            if ( i!=j && in(i,j)>1e-15)
            {
                std::cout<<"RedIntConvStab::verifyMatrix error\n";
                std::cout<<_artifDiffZ<<std::endl;
                std::cout<<_convZ<<std::endl;
                exit(1);
            }
}

void RedIntConvStab::myAssembleJacobian()
{
    myAssembleConvection();
    myAssembleArtificialDiffusion();
    
    Real coeffX,coeffY,coeffZ;
    
        coeffX=std::fabs( _convX(0,0)/_artifDiffX(0,0) );
        coeffY=std::fabs( _convY(0,0)/_artifDiffY(0,0) );
        coeffZ=std::fabs( _convZ(0,0)/_artifDiffZ(0,0) );

    if (coeffX<0.0)
    {
        std::cout<<"coeffX <0.0"<<coeffX<<std::endl;
        exit(1);
    }
    if (coeffY<0.0)
    {
        std::cout<<"coeffY <0.0"<<coeffY<<std::endl;
        exit(1);
    }
    if (coeffZ<0.0)
    {
        std::cout<<"coeffZ <0.0"<<coeffZ<<std::endl;
        exit(1);
    }

//    std::cout<<_artifDiffX<<std::endl;
    std::cout<<_convX<<std::endl;
    
    _artifDiffX*=coeffX;
    _artifDiffY*=coeffY;
    _artifDiffZ*=coeffZ;
    
//    verifyMatrix(_artifDiffX);
//    verifyMatrix(_artifDiffY);
//    verifyMatrix(_artifDiffZ);

    _convX+=_artifDiffX;
    _convY+=_artifDiffY;
    _convZ+=_artifDiffZ;
    
//    verifyMatrix(_convX);
//    verifyMatrix(_convY);
//   verifyMatrix(_convZ);
    
    _local_ke+=_convX;
    _local_ke+=_convY;
    _local_ke+=_convZ;
    
//    verifyMatrix(_local_ke);
    
}
