//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AdvectionSUPG.h"
#include "Function.h"

registerMooseObject("parrot_realApp", AdvectionSUPG);

template <>
InputParameters
validParams<AdvectionSUPG>()
{
    InputParameters params = validParams<Kernel>();
    // params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
    //                            "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
    params.addRequiredParam<Real>("coef","stab coef");
    params.addRequiredParam<bool>("use_h","use h size");
    // params.addRequiredCoupledVar("p",
    //                              "The gradient of this variable will be used as "
    //                              "the velocity vector.");
    return params;
}

AdvectionSUPG::AdvectionSUPG(const InputParameters & parameters)
: Kernel(parameters),
_coef(getParam<Real>("coef")),
// _gradP(coupledGradient("p")),
_U(getMaterialProperty<RealVectorValue>("VelocityVector")),
_use_h(getParam<bool>("use_h"))
// _Hsupg(getMaterialProperty<RealValue>("Hsupg"))


{
}



Real
AdvectionSUPG::computeQpResidual()
{
    
    
    Real v_mod = _U[_qp].norm();

    Real stab =0.0;

    if (_use_h)  {

        Real h = _current_elem->hmax();

        stab = _coef * v_mod * h;

    }
    else {

         stab = _coef * v_mod;

    }


    return stab * _U[_qp] * _grad_test[_i][_qp] * _U[_qp] * _grad_u[_qp];
}


void
AdvectionSUPG::computeJacobian()
{
    
   
    
    Real v_mod = _U[_qp] .norm();


    Real stab =0.0;

    if (_use_h)  {

        Real h = _current_elem->hmax();

        stab = _coef * v_mod * h;
    }

    else {

         stab = _coef * v_mod;

    }


    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
    _local_ke.resize(ke.m(), ke.n());
    _local_ke.zero();

    
    precalculateJacobian();
    for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _phi.size(); _j++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++){
                _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * stab * _U[_qp] * _grad_test[_i][_qp] * _U[_qp] * _grad_phi[_j][_qp];
                //std::cout<<"qpoints = "<<_qrule->n_points()<<std::endl;
            }
    


    ke += _local_ke;
    
    
}




