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
    params.addRequiredParam<Real>("coef","stab coef");
    params.addRequiredParam<bool>("use_h","use h size");
    return params;
}

AdvectionSUPG::AdvectionSUPG(const InputParameters & parameters)
: Kernel(parameters),
_coef(getParam<Real>("coef")),
_U(getMaterialProperty<RealVectorValue>("VelocityVector")),
_use_h(getParam<bool>("use_h"))


{
}



Real
AdvectionSUPG::computeQpResidual()
{
    
    
    Real v_mod = _U[_qp].norm();

    //std::cout<<"param_res"<< v_mod <<std::endl;

    Real stab =0.0;

    if (_use_h)  {

        Real h = _current_elem->hmax();

        stab = 1.0 * _coef * h * 1./(v_mod);

    }
    else {

         stab = _coef;

    }


    return 1.0 * stab * _U[_qp] * _grad_test[_i][_qp] * _U[_qp] * _grad_u[_qp];
}

Real
AdvectionSUPG::computeQpJacobian()
{

    Real v_mod = _U[_qp].norm();


    Real stab =0.0;

    if (_use_h)  {

        Real h = _current_elem->hmax();

        stab = 1.0 * _coef * h * 1./(v_mod);

        //std::cout<<"param_jac"<< v_mod <<std::endl;

    }
    else {

         stab = _coef;

    }
    
    return stab * _U[_qp] * _grad_test[_i][_qp] * _U[_qp] * _grad_phi[_j][_qp];
}


void
AdvectionSUPG::computeJacobian()
{

    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
    _local_ke.resize(ke.m(), ke.n());
    _local_ke.zero();

    
    precalculateJacobian();

    for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < _phi.size(); _j++)
            for (_qp = 0; _qp < _qrule->n_points(); _qp++){
                _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
                //std::cout<<"qpoints = "<<_qrule->n_points()<<std::endl;
            }
    
    //std::cout<<"stabilized==>"<<_local_ke <<std::endl<<std::endl;

    ke += _local_ke;
    
    
}




