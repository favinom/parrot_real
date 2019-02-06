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
    //    params.addRequiredParam<RealVectorValue>("velocity", "Velocity vector");
    params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                               "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
    params.addRequiredParam<Real>("epsilon", "epsilon");
    params.addRequiredCoupledVar("p",
                                 "The gradient of this variable will be used as "
                                 "the velocity vector.");
    return params;
}

AdvectionSUPG::AdvectionSUPG(const InputParameters & parameters)
: Kernel(parameters),
_epsilon(getParam<Real>("epsilon")),
_gradP(coupledGradient("p")),
_K(getMaterialProperty<RealTensorValue>("conductivityTensor"))


{
}

Real
AdvectionSUPG::negSpeedQp() const
{
    std::cout<<"pressure"<< _gradP[_qp] <<std::endl;
    
    RealVectorValue _velocity = - 1.0  * _epsilon * _K[_qp] * _gradP[_qp];
    
    return - 1.0 * _grad_test[_i][_qp] * _velocity;
}

Real
AdvectionSUPG::computeQpResidual()
{
    RealVectorValue _velocity = - 1.0 * _epsilon * _K[_qp] * _gradP[_qp];
    
    Real v_mod = _velocity.norm();
    //std::sqrt(_velocity(0) * _velocity(0) + _velocity(1) * _velocity(1) + _velocity(2) * _velocity(2));
    
    Real h = _current_elem->hmax();
    
    Real coef = 1./(2.0 * v_mod) * h;
    
    return negSpeedQp() * _u[_qp] + coef * _velocity * _grad_test[_i][_qp] * _velocity * _grad_u[_qp];
    
    //-1.0 *  _k * _vel[_qp] * _grad_u[_qp] * _test[_i][_qp]  + coef * vel * _grad_test[_i][_qp] * vel * _grad_u[_qp] ; no integration by parts
    
    //negSpeedQp() * _u[_qp] + coef * vel * _grad_test[_i][_qp] * vel * _grad_u[_qp] ; with integration by parts
}

Real
AdvectionSUPG::computeQpJacobian()
{

    RealVectorValue _velocity = - 1.0 * _epsilon * _K[_qp] * _gradP[_qp];
    
    Real v_mod = _velocity.norm();
    
    //Real v_mod = std::sqrt(_velocity(0) * _velocity(0) + _velocity(1) * _velocity(1) + _velocity(2) * _velocity(2));
    
    Real h = _current_elem->hmax();
    
    Real coef = 1./(2.0 * v_mod) * h;
    
    return negSpeedQp() * _phi[_j][_qp] + coef * _velocity * _grad_test[_i][_qp] * _velocity * _grad_phi[_j][_qp];

    // -1.0 *  _k * _vel[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp] + coef * vel * _grad_test[_i][_qp] * vel * _grad_phi[_j][_qp];
    
    //negSpeedQp() * _phi[_j][_qp] + coef * vel * _grad_test[_i][_qp] * vel * _grad_phi[_j][_qp];
}




