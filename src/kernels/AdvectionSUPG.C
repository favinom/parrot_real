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
    params.addRequiredParam<RealVectorValue>("velocity", "Velocity vector");
    params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                               "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
//    params.addRequiredCoupledVar("p",
//                                 "The gradient of this variable will be used as "
//                                 "the velocity vector.");
    params.addRequiredParam<bool>("supg","false","SUPG");
    params.addRequiredParam<RealVectorValue>("velocity", "Velocity vector");
    return params;
}

AdvectionSUPG::AdvectionSUPG(const InputParameters & parameters)
: Kernel(parameters),
_supg(getParam<bool>("supg")),
_velocity(getParam<RealVectorValue>("velocity")),
//_gradP(coupledGradient("p")),
_K(getMaterialProperty<RealTensorValue>("conductivityTensor"))


{
}

Real
AdvectionSUPG::negSpeedQp() const
{
   // std::cout<<"pressure"<< _gradP[_qp] <<std::endl;
    
    //RealVectorValue _velocity = - 1.0  * _epsilon * _K[_qp] * _gradP[_qp];
    
    return - 1.0 * _grad_test[_i][_qp] * _velocity;
}

Real
AdvectionSUPG::computeQpResidual()
{
    //RealVectorValue _velocity = - 1.0 * _epsilon * _K[_qp] * _gradP[_qp];
    
    if (_supg){
        
        Real v_mod = _velocity.norm();
        
        Real h = _current_elem->hmax();
        
        Real coef = 1./(2.0 * v_mod) * h;
        
        return negSpeedQp() * _u[_qp] + coef * _velocity * _grad_test[_i][_qp] * _velocity * _grad_u[_qp];
        
    }
    else{
        
        return negSpeedQp() * _u[_qp];
    }
    
    //-1.0 *  _k * _vel[_qp] * _grad_u[_qp] * _test[_i][_qp]  + coef * vel * _grad_test[_i][_qp] * vel * _grad_u[_qp] ; no integration by parts
    
    //negSpeedQp() * _u[_qp] + coef * vel * _grad_test[_i][_qp] * vel * _grad_u[_qp] ; with integration by parts
}

Real
AdvectionSUPG::computeQpJacobian()
{

    //RealVectorValue _velocity = - 1.0 * _epsilon * _K[_qp] * _gradP[_qp];
 
    if (_supg){
        Real v_mod = _velocity.norm();
        
        Real h = _current_elem->hmax();
        
        Real coef = 1./(2.0 * v_mod) * h;
        
        return negSpeedQp() * _phi[_j][_qp] + coef * _velocity * _grad_test[_i][_qp] * _velocity * _grad_phi[_j][_qp];
    }
    else{
        
        return negSpeedQp() * _phi[_j][_qp];
    }
    // -1.0 *  _k * _vel[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp] + coef * vel * _grad_test[_i][_qp] * vel * _grad_phi[_j][_qp];
    
    //negSpeedQp() * _phi[_j][_qp] + coef * vel * _grad_test[_i][_qp] * vel * _grad_phi[_j][_qp];
}




