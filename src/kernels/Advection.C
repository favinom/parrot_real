//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Advection.h"
#include "Function.h"

registerMooseObject("parrot_realApp", Advection);

template <>
InputParameters
validParams<Advection>()
{
    InputParameters params = validParams<Kernel>();
    params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                               "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
    params.addRequiredParam<bool>("int_by_parts", "true if you want to integrate by parts");
    return params;
}

Advection::Advection(const InputParameters & parameters)
: Kernel(parameters),
_U(getMaterialProperty<RealVectorValue>("VelocityVector")),
_int_by_parts(getParam<bool>("int_by_parts"))
{}

Real
Advection::computeQpResidual()
{
    if(_int_by_parts) 
        return - 1.0 * _grad_test[_i][_qp] * _U[_qp]  * _u[_qp];
    else
         return 1.0 * _grad_u[_qp] * ( _U[_qp] * _test[_i][_qp] );

}


Real
Advection::computeQpJacobian()
{
    if(_int_by_parts) 
        return - 1.0 * _grad_test[_i][_qp] * _U[_qp] * _phi[_j][_qp];
    else
        return 1.0 * _grad_phi[_j][_qp] * ( _U[_qp] * _test[_i][_qp] );
}


void
Advection::computeJacobian()
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
    
//    std::cout<< "Advection==>" << _local_ke <<std::endl<<std::endl;

    ke += _local_ke;
    
    
}



