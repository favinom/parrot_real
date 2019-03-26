//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimeAdvectionSUPG.h"
#include "Function.h"

registerMooseObject("parrot_realApp", TimeAdvectionSUPG);

template <>
InputParameters
validParams<TimeAdvectionSUPG>()
{
    InputParameters params = validParams<TimeDerivative>();
    params.addRequiredParam<Real>("coef","stab coef");
    params.addRequiredParam<bool>("use_h","use h size");
    return params;
}

TimeAdvectionSUPG::TimeAdvectionSUPG(const InputParameters & parameters)
: TimeDerivative(parameters),
_coef(getParam<Real>("coef")),
_U(getMaterialProperty<RealVectorValue>("VelocityVector")),
_use_h(getParam<bool>("use_h")),
_rho(getMaterialProperty<Real>("Porosity"))


{ }


Real
TimeAdvectionSUPG::computeQpResidual()
{
    
    
    Real v_mod = _U[_qp].norm();

    Real stab = 0.0;

       if (_use_h)  
    {

        if(v_mod > 1e-6) {

            Real h = _current_elem->hmax();

            stab = 1.0 * _coef * h * 1./(v_mod);
        }
         

    }

    else 
    {

         stab = _coef;

    }


    return stab * _rho[_qp] * _U[_qp] * _grad_test[_i][_qp] * TimeDerivative::computeQpResidual();
}

Real
TimeAdvectionSUPG::computeQpJacobian()
{

    Real v_mod = _U[_qp].norm();


    Real stab =0.0;

    if (_use_h)  {

        Real h = _current_elem->hmax();

        stab = 1.0 * _coef * h * 1./(v_mod);

    }
    else {

         stab = _coef;

    }
    
    return stab * _rho[_qp] * ( _U[_qp] * _grad_test[_i][_qp] ) * TimeDerivative::computeQpJacobian();
}

