//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
// preapred by Maria Nestola

#include "AdvectionUpwind.h"

registerMooseObject("parrot_realApp", AdvectionUpwind);

template <>
InputParameters
validParams<AdvectionUpwind>()
{
    InputParameters params = validParams<Kernel>();
    params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                               "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");
    MooseEnum upwinding_type("none full quick", "none");
    params.addParam<MooseEnum>("upwinding_type",
                               upwinding_type,
                               "Type of upwinding used.  None: Typically results in overshoots and "
                               "undershoots, but numerical diffusion is minimized.  Full: Overshoots "
                               "and undershoots are avoided, but numerical diffusion is large");
    return params;
}

AdvectionUpwind::AdvectionUpwind(const InputParameters & parameters)
: Kernel(parameters),
_upwinding(getParam<MooseEnum>("upwinding_type").getEnum<UpwindingType>()),
_gradP(coupledGradient("p")),
_K(getMaterialProperty<RealTensorValue>("conductivityTensor")),
_u_nodal(_var.dofValues()),
_upwind_node(0),
_dtotal_mass_out(0)
{
}



Real
AdvectionUpwind::computeQpResidual()
{
    RealVectorValue _velocity = - _K[_qp] * _gradP[_qp];
    
    return (- 1.0 * _grad_test[_i][_qp] * _velocity) * _u[_qp];
}

Real
AdvectionUpwind::computeQpJacobian()
{
    
    RealVectorValue _velocity = - _K[_qp] * _gradP[_qp];
    
    return (- 1.0 * _grad_test[_i][_qp] * _velocity) * _phi[_j][_qp];
}

void
AdvectionUpwind::computeResidual()
{
    switch (_upwinding)
    {
        case UpwindingType::none:
            Kernel::computeResidual();
            break;
        case UpwindingType::full:
            fullUpwind(JacRes::CALCULATE_RESIDUAL);
            break;
    }
}

void
AdvectionUpwind::computeJacobian()
{
    switch (_upwinding)
    {
        case UpwindingType::none:
            Kernel::computeJacobian();
            break;
        case UpwindingType::full:
            fullUpwind(JacRes::CALCULATE_JACOBIAN);
            break;
    }
}


void
AdvectionUpwind::fullUpwind(JacRes res_or_jac)
{

    RealVectorValue _velocity = - 1.0 * _K[_qp] * _gradP[_qp];
    
    const unsigned int num_nodes = _test.size();

    prepareVectorTag(_assembly, _var.number());
    
    if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
        prepareMatrixTag(_assembly, _var.number(), _var.number());
    
    _upwind_node.resize(num_nodes);
    for (_i = 0; _i < num_nodes; ++_i)
    {
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
            _local_re(_i) += _JxW[_qp] * _coord[_qp] * (- 1.0 * _grad_test[_i][_qp] * _velocity);
        _upwind_node[_i] = (_local_re(_i) >= 0.0);
    }
    
    Real total_mass_out = 0.0;
    Real total_in = 0.0;
    if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
        _dtotal_mass_out.assign(num_nodes, 0.0);
    
    for (unsigned int n = 0; n < num_nodes; ++n)
    {
        if (_upwind_node[n])
        {
            if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
            {
                if (_test.size() == _phi.size())
                    
                    _local_ke(n, n) += _local_re(n);
                
                _dtotal_mass_out[n] += _local_ke(n, n);
            }
            
            _local_re(n) *= _u_nodal[n];
            
            total_mass_out += _local_re(n);
        }
        else
            total_in -= _local_re(n);
    }
    
    for (unsigned int n = 0; n < num_nodes; ++n)
    {
        if (!_upwind_node[n])
        {
            if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
                for (_j = 0; _j < _phi.size(); _j++)
                    _local_ke(n, _j) += _local_re(n) * _dtotal_mass_out[_j] / total_in;
            _local_re(n) *= total_mass_out / total_in;
        }
    }
    
    if (res_or_jac == JacRes::CALCULATE_RESIDUAL)
    {
        accumulateTaggedLocalResidual();
        
        if (_has_save_in)
        {
            Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
            for (const auto & var : _save_in)
                var->sys().solution().add_vector(_local_re, var->dofIndices());
        }
    }
    
    if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    {
        accumulateTaggedLocalMatrix();
        
        if (_has_diag_save_in)
        {
            unsigned int rows = _local_ke.m();
            DenseVector<Number> diag(rows);
            for (unsigned int i = 0; i < rows; i++)
                diag(i) = _local_ke(i, i);
            
            Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
            for (const auto & var : _diag_save_in)
                var->sys().solution().add_vector(diag, var->dofIndices());
        }
    }
}
