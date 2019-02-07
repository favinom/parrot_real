#include "HydraulicConductivity.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseMesh.h"

registerMooseObject("parrot_realApp", HydraulicConductivity);

template <>
InputParameters
validParams<HydraulicConductivity>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("conductivity","conductivity");
  params.addCoupledVar("pressure",
                                 "The gradient of this variable will be used as "
                                 "the velocity vector.");
  return params;
}

HydraulicConductivity::HydraulicConductivity(const InputParameters &parameters) :
Material(parameters),
_cond(getParam<Real>("conductivity")),
_K(declareProperty<RealTensorValue>("conductivityTensor")),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"): _grad_zero),
_U(declareProperty<RealVectorValue>("VelocityVector"))
// _Hsupg(declareProperty<Real>("Hsupg"))
{}

void
HydraulicConductivity::computeQpProperties()
{


_K[_qp]=RealTensorValue(1.0,0.0,0.0,
                        0.0,1.0,0.0,
                        0.0,0.0,1.0);

    _K[_qp]=_cond*_K[_qp];
    
    _U[_qp] =  -1.0 * _K[_qp] * _gradP[_qp];

    // Real v_mod = _U[_qp] .norm();

    // _Hsupg[_qp] = _current_elem->hmax() * v_mod;

}
