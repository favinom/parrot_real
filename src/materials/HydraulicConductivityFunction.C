#include "HydraulicConductivityFunction.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseMesh.h"
#include "Function.h"
registerMooseObject("parrot_realApp", HydraulicConductivityFunction);

template <>
InputParameters
validParams<HydraulicConductivityFunction>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("conductivity","conductivity");
  params.addCoupledVar("pressure",
                                 "The gradient of this variable will be used as "
                                 "the velocity vector.");
    params.addParam<FunctionName>("function",
                                  "The function that describes the pressure");
  return params;
}

HydraulicConductivityFunction::HydraulicConductivityFunction(const InputParameters &parameters) :
Material(parameters),
_cond(getParam<Real>("conductivity")),
_K(declareProperty<RealTensorValue>("conductivityTensor")),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"): _grad_zero),
_U(declareProperty<RealVectorValue>("VelocityVector")),
_function(getFunction("function"))
// _Hsupg(declareProperty<Real>("Hsupg"))
{}

void
HydraulicConductivityFunction::computeQpProperties()
{


_K[_qp]=RealTensorValue(1.0,0.0,0.0,
                        0.0,1.0,0.0,
                        0.0,0.0,1.0);

    _K[_qp]=_cond *_K[_qp] * _function.value(_t, _q_point[_qp]);
    
    _U[_qp] =  -1.0 * _K[_qp] * _gradP[_qp];

    // Real v_mod = _U[_qp] .norm();

    // _Hsupg[_qp] = _current_elem->hmax() * v_mod;

}
