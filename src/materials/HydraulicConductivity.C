#include "HydraulicConductivity.h"

registerMooseObject("parrot_realApp", HydraulicConductivity);

template <>
InputParameters
validParams<HydraulicConductivity>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("conductivity","conductivity");
  return params;
}

HydraulicConductivity::HydraulicConductivity(const InputParameters &parameters) :
Material(parameters),
_cond(getParam<Real>("conductivity")),
_K(declareProperty<RealTensorValue>("conductivityTensor"))
{}

void
HydraulicConductivity::computeQpProperties()
{


_K[_qp]=RealTensorValue(1.0,0.0,0.0,
                        0.0,1.0,0.0,
                        0.0,0.0,1.0);

    _K[_qp]=_cond*_K[_qp];

}
