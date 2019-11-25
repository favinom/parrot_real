#include "PorosityFunction.h"

// MOOSE includes
//#include "Assembly.h"
//#include "MooseMesh.h"
#include "Function.h"
registerMooseObject("parrot_realApp", PorosityFunction);

template <>
InputParameters
validParams<PorosityFunction>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("phi","phi");
  params.addParam<FunctionName>("function",
                                  "The function that describes the pressure");

  return params;
}

PorosityFunction::PorosityFunction(const InputParameters &parameters) :
Material(parameters),
_phiInput(getParam<Real>("phi")),
_phi(declareProperty<Real>("Porosity")),
_function(getFunction("function"))
{}

void
PorosityFunction::computeQpProperties()
{
    _phi[_qp]=_phiInput * _function.value(_t, _q_point[_qp]);
}
