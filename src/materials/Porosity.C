#include "Porosity.h"

// MOOSE includes
//#include "Assembly.h"
//#include "MooseMesh.h"

registerMooseObject("parrot_realApp", Porosity);

template <>
InputParameters
validParams<Porosity>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("phi","phi");

  return params;
}

Porosity::Porosity(const InputParameters &parameters) :
Material(parameters),
_phiInput(getParam<Real>("phi")),
_poro(declareProperty<Real>("Porosity"))
{}

void
Porosity::computeQpProperties()
{
    _poro[_qp]=_phiInput;


     //std::cout<<"poro_material naterial ==>"<<_poro[_qp]<<std::endl;

  

}
