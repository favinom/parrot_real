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
_phi(declareProperty<Real>("Porosity"))
{}

void
Porosity::computeQpProperties()
{
    
     //std::cout<< _q_point[_qp] <<std::endl; 
    _phi[_qp]=_phiInput;

}
