#include "HydraulicConductivityFracture.h"

registerMooseObject("parrot_realApp", HydraulicConductivityFracture);

template <>
InputParameters
validParams<HydraulicConductivityFracture>()
{
    InputParameters params = validParams<Material>();
    params.addRequiredParam<std::vector<Real>>("conductivity","conductivity");
    params.addRequiredParam<std::vector<Real>>("theta","theta");
    params.addCoupledVar("pressure",
                                 "The gradient of this variable will be used as "
                                 "the velocity vector.");
    return params;
}

HydraulicConductivityFracture::HydraulicConductivityFracture
(const InputParameters &parameters) :
Material(parameters),
_cond(getParam<std::vector<Real>>("conductivity")),
_theta(getParam<std::vector<Real>>("theta")),
_K(declareProperty<RealTensorValue>("conductivityTensor")),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"):_grad_zero),
_U(declareProperty<RealVectorValue>("VelocityVector"))
{}

void
HydraulicConductivityFracture::computeQpProperties()
{
    
    Real pi=4.0*std::atan(1.0);
    
    _K[_qp]=RealTensorValue(1.0,0.0,0.0,
                            0.0,1.0,0.0,
                            0.0,0.0,1.0);
    _K[_qp](0,0)=_cond.at(0);
    _K[_qp](1,1)=_cond.at(1);
    _K[_qp](2,2)=_cond.at(2);
    
    RealVectorValue angles;
    for (int i=0; i<3; ++i)
        angles(i)=pi*_theta.at(i)/180.0;
    
    RealTensorValue R1;
    RealTensorValue R2;
    RealTensorValue R3;
    
    
    
    R1(0,0)=std::cos(angles(0));
    R1(0,1)=-std::sin(angles(0));
    R1(0,2)=0.0;
    R1(1,0)=std::sin(angles(0));
    R1(1,1)=std::cos(angles(0));
    R1(1,2)=0.0;
    R1(2,0)=0.0;
    R1(2,1)=0.0;
    R1(2,2)=1.0;
    
    R2(0,0)=std::cos(angles(1));
    R2(0,1)=0.0;
    R2(0,2)=-std::sin(angles(1));
    R2(1,0)=0.0;
    R2(1,1)=1.0;
    R2(1,2)=0.0;
    R2(2,0)=std::sin(angles(1));
    R2(2,1)=0.0;
    R2(2,2)=std::cos(angles(1));
    
    R3(0,0)=1.0;
    R3(0,1)=0.0;
    R3(0,2)=0.0;
    R3(1,0)=0.0;
    R3(1,1)=std::cos(angles(2));
    R3(1,2)=-std::sin(angles(2));
    R3(2,0)=0.0;
    R3(2,1)=std::sin(angles(2));
    R3(2,2)=std::cos(angles(2));
    
    RealTensorValue V=R1*R2*R3;
    
    _K[_qp]=V*_K[_qp]*V.transpose();

    _U[_qp] = -1.0 * _K[_qp] * _gradP[_qp];
}
