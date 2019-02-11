/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "CondactivityAux.h"
registerMooseObject("parrot_realApp", CondactivityAux);


template <>
InputParameters
validParams<CondactivityAux>()
{
  // inherit the parameters of AuxKernel:
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned>("comp_i", "component");
  // params.addRequiredParam<unsigned>("comp_j", "component");
  // params.addRequiredParam<bool>("fracture", "fracture");

  
  
  return params;
}

CondactivityAux::CondactivityAux(const InputParameters &parameters): 

    // inherit the parameters of AuxKernel:
   AuxKernel(parameters),
   _comp_i(getParam<unsigned>("comp_i")),
   // _comp_j(getParam<unsigned>("comp_j")),
   // _fracture(getParam<bool>("fracture")),
   _U(getMaterialProperty<RealTensorValue>("VelocityVector"))
{
}

Real
CondactivityAux::computeValue()
{

    // RealTensorValue R1, R2, R3, V;
    // Real theta = 0.5404;
    // RealVectorValue angles(0.0, theta, 0.0);

    // R1(0,0)=std::cos(angles(0));
    // R1(0,1)=-std::sin(angles(0));
    // R1(0,2)=0.0;
    // R1(1,0)=std::sin(angles(0));
    // R1(1,1)=std::cos(angles(0));
    // R1(1,2)=0.0;
    // R1(2,0)=0.0;
    // R1(2,1)=0.0;
    // R1(2,2)=1.0;
    
    // R2(0,0)=std::cos(angles(1));
    // R2(0,1)=0.0;
    // R2(0,2)=-1.0*std::sin(angles(1));
    // R2(1,0)=0.0;
    // R2(1,1)=1.0;
    // R2(1,2)=0.0;
    // R2(2,0)=std::sin(angles(1));
    // R2(2,1)=0.0;
    // R2(2,2)=std::cos(angles(1));
    
    // R3(0,0)=1.0;
    // R3(0,1)=0.0;
    // R3(0,2)=0.0;
    // R3(1,0)=0.0;
    // R3(1,1)=std::cos(angles(2));
    // R3(1,2)=-std::sin(angles(2));
    // R3(2,0)=0.0;
    // R3(2,1)=std::sin(angles(2));
    // R3(2,2)=std::cos(angles(2));

    // V=R1*R2*R3;

    // RealTensorValue _K_inv=V.transpose()*_K[_qp]*V;

    // if(_fracture){
    //   return _K_inv(_comp_i,_comp_j);
    // }
    // else{
    //   return _K[_qp](_comp_i,_comp_j);
    // }

    return _U(_comp_i);
}
