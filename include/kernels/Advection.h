//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADVECTION_H
#define ADVECTION_H

#include "Kernel.h"

// Forward Declaration
class Advection;

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
template <>
InputParameters validParams<Advection>();

class Advection : public Kernel
{
public:
    Advection(const InputParameters & parameters);
    
protected:
    virtual Real computeQpResidual() override;
    
    virtual Real computeQpJacobian() override;

    virtual void computeJacobian() override;
    
    const MaterialProperty<RealVectorValue> &_U;

    bool _int_by_parts;
    
    

};

#endif

brick x 0.49995 y 0.49995 z 0.0001
brick x 0.49995 y 0.49995 z 0.0001
brick x 0.49995 y 0.49995 z 0.0001
brick x 0.49995 y 0.49995 z 0.0001
move Volume 8 9 10 11 x 0.250025 y 0.250025
rotate Volume 9 angle 90  about Z 
rotate Volume 10 angle 180  about Z 
rotate Volume 11 angle 270 about Z
brick x 0.49995 y 0.0001 z 0.49995 
brick x 0.49995 y 0.0001 z 0.49995 
brick x 0.49995 y 0.0001 z 0.49995 
brick x 0.49995 y 0.0001 z 0.49995 
move Volume 12 13 14 15 x 0.250025 z 0.250025
rotate Volume 13 angle 90  about Y 
rotate Volume 14 angle 180  about Y 
rotate Volume 15 angle 270 about Y
brick x 0.0001 y 0.49995 z 0.49995 
brick x 0.0001 y 0.49995 z 0.49995 
brick x 0.0001 y 0.49995 z 0.49995 
brick x 0.0001 y 0.49995 z 0.49995 
move Volume 16 17 18 19 y 0.250025 z 0.250025
rotate Volume 17 angle 90  about X 
rotate Volume 18 angle 180  about X 
rotate Volume 19 angle 270 about X
merge volume all
unite volume all
Volume 1 copy
Volume 20 scale 0.5
move Volume 20 x 0.25 y 0.25 z 0.25

Volume 20 copy
Volume 20 copy
Volume 20 copy
Volume 20 copy
Volume 20 copy
Volume 20 copy
Volume 20 copy
Volume 20 copy
Volume 20 copy
mv 20 x -0.5
move Volume 22 x 0.5
move Volume 22 x -0.5
move Volume 22 x -0.5
move Volume 20 y -0.5
move Volume 20 y -0.5
move Volume 20 y 0.5
move Volume 23 y 0.5
move Volume 23 y -0.5
move Volume 23 y -0.5
move Volume 23 x -0.5
move Volume 24 z -0.5
move Volume 25 z -0.5
move Volume 25 y -0.5
move Volume 26 x -0.5
move Volume 26 x 0.5
move Volume 26 z -0.5
move Volume 27 z -0.5
move Volume 27 x -0.5
move Volume 21 z -0.5
move Volume 21 y -0.5
move Volume 21 x -0.5
