//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PARROTPROBLEM_H
#define PARROTPROBLEM_H

#include "FEProblem.h"


class ParrotProblem;

template <>
InputParameters validParams<ParrotProblem>();

class ParrotProblem : public FEProblem
{
public:
    ParrotProblem(const InputParameters & parameters);
};


#endif
