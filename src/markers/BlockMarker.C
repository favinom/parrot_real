//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlockMarker.h"

registerMooseObject("parrot_realApp", BlockMarker);

template <>
InputParameters
validParams<BlockMarker>()
{
  InputParameters params = validParams<Marker>();
params.addRequiredParam<int>("blockNum", "which block do you mark?");
  return params;
}

BlockMarker::BlockMarker(const InputParameters & parameters) :
Marker(parameters),
_block(getParam<int>("blockNum"))
{
}

Marker::MarkerValue
BlockMarker::computeElementMarker()
{
    if ((*_current_elem).subdomain_id() == _block)
        return REFINE;
    else
        return DO_NOTHING;
}
