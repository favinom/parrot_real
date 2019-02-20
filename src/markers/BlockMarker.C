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
    params.addRequiredParam<std::vector<int>>("blockNum", "which block do you mark?");
    return params;
}

BlockMarker::BlockMarker(const InputParameters & parameters) :
Marker(parameters),
_block(getParam<std::vector<int>>("blockNum"))
{
}

Marker::MarkerValue
BlockMarker::computeElementMarker()
{
    for (int j=0; j<_block.size(); j++){
        if ((*_current_elem).subdomain_id() == _block[j]){
            return REFINE;
        }
    }
    return DO_NOTHING;

}
