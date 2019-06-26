//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef BLOCKMARKER_H
#define BLOCKMARKER_H

#include "Marker.h"

class BlockMarker;

template <>
InputParameters validParams<BlockMarker>();

class BlockMarker : public Marker
{
public:
  BlockMarker(const InputParameters & parameters);

protected:
  virtual MarkerValue computeElementMarker() override;

  std::vector<int> _block;
};

#endif /* UNIFORMMARKER_H */
