//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef PARROT_REALAPP_H
#define PARROT_REALAPP_H

#include "MooseApp.h"

class parrot_realApp;

template <>
InputParameters validParams<parrot_realApp>();

class parrot_realApp : public MooseApp
{
public:
  parrot_realApp(InputParameters parameters);
  virtual ~parrot_realApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void registerObjectDepends(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
  static void associateSyntaxDepends(Syntax & syntax, ActionFactory & action_factory);
  static void registerExecFlags(Factory & factory);
};

#endif /* PARROT_REALAPP_H */
