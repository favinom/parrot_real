//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "parrot_realTestApp.h"
#include "parrot_realApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<parrot_realTestApp>()
{
  InputParameters params = validParams<parrot_realApp>();
  return params;
}

parrot_realTestApp::parrot_realTestApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  parrot_realApp::registerObjectDepends(_factory);
  parrot_realApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  parrot_realApp::associateSyntaxDepends(_syntax, _action_factory);
  parrot_realApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  parrot_realApp::registerExecFlags(_factory);

  bool use_test_objs = getParam<bool>("allow_test_objects");
  if (use_test_objs)
  {
    parrot_realTestApp::registerObjects(_factory);
    parrot_realTestApp::associateSyntax(_syntax, _action_factory);
    parrot_realTestApp::registerExecFlags(_factory);
  }
}

parrot_realTestApp::~parrot_realTestApp() {}

void
parrot_realTestApp::registerApps()
{
  registerApp(parrot_realApp);
  registerApp(parrot_realTestApp);
}

void
parrot_realTestApp::registerObjects(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new test objects here! */
}

void
parrot_realTestApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
  /* Uncomment Syntax and ActionFactory parameters and register your new test objects here! */
}

void
parrot_realTestApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execute flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
parrot_realTestApp__registerApps()
{
  parrot_realTestApp::registerApps();
}

// External entry point for dynamic object registration
extern "C" void
parrot_realTestApp__registerObjects(Factory & factory)
{
  parrot_realTestApp::registerObjects(factory);
}

// External entry point for dynamic syntax association
extern "C" void
parrot_realTestApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  parrot_realTestApp::associateSyntax(syntax, action_factory);
}

// External entry point for dynamic execute flag loading
extern "C" void
parrot_realTestApp__registerExecFlags(Factory & factory)
{
  parrot_realTestApp::registerExecFlags(factory);
}
