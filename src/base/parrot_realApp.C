#include "parrot_realApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<parrot_realApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

parrot_realApp::parrot_realApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  parrot_realApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  parrot_realApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  ModulesApp::registerExecFlags(_factory);
  parrot_realApp::registerExecFlags(_factory);
}

parrot_realApp::~parrot_realApp() {}

void
parrot_realApp::registerApps()
{
  registerApp(parrot_realApp);
}

void
parrot_realApp::registerObjects(Factory & factory)
{
    Registry::registerObjectsTo(factory, {"parrot_realApp"});
}

void
parrot_realApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & action_factory)
{
  Registry::registerActionsTo(action_factory, {"parrot_realApp"});

  /* Uncomment Syntax parameter and register your new production objects here! */
}

void
parrot_realApp::registerObjectDepends(Factory & /*factory*/)
{
}

void
parrot_realApp::associateSyntaxDepends(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}

void
parrot_realApp::registerExecFlags(Factory & /*factory*/)
{
  /* Uncomment Factory parameter and register your new execution flags here! */
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
parrot_realApp__registerApps()
{
  parrot_realApp::registerApps();
}

extern "C" void
parrot_realApp__registerObjects(Factory & factory)
{
  parrot_realApp::registerObjects(factory);
}

extern "C" void
parrot_realApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  parrot_realApp::associateSyntax(syntax, action_factory);
}

extern "C" void
parrot_realApp__registerExecFlags(Factory & factory)
{
  parrot_realApp::registerExecFlags(factory);
}
