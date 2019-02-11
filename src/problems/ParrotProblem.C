//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ParrotProblem.h"

registerMooseObject("parrot_realApp", ParrotProblem);

template <>
InputParameters
validParams<ParrotProblem>()
{
  InputParameters params = validParams<FEProblem>();
  return params;
}

ParrotProblem::ParrotProblem(const InputParameters & parameters) :
FEProblem(parameters)
{
    _const_jacobian=true;
}

void ParrotProblem::initialSetup()
{
    FEProblem::initialSetup();
    
    Moose::PetscSupport::petscSetOptions(*this);
    
    PetscErrorCode ierr;
    
    PetscNonlinearSolver<Number> * petsc_solver =
    dynamic_cast<PetscNonlinearSolver<Number> *>((*_nl_sys).nonlinearSolver());

    SNES snes = petsc_solver->snes();
    KSP ksp;
    SNESType ttype;
    KSPType type;
    ierr = SNESGetKSP(snes,&ksp);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESGetType(snes, &ttype);
    std::cout<< "SNES type from PARROTPROBLEM exec:"<<ttype<<std::endl;
    ierr = KSPGetType(ksp, &type);
    std::cout<< "KSP  type from PARROTPROBLEM exec:"<<type<<std::endl;
    _factorized = 0 ;
};

void ParrotProblem::timestepSetup()
{
    FEProblem::timestepSetup();
    
    PetscErrorCode ierr;
    Moose::PetscSupport::petscSetOptions(*this);

    PetscNonlinearSolver<Number> * petsc_solver =
    dynamic_cast<PetscNonlinearSolver<Number> *>((*_nl_sys).nonlinearSolver());
    SNES snes = petsc_solver->snes();
    KSP ksp;
    SNESType ttype;
    KSPType type;

    ierr = SNESGetKSP(snes,&ksp);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESGetType(snes, &ttype);
    ierr = KSPGetType(ksp, &type);
    
    std::cout<< "SNES type da passo steady exec:"<<ttype<<std::endl;
    std::cout<< "KSP type da passo steady exec:"<<type<<std::endl;
    
    _ksp_ptr = (KSP_PARROT *)ksp->data;
    (_ksp_ptr[0].local_pc)=&_problem_PC;
    (*_ksp_ptr).factorized=&_factorized;
};
