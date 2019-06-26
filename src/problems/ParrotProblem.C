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
    PCCreate(PETSC_COMM_WORLD, &_problem_PC);
    PCSetType(_problem_PC,PCLU);
}

void ParrotProblem::initialSetup()
{
    std::cout<<"BEGIN ParrotProblem::initialSetup"<<std::endl;
    
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
    
    std::cout<<"END ParrotProblem::initialSetup"<<std::endl;
    
};

void ParrotProblem::timestepSetup()
{
    std::cout<<"BEGIN ParrotProblem::timestepSetup"<<std::endl;
    
    FEProblem::timestepSetup();
    
    PetscErrorCode ierr;
    Moose::PetscSupport::petscSetOptions(*this);

    PetscNonlinearSolver<Number> * petsc_solver =
    dynamic_cast<PetscNonlinearSolver<Number> *>((*_nl_sys).nonlinearSolver());
    SNES snes = petsc_solver->snes();
    KSP ksp;
    PC pc;
    SNESType ttype;
    KSPType type;
    PCType ptype;
    
    ierr = SNESGetKSP(snes,&ksp);
    ierr = KSPGetPC(ksp,&pc);
    ierr = SNESSetFromOptions(snes);
    ierr = SNESGetType(snes, &ttype);
    ierr = KSPGetType(ksp, &type);
    ierr = PCGetType(pc, &ptype);
    
    std::cout<< "SNES type da passo steady exec:"<<ttype<<std::endl;
    std::cout<< "KSP type da passo steady exec:"<<type<<std::endl;
    std::cout<< "PC type da passo steady exec:"<<ptype<<std::endl;
    
    _ksp_ptr = (KSP_PARROT *)ksp->data;
    (_ksp_ptr[0].local_pc)=&_problem_PC;
    (*_ksp_ptr).factorized=&_factorized;
    
    _ksp_ptr->_fe_problem = this;
    
    std::cout<<"END ParrotProblem::timestepSetup"<<std::endl;
};

void ParrotProblem::solve()
{
    std::cout<<"BEGIN ParrotProblem::solve"<<std::endl;
    
    Moose::perf_log.push("solve()", "Execution");
    
#ifdef LIBMESH_HAVE_PETSC
    Moose::PetscSupport::petscSetOptions(*this); // Make sure the PETSc options are setup for this app
#endif
    
    Moose::setSolverDefaults(*this);
    
    // Setup the output system for printing linear/nonlinear iteration information
    initPetscOutput();
    
    possiblyRebuildGeomSearchPatches();
    
    // reset flag so that linear solver does not use
    // the old converged reason "DIVERGED_NANORINF", when
    // we throw  an exception and stop solve
    //_fail_next_linear_convergence_check = false;
    
    if (_solve)
    {
        std::cout<<"before solve"<<std::endl;
        _nl->solve();
                std::cout<<"after solve"<<std::endl;
    }
    
    if (_solve)
        _nl->update();
    
    // sync solutions in displaced problem
    //if (_displaced_problem)
    //    _displaced_problem->syncSolutions();
    
    Moose::perf_log.pop("solve()", "Execution");
    
    std::cout<<"END ParrotProblem::solve"<<std::endl;
}
