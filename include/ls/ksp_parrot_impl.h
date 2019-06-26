#ifndef _KSPPARROTIMPL_H
#define _KSPPARROTIMPL_H
#include "FEProblem.h"
#include <petsc/private/kspimpl.h>
//#include <petsc/private/snesimpl.h>

typedef struct {
    PC * local_pc;
    FEProblem * _fe_problem;
    int * factorized;
} KSP_PARROT;

#endif
