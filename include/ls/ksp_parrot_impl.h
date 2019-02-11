#ifndef _KSPPARROTIMPL_H
#define _KSPPARROTIMPL_H

#include <petsc/private/kspimpl.h>
//#include <petsc/private/snesimpl.h>

typedef struct {
    PC * local_pc;
    int * factorized;
} KSP_PARROT;

#endif
