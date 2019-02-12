
#include "Parrot_preonly.h"
#include "ksp_parrot_impl.h"
#include "iostream"
#include "chrono"
//static
PetscErrorCode KSPSetUp_Parrot_PREONLY(KSP ksp)
{
    std::cout<<"called SetUp\n";
    PetscFunctionBegin;
    std::cout<<"done SetUp\n";
    PetscFunctionReturn(0);
}

//static
PetscErrorCode  KSPSolve_Parrot_PREONLY(KSP ksp)
{
    std::cout<<"called Solve\n";
  PetscErrorCode ierr;
  PetscBool      diagonalscale;
  PCFailedReason pcreason;
    
  PetscFunctionBegin;
  ierr = PCGetDiagonalScale(ksp->pc,&diagonalscale);CHKERRQ(ierr);
  if (diagonalscale) SETERRQ1(PetscObjectComm((PetscObject)ksp),PETSC_ERR_SUP,"Krylov method %s does not support diagonal scaling",((PetscObject)ksp)->type_name);
  if (!ksp->guess_zero) SETERRQ(PetscObjectComm((PetscObject)ksp),PETSC_ERR_USER,"Running KSP of preonly doesn't make sense with nonzero initial guess\n\
               you probably want a KSP type of Richardson");
  ksp->its = 0;
    
    KSP_PARROT * _ksp_ptr;
    _ksp_ptr = (KSP_PARROT *)ksp->data;
    int factorized=_ksp_ptr[0].factorized[0];
    std::cout<<factorized<<std::endl;
    Mat Hmat,Pmat;
    
    ierr = KSPGetOperators(ksp,&Hmat,&Pmat);CHKERRQ(ierr);
    
    if (factorized==0)
    {
        _ksp_ptr[0].factorized[0]=1;
    
    
    //PC pc_lu;
//    PCCreate(PetscObjectComm((PetscObject)ksp), _ksp_ptr[0].local_pc);
//    PCSetType(_ksp_ptr[0].local_pc[0],PCLU);
    PCSetOperators(_ksp_ptr[0].local_pc[0],Hmat,Pmat);
        PCFactorSetMatSolverPackage(_ksp_ptr[0].local_pc[0],MATSOLVERSUPERLU_DIST);
    std::cout<<"start factorizing?\n";
    auto t_start = std::chrono::high_resolution_clock::now();
    PCSetUp(_ksp_ptr[0].local_pc[0]);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout<<"done factorizing?\n";
        std::cout<<"fact time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";
    }
    PCSetReusePreconditioner(_ksp_ptr[0].local_pc[0], PETSC_TRUE);
    std::cout<<"start solving?\n";
    auto t_start = std::chrono::high_resolution_clock::now();
    PCApply(_ksp_ptr[0].local_pc[0],ksp->vec_rhs,ksp->vec_sol);
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout<<"done solving?\n";
    std::cout<<"solve time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";

    _ksp_ptr[0].local_pc=NULL;
    
//    std::cout<<"start solving?\n";
//     t_start = std::chrono::high_resolution_clock::now();
//    PCApply(_ksp_ptr[0].local_pc[0],ksp->vec_rhs,ksp->vec_sol);
//     t_end = std::chrono::high_resolution_clock::now();
//    std::cout<<"done solving?\n";
//    std::cout<<"solve time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";
    
    Vec r;
    VecDuplicate(ksp->vec_rhs,&r);
    MatResidual(Hmat,ksp->vec_rhs,ksp->vec_sol,r);
    PetscReal norm;
    VecNorm(r,NORM_2,&norm);
    std::cout<<"qui "<<norm<<std::endl;
    PetscPrintf(PETSC_COMM_WORLD,"   %14.12e \n", norm);

    
//  ierr     = KSP_PCApply(ksp,ksp->vec_rhs,ksp->vec_sol);CHKERRQ(ierr);
//  ierr     = PCGetSetUpFailedReason(ksp->pc,&pcreason);CHKERRQ(ierr);
//  if (pcreason) {
//    ksp->reason = KSP_DIVERGED_PCSETUP_FAILED;
//  } else {
    ksp->its    = 1;
    ksp->reason = KSP_CONVERGED_ITS;
//  }
    
  PetscFunctionReturn(0);
}

/*MC
     KSPPREONLY - This implements a stub method that applies ONLY the preconditioner.
                  This may be used in inner iterations, where it is desired to
                  allow multiple iterations as well as the "0-iteration" case. It is
                  commonly used with the direct solver preconditioners like PCLU and PCCHOLESKY

   Options Database Keys:
.   -ksp_type preonly

   Level: beginner

   Notes: Since this does not involve an iteration the basic KSP parameters such as tolerances and iteration counts
          do not apply

   Developer Notes: Even though this method does not use any norms, the user is allowed to set the KSPNormType to any value.
    This is so the users does not have to change KSPNormType options when they switch from other KSP methods to this one.

.seealso:  KSPCreate(), KSPSetType(), KSPType (for list of available types), KSP

M*/

PETSC_EXTERN PetscErrorCode KSPCreate_Parrot_PREONLY(KSP ksp)
{
    std::cout<<"called Create\n";
  PetscErrorCode ierr;
    
    KSP_PARROT       *ksp_parrot;

  PetscFunctionBegin;
    
    ierr = PetscNewLog(ksp,&ksp_parrot);CHKERRQ(ierr);
    
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_LEFT,3);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_RIGHT,2);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_LEFT,2);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_RIGHT,2);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_LEFT,2);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_RIGHT,2);CHKERRQ(ierr);
  ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NATURAL,PC_LEFT,2);CHKERRQ(ierr);

  ksp->data                = ksp_parrot;
  ksp->ops->setup          = KSPSetUp_Parrot_PREONLY;
  ksp->ops->solve          = KSPSolve_Parrot_PREONLY;
  ksp->ops->destroy        = KSPDestroyDefault;
  ksp->ops->buildsolution  = KSPBuildSolutionDefault;
  ksp->ops->buildresidual  = KSPBuildResidualDefault;
  ksp->ops->setfromoptions = 0;
  ksp->ops->view           = 0;
std::cout<<"done Create\n";
    PetscFunctionReturn(0);
}
