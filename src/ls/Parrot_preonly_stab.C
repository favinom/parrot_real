
#include "Parrot_preonly_stab.h"
#include "ksp_parrot_impl.h"
#include "iostream"
#include "chrono"
#include "utopia.hpp"
#include "Test_Assembly.h"
#include "StoreTransferOperators.h"
#include "libmesh/petsc_matrix.h"
#include <petscksp.h>
#include <petscmat.h>
#include "FEProblem.h"
#include "PetscSupport.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "Transient.h"
#include "NonlinearSystemBase.h"
//static
PetscErrorCode KSPSetUp_Parrot_PREONLY_STAB(KSP ksp)
{
    std::cout<<"called SetUp\n";
    PetscFunctionBegin;
    std::cout<<"done SetUp\n";
    PetscFunctionReturn(0);
}

//static
PetscErrorCode  KSPSolve_Parrot_PREONLY_STAB(KSP ksp)
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
    
    std::string userobject_name_1 = "FractureApp";
    
    std::string userobject_name_2 = "Matrix_System";
    
    std::cout<<"FACTORIZED:"<<factorized<<std::endl;

    utopia::DSMatrixd A_m_t, D_stab,  A_m_tot, mass_p, mass_lumped_p;
    
    utopia::DVectord res_old, res, c_m, c_m_old, mass_c_m_old, mass_c_m_old_dot, rhs_m_t;
    
    if (factorized==0)
    {
        
        Mat Hmat,Pmat;
        
        ierr = KSPGetOperators(ksp,&Hmat,&Pmat);CHKERRQ(ierr);
        
        _ksp_ptr[0].factorized[0]=1;
    
        utopia::convert(Hmat, A_m_t);
        
        D_stab = A_m_t;
        
        D_stab*=0;
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).stabilize_A_matrix(*_ksp_ptr->_fe_problem, D_stab);
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).assemble_poro_mass_matrix(*_ksp_ptr->_fe_problem, mass_p, mass_lumped_p);
        
        utopia::DSMatrixd A_stab = D_stab + A_m_t;
        
        Real dt = static_cast<Transient*>(_ksp_ptr->_fe_problem->getMooseApp().getExecutioner())->getDT();
        
        double inv_dt= 1.0/dt;
    
        A_m_tot =  1.0 *  A_stab + mass_lumped_p * inv_dt;
        
        disp(A_m_tot);
        
        const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).setMatrix() = std::make_shared<utopia::DSMatrixd>(A_m_tot);
    
        utopia::convert(ksp->vec_rhs, res_old);
        
        res = res_old;
        
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).constraint_mat(res, A_m_tot, true);
        
        
        
        std::cout<<"start factorizing?\n";
        
        auto t_start = std::chrono::high_resolution_clock::now();
        
        const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).setRHS() = std::make_shared<utopia::DVectord>(res_old);
        
        auto t_end = std::chrono::high_resolution_clock::now();
        std::cout<<"done factorizing?\n";
        
        std::cout<<"fact time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";
    }
    
    if (_ksp_ptr->_fe_problem->timeStep()>1)
    {
        
        
        
        auto sol_m_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_ksp_ptr->_fe_problem->es().get_system<TransientNonlinearImplicitSystem>("nl0").old_local_solution.get())->vec();
        
        utopia::convert(sol_m_old, c_m_old);
        
        utopia::disp(c_m_old);
        
        Real dt = static_cast<Transient*>(_ksp_ptr->_fe_problem->getMooseApp().getExecutioner())->getDT();
        
        double inv_dt= 1.0/dt;
    
        auto rhs_m_c = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getRHS();
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).assemble_poro_mass_matrix(*_ksp_ptr->_fe_problem, mass_p, mass_lumped_p);
        
        utopia::DSMatrixd mass_lumped_pc =  mass_lumped_p;
        
        std::cout<<"done factorizing?\n";
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).constraint_mat(*rhs_m_c,  mass_lumped_pc, true);
        
        mass_c_m_old = mass_lumped_pc * c_m_old;
        
        mass_c_m_old_dot = mass_c_m_old * inv_dt;
        
        rhs_m_t =  - 1.0 * mass_c_m_old_dot;
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).boundary_constraint_vec(*rhs_m_c,  rhs_m_t, true);
        
        res = rhs_m_t;
        
        std::cout<<"time_step"<<_ksp_ptr->_fe_problem->timeStep()<<std::endl;

    }
    
    auto op = std::make_shared<utopia::Factorization<utopia::DSMatrixd, utopia::DVectord> >(MATSOLVERMUMPS,PCLU);

    c_m =  utopia::local_zeros(local_size(res));
    
    auto A_store = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getMatrix();
    
    A_m_tot = * A_store;

    op->update(make_ref(A_m_tot));

    op->apply(res, c_m);
    
//    utopia::disp(c_m);

    const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).CopyMatrixSolution(c_m);
    
    ksp->its    = 1;
    ksp->reason = KSP_CONVERGED_ITS;
    
    PetscFunctionReturn(0);

}


PETSC_EXTERN PetscErrorCode KSPCreate_Parrot_PREONLY_STAB(KSP ksp)
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
  ksp->ops->setup          = KSPSetUp_Parrot_PREONLY_STAB;
  ksp->ops->solve          = KSPSolve_Parrot_PREONLY_STAB;
  ksp->ops->destroy        = KSPDestroyDefault;
  ksp->ops->buildsolution  = KSPBuildSolutionDefault;
  ksp->ops->buildresidual  = KSPBuildResidualDefault;
  ksp->ops->setfromoptions = 0;
  ksp->ops->view           = 0;
  std::cout<<"done Create\n";
  PetscFunctionReturn(0);
}
