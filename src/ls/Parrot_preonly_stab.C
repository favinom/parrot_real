
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
    
    std::string userobject_name_2 = "Matrix_Storage_1";

    std::string userobject_name_3 = "Matrix_Storage_2";

    std::string userobject_name_4 = "Matrix_Storage_3";
    
    std::string userobject_name_5 = "Matrix_Storage_4";
    
    std::cout<<"FACTORIZED:"<<factorized<<std::endl;

    utopia::DSMatrixd A_m_t, D_stab,  A_m_tot, mass_p, mass_lumped_p, mass, mass_lumped, J_tot, A_stab;
    
    utopia::DVectord res_old, res, c_m, c_m_old, mass_c_m_old, mass_c_m_old_dot, rhs_m_t;

    
    if (factorized==0)
    {
        
        Mat Hmat,Pmat;


        ierr = KSPGetOperators(ksp,&Hmat,&Pmat);CHKERRQ(ierr);
        
        _ksp_ptr[0].factorized[0]=1;
    
        utopia::convert(Hmat, A_m_t);
        
        D_stab = A_m_t;
        
        D_stab*=0;
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).stabilize_A_matrix(*_ksp_ptr->_fe_problem, A_m_t, D_stab);
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).assemble_poro_mass_matrix(*_ksp_ptr->_fe_problem, mass_p, mass_lumped_p);
        
        utopia::DSMatrixd A_stab = D_stab + A_m_t;
        
        Real dt = static_cast<Transient*>(_ksp_ptr->_fe_problem->getMooseApp().getExecutioner())->getDT();
        
        double inv_dt= 1.0/dt;

        utopia::DSMatrixd A_ns_tot = mass_lumped_p  + A_m_t;
    
        A_m_tot =  1.0 *  A_stab + mass_lumped_p * inv_dt;

        utopia::convert(ksp->vec_rhs, res_old);
        
        res = res_old;
        
        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).constraint_mat(res, A_m_tot, true);

        utopia::DSMatrixd mass_temp = mass_lumped_p;

        const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).constraint_mat(res, mass_temp, true);
        
        const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).setMatrix() = std::make_shared<utopia::DSMatrixd>(A_m_tot);
        
        const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_3)).setMatrix() = std::make_shared<utopia::DSMatrixd>(D_stab);
        
        const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_4)).setMatrix() = std::make_shared<utopia::DSMatrixd>(A_stab);
        
        const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_5)).setMatrix() = std::make_shared<utopia::DSMatrixd>(mass_lumped_p);

        // utopia::write("Old_J_ns.m",A_ns_tot);

        // utopia::write("Old_J_s.m", A_m_tot);

        // utopia::write("Old_M.m", mass_temp);

        // utopia::write("A_m_t.m", A_m_t);  

        // utopia::write("A_stab.m",  A_stab); 

        // utopia::write("mass_lumped_p.m",   mass_lumped_p); 
        
        std::cout<<"start factorizing?\n";
        
        auto t_start = std::chrono::high_resolution_clock::now();
        
        const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).setRHS() = std::make_shared<utopia::DVectord>(res_old);
        
        auto t_end = std::chrono::high_resolution_clock::now();

        std::cout<<"done factorizing?\n";
        
        std::cout<<"fact time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";

    }
    

    auto op = std::make_shared<utopia::Factorization<utopia::DSMatrixd, utopia::DVectord> >(MATSOLVERSUPERLU,PCLU);

    //auto op = std::make_shared<utopia::GMRES<utopia::DSMatrixd, utopia::DVectord> >();

    auto J_store = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getMatrix();
    
    J_tot = * J_store;

    auto D_store = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_3)).getMatrix();
    
    D_stab = * D_store;

    auto A_store = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_4)).getMatrix();
    
    A_stab = * A_store;
    
    auto mass_lumped_p_store = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_5)).getMatrix();
    
    mass_lumped_p = * mass_lumped_p_store;
    
    c_m =  utopia::local_zeros(size(J_tot).get(0));

    c_m_old =  utopia::local_zeros(size(J_tot).get(0));
    
    auto sol_m_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_ksp_ptr->_fe_problem->es().get_system<TransientNonlinearImplicitSystem>("nl0").old_local_solution.get())->vec();
    
    utopia::convert(sol_m_old, c_m_old);
    
    // utopia::disp(c_m_old);
    
    Real dt = static_cast<Transient*>(_ksp_ptr->_fe_problem->getMooseApp().getExecutioner())->getDT();
    
    double inv_dt= 1.0/dt;

    auto rhs_m_c = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getRHS();
    
    //const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).assemble_poro_mass_matrix(*_ksp_ptr->_fe_problem, mass_p, mass_lumped_p);
  
    utopia::DSMatrixd mass_lumped_pc =  mass_lumped_p;
    
    std::cout<<"done factorizing?\n";
    
    const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).constraint_mat(*rhs_m_c,  mass_lumped_pc, true);
    
    mass_c_m_old = mass_lumped_pc * c_m_old;
    
    mass_c_m_old_dot = mass_c_m_old * inv_dt;
    
    rhs_m_t =  - 1.0 * mass_c_m_old_dot;
    
    const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).boundary_constraint_vec(*rhs_m_c,  rhs_m_t, true);
    
    res = rhs_m_t;
    
    /*utopia::write("J_tot.m", J_tot); utopia::write("res.m", res);*/
 
    std::cout<<"time_step"<<_ksp_ptr->_fe_problem->timeStep()<<std::endl;

    // utopia::disp(c_m_old);
  
    op->update(make_ref(J_tot));

    // utopia::DVectord res_stab = 1.0 * A_stab * c_m;

    // auto rhs_m_c = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getRHS();

    // const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).unconstraint_concentration_vec(*rhs_m_c, res_stab, true);

    // res+=res_stab;

    op->apply(res, c_m);  //utopia::write("c_m.m",c_m);

  /*  const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).assemble_mass_matrix(*_ksp_ptr->_fe_problem, mass, mass_lumped);

  

    utopia::DVectord _f_s;

    utopia::DSMatrixd mass_lumped_c = mass_lumped;

    const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).constraint_mat(*rhs_m_c,  mass_lumped_c, true);
    


    utopia::DVectord diag_elem = 1./sum((mass_lumped_c),1);

    utopia::DSMatrixd inv_mass_lumped_mat = diag(diag_elem);


    utopia::DVectord c_dot = 1.0 * inv_mass_lumped_mat * A_stab * c_m;

    c_dot*=-1.0;

    const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).stabilize_coeffiecient(*rhs_m_c, *_ksp_ptr->_fe_problem, c_m, c_dot, D_stab, mass, _f_s);

    // utopia::DVectord _f_s = (mass - mass_lumped) * (c_m_old - c_m) * inv_dt + 1.0 * D_stab * c_m;

    // auto rhs_m_c = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getRHS();

    //const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).unconstraint_concentration_vec(*rhs_m_c, _f_s, true);
   

    utopia::DVectord diag_elem_c = 1./sum((mass_lumped_c),1);

    utopia::DSMatrixd inv_mass_lumped_mat_c = diag(diag_elem);

    utopia::DVectord c_m_tot = c_m + 1.0 * inv_mass_lumped_mat_c * _f_s * dt;

   */

    //utopia::DVectord r_m =  utopia::local_zeros(size(J_tot).get(0));

    //const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).CopyMatrixSolution(c_m);
    
    ksp->its    = 1;

    ksp->reason = KSP_CONVERGED_ITS;

   // utopia::DVectord r_m =  utopia::local_zeros(size(J_tot).get(0));

    Vec r;

    VecDuplicate(ksp->vec_rhs,&r);

    utopia::DVectord r_m =  utopia::local_zeros(size(J_tot).get(0));

    convert(r, r_m);

    r_m*=0;

    convert(r_m, ksp->vec_sol);

    const_cast<FractureAppConforming&>(_ksp_ptr->_fe_problem->getUserObject<FractureAppConforming>(userobject_name_1)).CopyMatrixSolution(c_m);

    // Vec r;
    // VecDuplicate(ksp->vec_rhs,&r);

    // utopia::DVectord r_m =  utopia::local_zeros(size(J_tot).get(0));

    // convert(r, r_m);

    // r_m*=0;

    // convert(r_m, ksp->vec_rhs);
    
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



// #include "Parrot_preonly_stab.h"
// #include "ksp_parrot_impl.h"
// #include "iostream"
// #include "chrono"
// #include "utopia.hpp"
// #include "Test_Assembly.h"
// #include "SystemInitialize.h"
// #include "StoreTransferOperators.h"
// #include "libmesh/petsc_matrix.h"
// #include <petscksp.h>
// #include <petscmat.h>
// #include "FEProblem.h"
// #include "PetscSupport.h"
// #include "libmesh/nonlinear_solver.h"
// #include "libmesh/petsc_nonlinear_solver.h"
// #include "libmesh/petsc_vector.h"
// #include "libmesh/petsc_matrix.h"
// #include "Transient.h"
// #include "NonlinearSystemBase.h"
// //static
// PetscErrorCode KSPSetUp_Parrot_PREONLY_STAB(KSP ksp)
// {
//     std::cout<<"called SetUp\n";
//     PetscFunctionBegin;
//     std::cout<<"done SetUp\n";
//     PetscFunctionReturn(0);
// }

// //static
// PetscErrorCode  KSPSolve_Parrot_PREONLY_STAB(KSP ksp)
// {
//   std::cout<<"called Solve\n";
//   PetscErrorCode ierr;
//   PetscBool      diagonalscale;
//   PCFailedReason pcreason;
    
//   PetscFunctionBegin;
//   ierr = PCGetDiagonalScale(ksp->pc,&diagonalscale);CHKERRQ(ierr);
//   if (diagonalscale) SETERRQ1(PetscObjectComm((PetscObject)ksp),PETSC_ERR_SUP,"Krylov method %s does not support diagonal scaling",((PetscObject)ksp)->type_name);
//   if (!ksp->guess_zero) SETERRQ(PetscObjectComm((PetscObject)ksp),PETSC_ERR_USER,"Running KSP of preonly doesn't make sense with nonzero initial guess\n\
//                you probably want a KSP type of Richardson");
//   ksp->its = 0;
    
//     KSP_PARROT * _ksp_ptr;
//     _ksp_ptr = (KSP_PARROT *)ksp->data;
//     int factorized=_ksp_ptr[0].factorized[0];
//     std::cout<<factorized<<std::endl;
//     Mat Hmat,Pmat;
    
//     ierr = KSPGetOperators(ksp,&Hmat,&Pmat);CHKERRQ(ierr);

//     utopia::DSMatrixd H_m_t, P_m_t, J_tot, P_tot, D_stab, P_stab, mass, mass_lumped;

//      utopia::DVectord res, res_old;

//     std::string userobject_name_1 = "FractureApp";

//     std::string userobject_name_2 = "Matrix_Storage_1";

//     std::string userobject_name_3 = "Matrix_Storage_2";
    
//     if (factorized==0)
//     {
//         _ksp_ptr[0].factorized[0]=1;
    
  
//     utopia::convert(Hmat, H_m_t);

//     utopia::convert(Pmat, P_m_t);
    
//     D_stab = H_m_t;

//     //utopia::write("New_M.m",H_m_t);
    
//     D_stab*=0;

//     P_stab = P_m_t;
    
//     P_stab*=0;

//     utopia::convert(ksp->vec_rhs, res);

//     const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).unconstraint_mat(res, H_m_t, true);

//     const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).stabilize_A_matrix(*_ksp_ptr->_fe_problem, H_m_t, D_stab);
//     //const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).stabilize_A_matrix(*_ksp_ptr->_fe_problem, P_m_t, P_stab);
   
//     J_tot = H_m_t + D_stab;

//     const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_3)).setMatrix() = std::make_shared<utopia::DSMatrixd>(D_stab);

//     const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).constraint_mat(res, J_tot, true);
    
//     utopia::write("J_tot.m",J_tot);

//     res_old = res;

//     const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).setMatrix() = std::make_shared<utopia::DSMatrixd>(J_tot);
    
//     const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).setRHS() = std::make_shared<utopia::DVectord>(res_old);


//     PCSetOperators(_ksp_ptr[0].local_pc[0], utopia::raw_type(J_tot),utopia::raw_type(J_tot));

//     PCSetType(_ksp_ptr[0].local_pc[0],PCLU);

//     PCFactorSetMatSolverPackage(_ksp_ptr[0].local_pc[0],MATSOLVERMUMPS);

//     std::cout<<"start factorizing?\n";
//     auto t_start = std::chrono::high_resolution_clock::now();
//     PCSetUp(_ksp_ptr[0].local_pc[0]);
//     auto t_end = std::chrono::high_resolution_clock::now();
//     std::cout<<"done factorizing?\n";
//     std::cout<<"fact time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";
//     }

    
 

  
//     //NumericVector<libMesh::Number> & res_time = _ksp_ptr->_fe_problem->getNonlinearSystemBase().getResidualTimeVector();
   
//     auto A_store = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getMatrix();

//     J_tot = *A_store;

//     PCSetReusePreconditioner(_ksp_ptr[0].local_pc[0], PETSC_TRUE);
//     std::cout<<"start solving?\n";
//     auto t_start = std::chrono::high_resolution_clock::now();

//     if(_ksp_ptr->_fe_problem->timeStep()>1){

//        utopia::DVectord temp;

//        utopia::convert(ksp->vec_rhs, temp);

//        res =  J_tot * temp;
//     }

//     else{

//       utopia::DVectord temp;

//       utopia::convert(ksp->vec_rhs, temp);

//       res = temp;

//     }
    
//     // auto rhs_m_c = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getRHS();
//     // const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).boundary_constraint_vec(*rhs_m_c,  res, true);
//     PCApply(_ksp_ptr[0].local_pc[0],utopia::raw_type(res),ksp->vec_sol);
 
//     auto t_end = std::chrono::high_resolution_clock::now();
//     std::cout<<"done solving?\n";
//     std::cout<<"solve time: "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()<< " ms\n";

//     _ksp_ptr[0].local_pc=NULL;
    
    
//     Vec r;
//     VecDuplicate(ksp->vec_rhs,&r);

//     // utopia::DVectord c_m, c_m_old, c_m_oldest;

//     // if (const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getOldSol()!=NULL){
//     //      c_m_oldest = *const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getOldSol();
//     //    }

//     // else{

//     //      utopia::convert(ksp->vec_sol,c_m_oldest);
//     // }   
    

//     // utopia::convert(ksp->vec_sol,c_m_old);

//     // const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).setOldSol() = std::make_shared<utopia::DVectord>(c_m_old);

  

//     MatResidual(utopia::raw_type(J_tot),utopia::raw_type(res),ksp->vec_sol,r);

//     PetscReal norm;

//     VecNorm(r,NORM_2,&norm);

//     std::cout<<"qui "<<norm<<std::endl;

//     PetscPrintf(PETSC_COMM_WORLD,"   %14.12e \n", norm);

//     // const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).assemble_mass_matrix(*_ksp_ptr->_fe_problem, mass, mass_lumped);

//     // auto D_Store = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_3)).getMatrix();

//     // D_stab = *D_Store;

//     // utopia::convert(ksp->vec_rhs,c_m);

//     // const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).assemble_mass_matrix(*_ksp_ptr->_fe_problem, mass, mass_lumped);

//     // Real dt = static_cast<Transient*>(_ksp_ptr->_fe_problem->getMooseApp().getExecutioner())->getDT();
        

//     // double inv_dt= 1.0/dt;
    

//     // utopia::DVectord _f_s = (mass - mass_lumped) * (c_m + c_m_oldest) * inv_dt + 1.0 * D_stab * (c_m - c_m_old);

//     // auto rhs_m_c = const_cast<StoreTransferOperators&>(_ksp_ptr->_fe_problem->getUserObject<StoreTransferOperators>(userobject_name_2)).getRHS();

//     // const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).unconstraint_concentration_vec(*rhs_m_c, _f_s, true);

//     // utopia::DVectord diag_elem = 1./sum((mass_lumped),1);

//     // utopia::DSMatrixd inv_mass_lumped_mat = diag(diag_elem);

//     // c_m +=1.0 * inv_mass_lumped_mat * _f_s;
    
//     // c_m += _f_s; 


    
//     // auto op = std::make_shared<utopia::Factorization<utopia::DSMatrixd, utopia::DVectord> >(MATSOLVERMUMPS,PCLU);

//     // utopia::convert(c_m, ksp->vec_sol);

//     // utopia::disp(res);
     
//     // op->update(make_ref(J_tot));

//     // op->apply(res, c_m);

//     // const_cast<SystemInitialize&>(_ksp_ptr->_fe_problem->getUserObject<SystemInitialize>(userobject_name_1)).CopyMatrixSolution(c_m);

//     ksp->its    = 1;
      
//     ksp->reason = KSP_CONVERGED_ITS;

//     PetscFunctionReturn(0);

// }


// PETSC_EXTERN PetscErrorCode KSPCreate_Parrot_PREONLY_STAB(KSP ksp)
// {
//   std::cout<<"called Create\n";
//   PetscErrorCode ierr;
    
//   KSP_PARROT       *ksp_parrot;

//   PetscFunctionBegin;
    
//   ierr = PetscNewLog(ksp,&ksp_parrot);CHKERRQ(ierr);
    
//   ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_LEFT,3);CHKERRQ(ierr);
//   ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NONE,PC_RIGHT,2);CHKERRQ(ierr);
//   ierr = KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_LEFT,2);CHKERRQ(ierr);
//   ierr = KSPSetSupportedNorm(ksp,KSP_NORM_PRECONDITIONED,PC_RIGHT,2);CHKERRQ(ierr);
//   ierr = KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_LEFT,2);CHKERRQ(ierr);
//   ierr = KSPSetSupportedNorm(ksp,KSP_NORM_UNPRECONDITIONED,PC_RIGHT,2);CHKERRQ(ierr);
//   ierr = KSPSetSupportedNorm(ksp,KSP_NORM_NATURAL,PC_LEFT,2);CHKERRQ(ierr);

//   ksp->data                = ksp_parrot;
//   ksp->ops->setup          = KSPSetUp_Parrot_PREONLY_STAB;
//   ksp->ops->solve          = KSPSolve_Parrot_PREONLY_STAB;
//   ksp->ops->destroy        = KSPDestroyDefault;
//   ksp->ops->buildsolution  = KSPBuildSolutionDefault;
//   ksp->ops->buildresidual  = KSPBuildResidualDefault;
//   ksp->ops->setfromoptions = 0;
//   ksp->ops->view           = 0;
//   std::cout<<"done Create\n";
//   PetscFunctionReturn(0);
// }

