/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*    Immersed_Boundary- ICS Mechanical simulation framework    */
/*                Prepared by Maria Nestola,                    */
/*                  ICS, USI, 6900 Lugano                       */
/*                                                              */
/****************************************************************/

#include "utopia_ForwardDeclarations.hpp"
#include <algorithm>    // std::max
#include "utopia.hpp"

#include "Test_Assembly.h"
#include "FEProblem.h"
#include "AddVariableAction.h"
#include "NonlinearSystemBase.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"
#include "DisplacedProblem.h"
#include <opencl_adapter.hpp>
#include "libmesh/linear_partitioner.h"
#include "libmesh/petsc_matrix.h"
#include "utopia_assemble_volume_transfer.hpp"
#include <petscksp.h>
#include <petscmat.h>
#include "MooseVariable.h"
#include <utopia.hpp>


#include <queue>
#include "Transient.h"
#include <iostream>


#include "libmesh/boundary_info.h"
#include "libmesh/boundary_mesh.h"

#include "ksp_parrot_impl.h"



typedef utopia::DSMatrixd SparseMatT;
typedef utopia::DVectord VecT;


registerMooseObject("parrot_realApp",FractureAppConforming);

template <>
InputParameters
validParams<FractureAppConforming>()
{
    
    InputParameters params = validParams<GeneralUserObject>();
    params.addRequiredParam<VariableName>("matrix_variable",
                                          "The variable to transfer from.");
    params.addParam<bool>("constraint_m","false","put true if matrix has Dirichlet BC");
    return params;
    
}

FractureAppConforming::FractureAppConforming(const InputParameters & parameters):
GeneralUserObject(parameters),
_m_var_name(getParam<VariableName>("matrix_variable")),
//_operator_storage(getUserObject<StoreTransferOperators>("operator_userobject")),
//_transport(getParam<bool>("transport")),
_constraint_m(getParam<bool>("constraint_m"))

{
    
}




void
FractureAppConforming::initialize()
{
    _console << "Initial Setup of Fracture App " << std::endl;
    
  if(!_fe_problem.hasUserObject(_userobject_name_1)){
      std::string class_name = "StoreTransferOperators";
      auto params = _fe_problem.getMooseApp().getFactory().getValidParams(class_name);
      params.set<bool>("use_displaced_mesh") = false;
      params.set<ExecFlagEnum>("execute_on") = "initial";
      
     _fe_problem.addUserObject("StoreTransferOperators", _userobject_name_1, params);
  }
    
//    assemble_mass_matrix(_fe_problem,  mass, mass_lumped);
//
//    assemble_poro_mass_matrix(_fe_problem,  mass_p, mass_lumped_p);
    
    
    
}


void
FractureAppConforming::execute()
{
//    solve_stabilize_monolithic();
}




void
FractureAppConforming::boundary_constraint_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
{
    
    
    using namespace utopia;
    
    {
        Write<utopia::UVector> w_v(vec);
        
        Read<utopia::UVector> r_v(boundary);
        
        if(has_constaints) {
            
            Range r = range(vec);
            
            for(SizeType i = r.begin(); i < r.end(); ++i) {
                
                
                
                if(boundary.get(i)!=0) {
                    
                    
                    vec.set(i, boundary.get(i));
                    
                }
            }
        }
    }
    
    synchronize(vec);
    
}


void
FractureAppConforming::zero_constraint_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
{
    
    
    using namespace utopia;
    
    {
        Write<utopia::UVector> w_v(vec);
        
        Read<utopia::UVector> r_v(boundary);
        
        if(has_constaints) {
            
            Range r = range(vec);
            
            for(SizeType i = r.begin(); i < r.end(); ++i) {
                
                
                
                if(boundary.get(i)!=0) {
                    
                    
                    vec.set(i, 0.0);
                    
                }
            }
        }
    }
    
    synchronize(vec);
    
}

void
FractureAppConforming::one_constraint_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
{
    
    
    using namespace utopia;
    
    {
        Write<utopia::UVector> w_v(vec);
        
        Read<utopia::UVector> r_v(boundary);
        
        if(has_constaints) {
            
            Range r = range(vec);
            
            for(SizeType i = r.begin(); i < r.end(); ++i) {
                
                
                
                if(boundary.get(i)!=0) {
                    
                    
                    vec.set(i, 1.0);
                    
                }
            }
        }
    }
    
    synchronize(vec);
    
}


void
FractureAppConforming::constraint_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints)
{
    //    auto &V = space->space().last_subspace();
    
    using namespace utopia;
    
    typedef UTOPIA_SIZE_TYPE(UVector) SizeType;
    
    std::vector<SizeType> rows;
    
    {
        
        Read<utopia::UVector> r_v(boundary);
        
        rows.reserve(local_size(mat).get(0));
        
        if(has_constaints) {
            
            auto r = row_range(mat);
            
            for(SizeType i = r.begin(); i < r.end(); ++i) {
                
                //std::cout<<"value "<< boundary.get(i) <<std::endl;
                
                if(boundary.get(i)!=0) {
                    
                    //std::cout<<"i "<< i <<"valpos->second "<<boundary.get(i)<<std::endl;
                    
                    // if(valpos != rhs_values.end()) {
                    rows.push_back(i);
                    // }
                }
            }
        }
        
        set_zero_rows(mat, rows, 1.);
    }
}




void
FractureAppConforming::assemble_mass_matrix(FEProblemBase & _problem, utopia::USparseMatrix &mass_matrix, utopia::USparseMatrix &lumped_mass_matrix){
    
    _console << "Assemble_Mass_matrix() begin "  << std::endl;

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _problem.es().get_mesh();

    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to our system.
    auto & _system = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _system.get_dof_map().variable_type(0);

    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

    SparseMatrix<Number> & matrix_M = *_system.matrix;

    // SparseMatrix<Number> & matrix_M_p = *_system.matrix;

    //_console << "is matrix closed: " << matrix_A.closed() << std::endl;

    // The element mass matrix.
    DenseMatrix<Number> Me;

    //DenseMatrix<Number> Me_p;

    // A  Gauss quadrature rule for numerical integration.
    QGauss qrule (dim, fe_type.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);

    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> > & phi = fe->get_phi();

    const DofMap & dof_map = _problem.getNonlinearSystemBase().dofMap();

    std::vector<dof_id_type> dof_indices;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // first we need to manually zero the matrix
    matrix_M.zero();

    //matrix_M_p.zero();

    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;

        Elem * ele = *el;

        fe->reinit (elem);

        dof_map.dof_indices(elem, dof_indices);

        Me.resize (dof_indices.size(), dof_indices.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

            for (unsigned int i=0; i<phi.size(); i++){

                for (unsigned int j=0; j<phi.size(); j++){

                    Me(i,j) += JxW[qp] * phi[i][qp] * phi[j][qp];
                }
            }
        }


        dof_map.constrain_element_matrix(Me,dof_indices,false);

        matrix_M.add_matrix (Me, dof_indices);



    }



    matrix_M.close();


    //matrix_M_p.close();

    utopia::SizeType nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()); //dof_map.get_n_nz();

    utopia::USparseMatrix mass_temp = utopia::local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

    utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M), mass_temp);

    // utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M_p), mass_mp);

    mass_matrix =  mass_temp;

    lumped_mass_matrix = diag(sum(mass_matrix,1));


    _console << "Assemble_Mass_matrix() end "  << std::endl;

    
}


void
FractureAppConforming::assemble_poro_mass_matrix(FEProblemBase & _problem, utopia::USparseMatrix &mass_matrixp, utopia::USparseMatrix &lumped_mass_matrixp){
    
    _console << "Assemble_Poro_Mass_matrix() begin "  << std::endl;

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _problem.es().get_mesh();

    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to our system.
//    _problem.es().add_system<LinearImplicitSystem>("aux").add_variable("var",FIRST);
//
//    _problem.es().reinit();
//
//    // Get a reference to our system.
//    LinearImplicitSystem & _system = _problem.es().get_system<LinearImplicitSystem>("aux");
    
    auto & _system = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _system.get_dof_map().variable_type(0);

    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

    SparseMatrix<Number> & matrix_M_p = *_system.matrix;

    //_console << "is matrix closed: " << matrix_A.closed() << std::endl;

    // The element mass matrix.
    DenseMatrix<Number> Me_p;

    // A  Gauss quadrature rule for numerical integration.
    QGauss qrule (dim, fe_type.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (&qrule);

    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> & JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real> > & phi = fe->get_phi();

    const DofMap & dof_map = _system.get_dof_map();

    std::vector<dof_id_type> dof_indices;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();

    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // first we need to manually zero the matrix
    matrix_M_p.zero();

    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;

        Elem * ele = *el;

        fe->reinit (elem);

        dof_map.dof_indices(elem, dof_indices);

        Me_p.resize (dof_indices.size(), dof_indices.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

            for (unsigned int i=0; i<phi.size(); i++){

                for (unsigned int j=0; j<phi.size(); j++){

                    // Me_p(i,j) += _poro[qp] * JxW[qp] * phi[i][qp] * phi[j][qp];
                    if (ele->subdomain_id()==2) Me_p(i,j) += 0.40 * JxW[qp] * phi[i][qp] * phi[j][qp];
                    if (ele->subdomain_id()==4) Me_p(i,j) += 0.20 * JxW[qp] * phi[i][qp] * phi[j][qp];
                    if (ele->subdomain_id()==5) Me_p(i,j) += 0.20 * JxW[qp] * phi[i][qp] * phi[j][qp];
                    if (ele->subdomain_id()==6) Me_p(i,j) += 0.20 * JxW[qp] * phi[i][qp] * phi[j][qp];
                    if (ele->subdomain_id()==7) Me_p(i,j) += 0.25 * JxW[qp] * phi[i][qp] * phi[j][qp];
                }
            }
        }


        //std::cout<<"Me_before"<<Me_p<<std::endl;

        dof_map.constrain_element_matrix(Me_p,dof_indices, true);

        matrix_M_p.add_matrix (Me_p, dof_indices);

        //std::cout<<"Me_after"<<Me<<std::endl;


    }



    //matrix_M.close();


    matrix_M_p.close();

    utopia::SizeType nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()); //dof_map.get_n_nz();

    utopia::USparseMatrix mass_temp = utopia::local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

    utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M_p), mass_temp);

    mass_matrixp = mass_temp;

    lumped_mass_matrixp = diag(sum(mass_matrixp,1));


    _console << "Assemble_Poro_Mass_matrix() end "  << std::endl;
    
    
}

void
FractureAppConforming::solve_stabilize_monolithic(){
    
    using namespace utopia;
    
    _console << "transport_monolithic_stabilized"  << std::endl;
    
    FEProblemBase & _m_problem = _fe_problem;
    
    auto &_m_sys = _m_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_m = _m_problem.getNonlinearSystemBase();
    

    if  (_fe_problem.timeStep()==1)
    {
    
        //libMesh::PetscMatrix<libMesh::Number> *petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_m_sys.matrix);
        
        PetscErrorCode ierr;
        
        KSP ksp;

        PetscNonlinearSolver<libMesh::Number> * petsc_solver =
        dynamic_cast<PetscNonlinearSolver<libMesh::Number> *>((_nl_m).nonlinearSolver());
        
        SNES snes = petsc_solver->snes();
        
        ierr = SNESSetFromOptions(snes);
        
        ierr = SNESGetKSP(snes,&ksp);
        
        Mat Hmat,Pmat;
        
        ierr = KSPGetOperators(ksp,&Hmat,&Pmat);

        //_m_problem.computeJacobianSys(_m_sys, *_nl_m.currentSolution(), *petsc_mat_m);

        utopia::convert(Hmat, A_m_t);
        
        utopia::disp(local_size(A_m_t).get(0));
    
        std::cout<<"inizio"<<std::endl;

//       _m_problem.computeResidualSys(_m_sys, *_nl_m.currentSolution(), _nl_m.RHS());
//
//       utopia::convert(ksp->vec_rhs, rhs_m_t);
//
        rhs_m_c = rhs_m_t;
//
//        disp(rhs_m_c);

        Real dt = static_cast<Transient*>(_m_problem.getMooseApp().getExecutioner())->getDT();

        double inv_dt = 1.0/dt;

        USparseMatrix D_stab = A_m_t;

        D_stab*=0;

        stabilize_A_matrix(_m_problem, D_stab);
    
    

        utopia::USparseMatrix A_stab = D_stab + A_m_t;
    
        std::cout<<"inizio2"<<std::endl;

        USparseMatrix A_m_tot =  1.0 *  A_stab + mass_lumped_p * inv_dt;

//        constraint_mat(rhs_m_c, A_m_tot, _constraint_m);
    
        const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).setMatrix() =  std::make_shared<USparseMatrix>(A_m_tot);
    
        std::cout<<"inizio3"<<std::endl;
        
        utopia::disp(A_m_t);
    
        std::cout<<"fine"<<std::endl;
        
        //exit(1);
    
        const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).setRHS() =  std::make_shared<UVector>(rhs_m_c);
        
        
//        auto op = std::make_shared<utopia::Factorization<utopia::USparseMatrix, utopia::UVector> >(MATSOLVERMUMPS,PCLU);

   
    }
    
    
    if (_fe_problem.timeStep()>1)
    {
        
        auto sol_m_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_m_sys.old_local_solution.get())->vec();

        utopia::convert(sol_m_old, c_m_old);

        Real dt = static_cast<Transient*>(_fe_problem.getMooseApp().getExecutioner())->getDT();

        double inv_dt = 1.0/dt;

        UVector  mass_c_m_old;

        utopia::USparseMatrix mass_lumped_pc =  mass_lumped_p;

        //constraint_mat(rhs_m_c,  mass_lumped_pc, _constraint_m);

        mass_c_m_old = mass_lumped_pc * c_m_old;

        mass_c_m_old_dot = mass_c_m_old * inv_dt;

        rhs_m_t =  - 1.0 * mass_c_m_old_dot;
        
        utopia::UVector ma= utopia::local_zeros(local_size(rhs_m_t));

        zero_constraint_vec(rhs_m_c,  rhs_m_t, _constraint_m);
        
        const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).setRHS() =  std::make_shared<UVector>(rhs_m_t);
        
//        auto op = std::make_shared<utopia::Factorization<utopia::USparseMatrix, utopia::UVector> >(MATSOLVERMUMPS,PCLU);
        
//        auto  A_m_s = const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).getMatrix();
        
//        c_m =  utopia::local_zeros(local_size(rhs_m_t));
//
//        op->update(make_ref(*A_m_s));
//
//        op->apply(rhs_m_t, c_m);
//
//        _console << "CIAO NO"  << std::endl;
//
//        CopyMatrixSolution(c_m);


    }
    
}

void
FractureAppConforming::stabilize_A_matrix(FEProblemBase & _problem, utopia::USparseMatrix &S_matrix)
{
    
    
    _console << "Stabilize A matrix:: begin  "  << std::endl;
    
    auto &_sys = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl = _problem.getNonlinearSystemBase();
    
    libMesh::PetscMatrix<libMesh::Number> *petsc_mat = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_sys.matrix);
    
    _problem.computeJacobianSys(_sys, *_nl.currentSolution(), *petsc_mat);
    
    utopia::USparseMatrix A_0, A_0_t;
    
    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat).mat(), A_0);
    
    A_0_t = utopia::transpose(A_0);
    
    
    
    
    {
        utopia::Read<utopia::USparseMatrix>  r_s(A_0), r_s_t(A_0_t);
        utopia::Write<utopia::USparseMatrix> w_s(S_matrix);
        utopia::each_read(A_0, [&](const utopia::SizeType i, const utopia::SizeType j, double value){
            if(i!=j)
            {
                double value_1 = 1.0 * value;
                
                //utopia::disp(value_1);
                
                double value_2 = 1.0 * A_0_t.get(i,j);
                
                //utopia::disp(value_2);
                
                if(value_1 > value_2){
                    
                    double max = - 1.0 * std::max(0.0, value_1);
                    
                    S_matrix.set(i,j, max);
                }
                
                else{
                    
                    double max = - 1.0 * std::max(0.0, value_2);
                    
                    S_matrix.set(i,j, max);
                }
            }
            else{
                
                S_matrix.set(i,i,0);
            }
        });
    }
    
    
    
    
    utopia::UVector diag_elem = -1.0 * sum(S_matrix,1);
    
    utopia::USparseMatrix S_diag=diag(diag_elem);
    
    // S_matrix += transpose(S_matrix);
    
    // S_matrix *=0.5;
    
    S_matrix+=S_diag;
    
    _console << "Stabilize A matrix:: end  "  << std::endl;
    
    // utopia::write("S.m", S_matrix);
    
    
}

void
FractureAppConforming::CopyMatrixSolution(utopia::UVector _sol_m)
{
    
    FEProblemBase & _problem_m = _fe_problem;
    
    MooseVariable & _var_m = _problem_m.getStandardVariable(0, _m_var_name);
    
    System & _sys_m = _var_m.sys().system();
    
    MeshBase *_mesh_m = &_problem_m.mesh().getMesh();
    
    NumericVector<Number> * _solution_m = _sys_m.solution.get();
    
    PetscInt       rstart,rend;
    
    PetscScalar    tmp_sol;
    
    VecGetOwnershipRange(utopia::raw_type(_sol_m),&rstart,&rend);
    
    _sol_m = _sol_m * (-1);
    
    utopia::Read<VecT> w_d(_sol_m);
    
    
    {
        MeshBase::const_node_iterator it = _mesh_m->local_nodes_begin();
        const MeshBase::const_node_iterator end_it = _mesh_m->local_nodes_end();
        for ( ; it != end_it; ++it)
        {
            const Node * node = *it;
            
            for (unsigned int comp = 0;comp < node->n_comp(_sys_m.number(), _var_m.number()); comp++)
            {
                const dof_id_type from_index = node->dof_number(_sys_m.number(), _var_m.number(), comp);
                
                _solution_m->set(from_index, _sol_m.get(from_index));
                
            }
        }
    }
    
    {
        MeshBase::const_element_iterator it = _mesh_m->active_local_elements_begin();
        const MeshBase::const_element_iterator end_it = _mesh_m->active_local_elements_end();
        for ( ; it != end_it; ++it)
        {
            const Elem * elem = *it;
            for (unsigned int comp = 0; comp < elem->n_comp(_sys_m.number(), _var_m.number()); comp++)
            {
                const dof_id_type from_index = elem->dof_number(_sys_m.number(), _var_m.number(), comp);
                
                _solution_m->set(from_index, _sol_m.get(from_index));
            }
        }
    }
    
    
    _solution_m->close();
    _sys_m.update();
    
    ExodusII_IO (*_mesh_m).write_equation_systems("matrix_c.e", _var_m.sys().system().get_equation_systems());
    
    
}
