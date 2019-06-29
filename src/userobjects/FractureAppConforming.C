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
#include "utopia_Socket.hpp"
#include "utopia_FractureFlowUtils.hpp"
#include "utopia_TransferUtils.hpp"
#include "utopia_TransferAssembler.hpp"
#include "FractureAppConforming.h"
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
#include "utopia_MeshTransferOperator.hpp"
#include "MultiApp.h"
#include "utopia_libmesh_FunctionSpace.hpp"
#include <queue>
#include "Transient.h"
#include <iostream>
#include "TransientMultiApp.h"
#include "utopia_assemble_contact.hpp"
#include "libmesh/boundary_info.h"
#include "libmesh/boundary_mesh.h"
#include "Porosity.h"

#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"

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
    params.addRequiredParam<VariableName>("fracture_variable",
                                          "The variable to transfer from.");
    params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    params.addParam<std::string>("operator_type","opeartor_type" /*INTERPOLATION| L2_PROJECTION| PSEUDO_L2_PROJECTION | APPROX_L2_PROJECTION*/);
    params.addParam<bool>("solve_cg","use cg and solve for Schur complement");
    params.addParam<bool>("solve_mg","solve the master problem with MG");
    params.addParam<bool>("pressure","false","put true if you solve pressure system");
    params.addParam<bool>("transport","false","put true if you transport system");
    params.addParam<bool>("boundary","false","put true for surface to volume coupling");
    params.addParam<bool>("constraint_m","false","put true if matrix has Dirichlet BC");
    params.addParam<bool>("constraint_f","false","put true if fibres has Dirichlet BC");
    params.addParam<int>("id_slave","The boundary ID associated with the master side");
    params.addParam<double>("porosity_m","posrosity matrix");
    params.addParam<double>("porosity_f","posrosity fracture");
    return params;
    
}

FractureAppConforming::FractureAppConforming(const InputParameters & parameters):
GeneralUserObject(parameters),
_poro(getMaterialProperty<Real>("Porosity")),
_slave_id(parameters.get<int>("id_slave")),
_f_var_name(getParam<VariableName>("fracture_variable")),
_m_var_name(getParam<VariableName>("matrix_variable")),
// _operator_storage(getUserObject<StoreTransferOperators>("operator_userobject")),
_operator_type(getParam<std::string>("operator_type")),
_solve_cg(getParam<bool>("solve_cg")),
_solve_mg(getParam<bool>("solve_mg")),
_pressure(getParam<bool>("pressure")),
_transport(getParam<bool>("transport")),
_boundary(getParam<bool>("boundary")),
_constraint_m(getParam<bool>("constraint_m")),
_constraint_f(getParam<bool>("constraint_f")),
_porosity_m(getParam<double>("porosity_m")),
_porosity_f(getParam<double>("porosity_f"))

{

}




void
FractureAppConforming::initialize()
{
    _console << "Initial Setup of Fracture App " << std::endl;
    
    if (_transport){

        assemble_mass_matrix(_porosity_m, _fe_problem,  mass_m, mass_lumped);

        assemble_poro_mass_matrix(_porosity_m, _fe_problem,  mass_mp, mass_lumpedp);
    }

    if(!_fe_problem.hasUserObject(_userobject_name_1)){

            _console << "Add UserObject" << std::endl;

            std::string class_name = "StoreTransferOperators";
            auto params = _fe_problem.getMooseApp().getFactory().getValidParams(class_name);
            params.set<bool>("use_displaced_mesh") = false;
            params.set<ExecFlagEnum>("execute_on") = "initial";

            _fe_problem.addUserObject("StoreTransferOperators", _userobject_name_1, params);
            _fe_problem.addUserObject("StoreTransferOperators", _userobject_name_2, params);
            _fe_problem.addUserObject("StoreTransferOperators", _userobject_name_3, params);
            _fe_problem.addUserObject("StoreTransferOperators", _userobject_name_4, params);


    }


}


bool
FractureAppConforming::solve(){
    
    if (_solve_cg){
        ok = solve_cg_dual();
    }
    else{
        ok = solve_monolithic();
    }
    
      return ok;
    
}


void
FractureAppConforming::execute()
{
   if(_pressure) 
   {
    solve(); 
    
    _console << "Solve with "  << _operator_type << std::endl;
        
    if (ok){
        CopyMatrixSolution(x_m);
        CopyFractureSolution(x_f);
        }
    }

    if(_transport){

        solve_stabilize_monolithic();
    }
}


bool FractureAppConforming::solve_monolithic()
{
    _console << "Solve_monolithic()"  << std::endl;
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
    
    FEProblemBase & _f_problem =_multi_app.problemBase();
    
    FEProblemBase & _m_problem = _multi_app.appProblemBase(0);
    
    MeshBase *_m_mesh = &_m_problem.mesh().getMesh();
    
    MeshBase *_f_mesh = &_f_problem.mesh().getMesh();
    
    MooseVariable & _m_var = _m_problem.getStandardVariable(0, _m_var_name);
    
    MooseVariable & _f_var   = _f_problem.getStandardVariable(0, _f_var_name);
    
    auto V_m = utopia::LibMeshFunctionSpace(utopia::make_ref(_m_var.sys().system().get_equation_systems()),_m_var.sys().system().number(), _m_var.number());
    
    auto V_f = utopia::LibMeshFunctionSpace(utopia::make_ref(_f_var.sys().system().get_equation_systems()),_f_var.sys().system().number(), _f_var.number());
    
    auto &_f_sys = _f_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_f = _f_problem.getNonlinearSystemBase();
    
    // Pointer to underlying PetscMatrix type
    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_f = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_f_sys.matrix);
    
    _f_problem.computeResidualSys(_f_sys, *_nl_f.currentSolution(), _nl_f.RHS());
    
    _f_problem.computeJacobianSys(_f_sys, *_nl_f.currentSolution(), *petsc_mat_f);

    auto &_m_sys = _m_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_m = _m_problem.getNonlinearSystemBase();


    // Pointer to underlying PetscMatrix type
    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_m_sys.matrix);
    
    _m_problem.computeResidualSys(_m_sys, *_nl_m.currentSolution(), _nl_m.RHS());
    
    _m_problem.computeJacobianSys(_m_sys, *_nl_m.currentSolution(), *petsc_mat_m);
    
    utopia::UVector rhs_f;

    utopia::USparseMatrix A_f;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_f.RHS()), rhs_f);

    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_f).mat(), A_f);
    
    utopia::UVector rhs_m;

    utopia::USparseMatrix A_m;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_m.RHS()), rhs_m);

    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_m).mat(), A_m);

    _console << "Solve_monolithic():: START SOLVING"  << std::endl;


    // utopia::disp(B.size());

    // utopia::disp(D.size());

    // utopia::disp(B_t.size());

    // utopia::disp(D_t.size());

    // utopia::disp(A_m.size());

    // utopia::disp(A_f.size());

   

    utopia::USparseMatrix A = utopia::Blocks<utopia::USparseMatrix>(3, 3,
                                            {
                                                utopia::make_ref(A_m), nullptr, utopia::make_ref(B_t),
                                                nullptr, utopia::make_ref(A_f), utopia::make_ref(D_t),
                                                utopia::make_ref(B), utopia::make_ref(D), nullptr
                                            });


    utopia::UVector z = utopia::local_zeros(local_size(D).get(0));

    utopia::UVector rhs = utopia::blocks(rhs_m, rhs_f, z);
    
    x_m  = utopia::local_zeros(local_size(rhs_m));

    x_f  = utopia::local_zeros(local_size(rhs_f));

    lagr = utopia::local_zeros(local_size(D).get(0));
    
    utopia::UVector sol = blocks(x_m, x_f, lagr);
    
    utopia::Factorization<utopia::USparseMatrix, utopia::UVector> op;
   
    op.update(make_ref(A));
    
    bool ok = op.apply(rhs, sol);

    _console << "Solve_monolithic():: FINISH SOLVING"  << std::endl;
    
    utopia::undo_blocks(sol, x_m, x_f, lagr);
    
    return ok;
}
bool FractureAppConforming::solve_cg_dual(){
    
    _console << "Solve_cg_dual()  "  << std::endl;
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
    
    FEProblemBase & _f_problem =_multi_app.problemBase();
    
    FEProblemBase & _m_problem = _multi_app.appProblemBase(0);
    
    MeshBase *_f_mesh = &_f_problem.mesh().getMesh();
    
    MeshBase * _m_mesh = &_m_problem.mesh().getMesh();
    
    MooseVariable & _m_var = _m_problem.getStandardVariable(0, _m_var_name);
    
    MooseVariable & _f_var   = _f_problem.getStandardVariable(0, _f_var_name);
    
    auto V_m = utopia::LibMeshFunctionSpace(utopia::make_ref(_m_var.sys().system().get_equation_systems()),_m_var.sys().system().number(), _m_var.number());
    
    auto V_f = utopia::LibMeshFunctionSpace(utopia::make_ref(_f_var.sys().system().get_equation_systems()),_f_var.sys().system().number(), _f_var.number());
    
    utopia::SPBlockConjugateGradient<utopia::USparseMatrix, utopia::UVector> solver;
    
    solver.verbose(true);
    
    solver.max_it(2000);
    
    solver.atol(1e-14);
    
    solver.use_simple_preconditioner();
    
    int mg_sweeps = 1;
    int mg_levels = 5;
    
    if(_solve_mg) {
        auto mg = make_mg_solver(V_m, mg_levels);
        solver.set_master_solver(mg);
        
        solver.set_master_sweeps(mg_sweeps);
        solver.set_master_max_it(mg->max_it());
    }
    
    //to_problem.es().print_info();
    
    auto &_f_sys = _f_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_f = _f_problem.getNonlinearSystemBase();

    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_f = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_f_sys.matrix);
    
    _f_problem.computeResidualSys(_f_sys, *_nl_f.currentSolution(), _nl_f.RHS());
    
    _f_problem.computeJacobianSys(_f_sys, *_nl_f.currentSolution(), *petsc_mat_f);
    
    auto &_m_sys = _m_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_m = _m_problem.getNonlinearSystemBase();
    
    // Pointer to underlying PetscMatrix type
    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_m_sys.matrix);
    
    _m_problem.computeResidualSys(_m_sys, *_nl_m.currentSolution(), _nl_m.RHS());
    
    _m_problem.computeJacobianSys(_m_sys, *_nl_m.currentSolution(), *petsc_mat_m);
    
    utopia::UVector rhs_f;
    utopia::USparseMatrix A_f;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_f.RHS()), rhs_f);
    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_f).mat(), A_f);
    //disp(mat_to);
    
    utopia::UVector rhs_m;
    utopia::USparseMatrix A_m;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_m.RHS()), rhs_m);
    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_m).mat(), A_m);
    //disp(mat_from);
    
    solver.update(utopia::make_ref(A_m),
                  utopia::make_ref(A_f),
                  utopia::make_ref(B),
                  utopia::make_ref(D),
                  utopia::make_ref(B_t),
                  utopia::make_ref(D_t)
                  );
   
//  rhs_f*=-1.0;
//  rhs_m*=-1.0;
    
  lagr = utopia::local_zeros(local_size(D).get(0));
  bool ok = solver.apply(rhs_m, rhs_f, x_m, x_f, lagr);
    
  return ok;
    
    
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

    // if(_pressure)
    //    ExodusII_IO (*_mesh_m).write_equation_systems("matrix_p.e", _var_m.sys().system().get_equation_systems());
    // else
    //    ExodusII_IO (*_mesh_m).write_equation_systems("matrix_c.e", _var_m.sys().system().get_equation_systems());
    
    
}

void
FractureAppConforming::CopyFractureSolution(utopia::UVector _sol_f)
{
    
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);

    FEProblemBase & _problem_f =_multi_app.problemBase();

    MooseVariable & _var_f   = _problem_f.getStandardVariable(0, _f_var_name);
    
    System & _sys_f = _var_f.sys().system();

    MeshBase * _mesh_f = &_problem_f.mesh().getMesh();
    
    NumericVector<Number> * _solution_f = _sys_f.solution.get();
    
    PetscInt       rstart,rend;
    
    PetscScalar    tmp_sol;
    
    VecGetOwnershipRange(utopia::raw_type(_sol_f),&rstart,&rend);
    
    _sol_f=_sol_f*(-1);
    
   // disp(to_sol);
    
    utopia::Read<VecT> w_d(_sol_f);
    
    
    {
        MeshBase::const_node_iterator it = _mesh_f->local_nodes_begin();
        const MeshBase::const_node_iterator end_it = _mesh_f->local_nodes_end();
        for ( ; it != end_it; ++it)
        {
            const Node * node = *it;
            
            for (unsigned int comp = 0;comp < node->n_comp(_sys_f.number(), _var_f.number()); comp++)
            {
                const dof_id_type to_index = node->dof_number(_sys_f.number(), _var_f.number(), comp);
                
                _solution_f->set(to_index,  _sol_f.get(to_index));
                
            }
        }
    }
    
    {
        MeshBase::const_element_iterator it = _mesh_f->active_local_elements_begin();
        const MeshBase::const_element_iterator end_it = _mesh_f->active_local_elements_end();
        for ( ; it != end_it; ++it)
        {
            const Elem * elem = *it;
            for (unsigned int comp = 0; comp < elem->n_comp(_sys_f.number(), _var_f.number()); comp++)
            {
                const dof_id_type to_index = elem->dof_number(_sys_f.number(), _var_f.number(), comp);
                _solution_f->set(to_index, _sol_f.get(to_index));
            }
        }
    }
    
    
    _solution_f->close();
    _sys_f.update();
    
    if(_pressure)
       ExodusII_IO (*_mesh_f).write_equation_systems("fracture_p.e", _var_f.sys().system().get_equation_systems());
    else
       ExodusII_IO (*_mesh_f).write_equation_systems("fracture_c.e", _var_f.sys().system().get_equation_systems());

    
}



void
FractureAppConforming::constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
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
FractureAppConforming::unconstraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
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
FractureAppConforming::one_constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
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
FractureAppConforming::constraint_concentration_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints)
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
FractureAppConforming::assemble_mass_matrix(double porosity, FEProblemBase & _problem, utopia::USparseMatrix &mass_matrix, utopia::USparseMatrix &lumped_mass_matrix){
    
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

        // Me_p.resize (dof_indices.size(), dof_indices.size());
        
        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

            for (unsigned int i=0; i<phi.size(); i++){

                for (unsigned int j=0; j<phi.size(); j++){

                        Me(i,j) += JxW[qp] * phi[i][qp] * phi[j][qp];
            }
        }
    }

        
        dof_map.constrain_element_matrix(Me,dof_indices,false);

        matrix_M.add_matrix (Me, dof_indices);
        //matrix_M_p.add_matrix (Me_p, dof_indices);

        
    }


    
    matrix_M.close();

    utopia::SizeType nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()); //dof_map.get_n_nz();
        
    utopia::USparseMatrix mass_temp = utopia::local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

    utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M), mass_temp);

    mass_matrix =  mass_temp;

    lumped_mass_matrix = diag(sum(mass_matrix,1));


    _console << "Assemble_Mass_matrix() end "  << std::endl;

    
}


void
FractureAppConforming::assemble_poro_mass_matrix(double porosity, FEProblemBase & _problem, utopia::USparseMatrix &mass_matrixp, utopia::USparseMatrix &lumped_mass_matrixp){
    
   _console << "Assemble_Poro_Mass_matrix() begin "  << std::endl;
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _problem.es().get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();    
    
    // Get a reference to our system.
    _problem.es().add_system<LinearImplicitSystem>("aux").add_variable("var",FIRST);
    
    _problem.es().reinit();
    
    // Get a reference to our system.
    LinearImplicitSystem & _system = _problem.es().get_system<LinearImplicitSystem>("aux");
    
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

        
        dof_map.constrain_element_matrix(Me_p,dof_indices,false);

        matrix_M_p.add_matrix (Me_p, dof_indices);

        
    }


    
    //matrix_M.close();


    matrix_M_p.close();

    utopia::SizeType nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()); //dof_map.get_n_nz();
        
    utopia::USparseMatrix mass_temp = utopia::local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

    utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M_p), mass_temp);

    // utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M_p), mass_mp);

    mass_matrixp = mass_temp;

    lumped_mass_matrixp = diag(sum(mass_matrixp,1));


    _console << "Assemble_Poro_Mass_matrix() end "  << std::endl;

    
}

void
FractureAppConforming::solve_stabilize_monolithic(){

  using namespace utopia;
    
    _console << "transport_monolithic_stabilized"  << std::endl;

    //MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
        
    FEProblemBase & _f_problem = _fe_problem;
    
    FEProblemBase & _m_problem = _fe_problem;

    auto &_m_sys = _m_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    NonlinearSystemBase & _nl_m = _m_problem.getNonlinearSystemBase();
    
    auto &_f_sys = _f_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    NonlinearSystemBase & _nl_f = _f_problem.getNonlinearSystemBase();



    if  (_fe_problem.timeStep()==1)
     {
              
        
        /// Pointer to underlying PetscMatrix type
        libMesh::PetscMatrix<libMesh::Number> *petsc_mat_f = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_f_sys.matrix);
        
        _f_problem.computeJacobianSys(_f_sys, *_nl_f.currentSolution(), *petsc_mat_f);
        
        // Pointer to underlying PetscMatrix type
        libMesh::PetscMatrix<libMesh::Number> *petsc_mat_m = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_m_sys.matrix);

        
        _m_problem.computeJacobianSys(_m_sys, *_nl_m.currentSolution(), *petsc_mat_m);
        
        
        utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_f).mat(), A_f_t);

        utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_m).mat(), A_m_t);

       
  
        _f_problem.computeResidualSys(_f_sys, *_nl_f.currentSolution(), _nl_f.RHS());

        utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_f.RHS()), rhs_f_t);

        utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_m.RHS()), rhs_m_t);


        rhs_m_c = rhs_m_t; 

        rhs_f_c = rhs_f_t; 


        Real dt = static_cast<Transient*>(_m_problem.getMooseApp().getExecutioner())->getDT();

        double inv_dt = 1.0/dt;

        _A_t_store   = std::make_shared<USparseMatrix>(A_m_t);



        auto op = std::make_shared<utopia::Factorization<utopia::USparseMatrix, utopia::UVector> >(MATSOLVERMUMPS,PCLU);

        const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_4)).getVoidPointer() = op;

        c_m  = utopia::local_zeros(local_size(rhs_m_t));


        
        USparseMatrix A_m_t_stab = A_m_t;

        A_m_t_stab*=0;

        _console << "solve_transport_monolithic_one::BEGIN "  << std::endl;

        stabilize_A_matrix(_fe_problem, A_m_t, A_m_t_stab);

        utopia::USparseMatrix A_stab = A_m_t_stab + A_m_t;

        USparseMatrix A_m_tot =  1.0 *  A_stab + mass_lumpedp * inv_dt;

        constraint_concentration_mat(rhs_m_c, A_m_tot, _constraint_m);

        op->update(make_ref(A_m_tot));
        
        op->apply(rhs_m_t, c_m);

        

        // utopia::UVector rhs_dot;

        // USparseMatrix mass_m_c = mass_m;

        // USparseMatrix mass_lumped_c = mass_lumped;

        // constraint_concentration_mat(rhs_m_c, mass_m_c, _constraint_m);

        // constraint_concentration_mat(rhs_m_c, mass_lumped_c, _constraint_m);

        // UVector diag_elem_c = 1./sum((mass_lumped_c),1);

        // utopia::USparseMatrix inv_mass_lumped_mat_c = diag(diag_elem_c);

        // UVector diag_elem = 1./sum((mass_lumped_c),1);

        // utopia::USparseMatrix inv_mass_lumped_mat = diag(diag_elem);


        // c_dot = utopia::local_zeros(utopia::local_size(c_m));

        // USparseMatrix A_stab_c = A_stab;

        // constraint_concentration_mat(rhs_m_c, A_stab_c, _constraint_m);

        // c_dot = 1.0 * inv_mass_lumped_mat * A_stab * c_m;

        // //c_dot *= -1.0;


        // USparseMatrix A_m_t_stab_c = A_m_t_stab;

        // constraint_concentration_mat(rhs_m_c, A_m_t_stab_c, _constraint_m);

        // constraint_concentration_mat(rhs_m_c, mass_m_c, _constraint_m);

        // utopia::UVector _f_s = local_zeros(local_size(c_dot));

        // stabilize_coeffiecient(_fe_problem, c_m, c_dot, A_m_t_stab, mass_m, _f_s);


        // utopia::UVector c_m_tot = c_m + 1.0 * dt * inv_mass_lumped_mat_c * _f_s;


        CopyMatrixSolution(c_m);

        const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_3)).setTransferOperator() =  std::make_shared<USparseMatrix>(A_m_tot);

        const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).setTransferOperator() =  std::make_shared<USparseMatrix>(A_m_t_stab);

        const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_2)).setTransferOperator() =  std::make_shared<USparseMatrix>(A_stab);

        _console << "solve_transport_monolithic_one:: END "  << std::endl;
        
    
    }


    if (_fe_problem.timeStep()>0)
    {

        auto sol_m_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_m_sys.old_local_solution.get())->vec();

        utopia::convert(sol_m_old, c_m_old);

        Real dt = static_cast<Transient*>(_fe_problem.getMooseApp().getExecutioner())->getDT();

        double inv_dt = 1.0/dt;

        UVector  mass_c_m_old;

        utopia::USparseMatrix mass_lumpedpc =  mass_lumpedp;

        constraint_concentration_mat(rhs_m_c,  mass_lumpedpc, _constraint_m);

        mass_c_m_old = mass_lumpedpc * c_m_old;

        mass_c_m_old_dot = mass_c_m_old * inv_dt;

        rhs_m_t = - 1.0 * mass_c_m_old_dot;
  
        constraint_concentration_vec(rhs_m_c,  rhs_m_t, _constraint_m);

        auto op = std::static_pointer_cast< Factorization<USparseMatrix, UVector> >(const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_4)).getVoidPointer());



        //utopia::write("A.m", A_stab);
        
        //USparseMatrix A_m_tot =  1.0 * A_stab  + mass_lumpedp * inv_dt;

        USparseMatrix A_m_tot = *const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_3)).getTransferOperator();

        USparseMatrix A_m_t_stab = *const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_1)).getTransferOperator();

        USparseMatrix A_stab = *const_cast<StoreTransferOperators&>(_fe_problem.getUserObject<StoreTransferOperators>(_userobject_name_2)).getTransferOperator();

         _console << "solve_transport_monolithic_one:: Solve  "  << std::endl;

        op->update(make_ref(A_m_tot));
        
        op->apply(rhs_m_t, c_m);

        utopia::UVector rhs_dot;

        USparseMatrix mass_m_c = mass_m;

        USparseMatrix mass_lumped_c = mass_lumped;

        // constraint_concentration_mat(rhs_m_c, mass_lumped, _constraint_m);

        constraint_concentration_mat(rhs_m_c, mass_m_c, _constraint_m);

        constraint_concentration_mat(rhs_m_c, mass_lumped_c, _constraint_m);

        UVector diag_elem_c = 1./sum((mass_lumped_c),1);

        utopia::USparseMatrix inv_mass_lumped_mat_c = diag(diag_elem_c);

        UVector diag_elem = 1./sum((mass_lumped),1);

        // disp(diag_elem_c);

        utopia::USparseMatrix inv_mass_lumped_mat = diag(diag_elem);

        // rhs_dot =  A_stab * c_m;

        c_dot = utopia::local_zeros(utopia::local_size(c_m));

        USparseMatrix A_stab_c = A_stab;

        //constraint_concentration_mat(rhs_m_c, A_stab_c, _constraint_m);

        c_dot = 1.0 * inv_mass_lumped_mat_c * A_stab_c * c_m;

        // utopia::UVector rhs_l = (*_A_t_store) * c_m;

        //utopia::UVector rhs_dot_c = local_zeros(local_size(rhs_dot));

        //constraint_concentration_vec(rhs_m_c,  c_dot, _constraint_m);

        // op->update(make_ref(mass_m_c));

        // op->apply(rhs_l, c_dot);

        c_dot*=1.0;


        // utopia::UVector c_dot_old = - 1.0 * inv_mass_lumped_mat_c * A_stab * c_m_old;

        // utopia::UVector c_m_tilde = c_m_old + dt/2.0 * c_dot_old;

        //constraint_concentration_vec(c_dot,  _f_con, _constraint_m);

        //utopia::UVector _f = inv_mass_lumped_mat * ((mass_lumped - mass_m) * _f_con  + 1.0 * A_m_t_stab * dt * c_m);

        USparseMatrix A_m_t_stab_c = A_m_t_stab;

        constraint_concentration_mat(rhs_m_c, A_m_t_stab_c, _constraint_m);

        constraint_concentration_mat(rhs_m_c, mass_m_c, _constraint_m);

        utopia::UVector _f_s = local_zeros(local_size(c_dot));

        // unstabilize_coeffiecient(_fe_problem, c_m, c_dot, A_m_t_stab, mass_m, _f_s);

        //unconstraint_concentration_vec(rhs_m_c,  _f_s, _constraint_m);


        stabilize_coeffiecient(_fe_problem, c_m, c_dot, A_m_t_stab, mass_m, _f_s);

        //utopia::disp(_f_s);


        // utopia::UVector _f_con = local_zeros(local_size(c_m));

        // _f_con = dt * f_s + rhs_m_t;

        // unconstraint_concentration_vec(rhs_m_c,  _f_s, _constraint_m);

        // _console << "transport_monolithic:: Correction"  << std::endl;

        //_f_s = (mass_m - mass_lumped) * (c_m_old - c_m) * inv_dt + 1.0 * A_m_t_stab * c_m;

        utopia::UVector c_m_tot = c_m + 1.0 * dt * inv_mass_lumped_mat_c * _f_s;

        //utopia::disp(corr);
        
        CopyMatrixSolution(c_m_tot);

        // rhs_m_t+= _f_s;

        //utopia::disp(_f_con);

        // c_m_tot = c_m - 1.0 * _f_con * dt;

        // op->update(make_ref(A_m_tot));
        
        // op->apply(_f_con, c_m);

        // CopyMatrixSolution(c_m);

    }

}



void
FractureAppConforming::stabilize_A_matrix(FEProblemBase & _problem, utopia::USparseMatrix &A_0, utopia::USparseMatrix &S_matrix)
{ 
    

    _console << "Stabilize A matrix:: begin  "  << std::endl;

    // auto &_sys = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    // NonlinearSystemBase & _nl = _problem.getNonlinearSystemBase();

    // libMesh::PetscMatrix<libMesh::Number> *petsc_mat = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_sys.matrix);
        
    // _problem.computeJacobianSys(_sys, *_nl.currentSolution(), *petsc_mat);

    utopia::USparseMatrix A_0_t;

    //utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat).mat(), A_0);

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
FractureAppConforming::stabilize_coeffiecient(FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs)
{

    _console << "stabilize_coeffiecient::begin"  << std::endl;

    utopia::USparseMatrix _f = _D + _M; 

    _f *=0;

    auto &_sys = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    DofMap & dof_map_1 = _sys.get_dof_map();

    MooseVariable & _m_var = _fe_problem.getStandardVariable(0, _m_var_name);

    std::vector<int> send_list_v;

    FractureAppConforming::send_list(_fe_problem, _M, send_list_v);

    std::sort(send_list_v.begin(), send_list_v.end());

    send_list_v.erase( unique( send_list_v.begin(), send_list_v.end() ),  send_list_v.end());


    utopia::UVector sol_ghosted_u = utopia::ghosted(dof_map_1.n_local_dofs(), dof_map_1.n_dofs(), send_list_v);
    sol_ghosted_u = _u;


    utopia::UVector sol_ghosted_u_dot = utopia::ghosted(dof_map_1.n_local_dofs(), dof_map_1.n_dofs(), send_list_v);
    sol_ghosted_u_dot = _u_dot;

    utopia::synchronize(sol_ghosted_u);
    utopia::synchronize(sol_ghosted_u_dot);


    assert(sol_ghosted_u.implementation().has_ghosts());
    assert(sol_ghosted_u_dot.implementation().has_ghosts());

    
    {    
        utopia::Read<utopia::UVector>  r_u(_u_dot), r_v(_u), r_s_u(sol_ghosted_u), r_s_d(sol_ghosted_u_dot);

        utopia::Read<utopia::USparseMatrix>  r_s(_D), r_m(_M);

        utopia::Write<utopia::USparseMatrix> w_m(_f);

        utopia::Range rr = utopia::row_range(_M);

        std::vector<int> index_v;


        for(auto i = rr.begin(); i != rr.end(); ++i) {

            utopia::RowView<const utopia::USparseMatrix> row_view(_M, i);

            decltype(i) n_values = row_view.n_values();

            index_v.clear();

            for(auto index = 0; index < n_values; ++index) {

                const decltype(i) j = row_view.col(index);

                const auto a_ij = row_view.get(index);

                if(std::abs(a_ij) > 1.e-14) {

                    index_v.push_back(j);
                }
            }

            if (!index_v.empty()) {

                std::vector<double> values(index_v.size());

                std::vector<double> values_d(index_v.size());

                sol_ghosted_u.get(index_v, values);

                sol_ghosted_u_dot.get(index_v, values_d);

                auto v_i = _u.get(i);

                auto v_d_i = _u_dot.get(i);

                auto num_values = values.size();

                for(int k=0; k<num_values; k++){

                    int j_ind = index_v.at(k);

                    if(i != j_ind) {

                        double entry = _M.get(i,j_ind) * (v_d_i - values_d.at(k)) - _D.get(i,j_ind) * (v_i - values.at(k));

                        _f.set(i, j_ind, entry);

                        double check = ( v_i - values.at(k) ) * entry;

                        if(check==0) _f.set(i, j_ind,0.0);

                        else _f.set(i, j_ind,entry);
                    }                   
                }
            }
        }

    }


    // utopia::write("F.m",_f);

  

    utopia::UVector P_plus = utopia::local_zeros(local_size(_D).get(0));


    utopia::UVector P_minus = utopia::local_zeros(local_size(_D).get(0));


    {    

        utopia::Write<utopia::UVector> w_p(P_plus), w_m(P_minus);

        utopia::each_read(_f, [&](const utopia::SizeType i, const utopia::SizeType j, double value){
              if(j!=i)
              {
                double entry_max = std::max(0.0, value);
                
                double entry_min = std::min(0.0, value);
        
                P_minus.add(i,entry_min);

                P_plus.add(i,entry_max);
            }
        });
    } 


    // utopia::disp("P_minus");

    // utopia::disp(P_minus);

    // utopia::disp("P_plus");

    // utopia::disp(P_plus);

    // utopia::write("P_minus.m", P_minus);

    // utopia::write("P_plus.m", P_plus);
    

    utopia::UVector Q_plus = utopia::local_zeros(local_size(_f).get(0));

    utopia::UVector Q_minus = utopia::local_zeros(local_size(_f).get(0));

    // auto &_sys = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    // const DofMap & dof_map = _sys.get_dof_map();

    Real dt = static_cast<Transient*>(_problem.getMooseApp().getExecutioner())->getDT();

   {    
        
        utopia::Range rr = utopia::row_range(_M);
        
        utopia::Write<utopia::UVector> w_p(Q_plus), w_m(Q_minus);

        utopia::Read<utopia::USparseMatrix> r_M(_M);

        utopia::UVector _m = sum(_M,1);

        std::vector<int> index_w;

        utopia::Read<utopia::UVector> r_g(sol_ghosted_u), r_u(_u), r_m(_m);

        for(auto i = rr.begin(); i != rr.end(); ++i) {

            utopia::RowView<const utopia::USparseMatrix> row_view(_M, i);
            
            decltype(i) n_values = row_view.n_values();

            index_w.clear();

            for(auto index = 0; index < n_values; ++index) {

                const decltype(i) j = row_view.col(index);

                const auto a_ij = row_view.get(index);

                //tirare fuori gli indici e loopare

                if(std::abs(a_ij) > 1.e-4) {

                    index_w.push_back(j);
                }
            }


            if (index_w.size() > 0){
            
                std::vector<double> values(index_w.size());

                sol_ghosted_u.get(index_w, values);

                auto v_i = _u.get(i);

                auto m_i = _m.get(i);

                //auto num_values = values.size();

                double max = m_i/dt*(*std::max_element(std::begin(values), end(values)) - v_i);

                double min = m_i/dt*(*std::min_element(std::begin(values), end(values)) - v_i);

                Q_plus.set(i,max);

                Q_minus.set(i,min);
                
            }
        } 
    }
    


    // utopia::disp("Q_minus");

    // utopia::disp(Q_minus);

    // utopia::disp("Q_plus");

    // utopia::disp(Q_plus);

     
    
    utopia::UVector R_plus = utopia::local_zeros(local_size(Q_plus));

    utopia::UVector R_minus = utopia::local_zeros(local_size(Q_minus));

    // utopia::UVector R_plus_c = utopia::local_zeros(local_size(Q_plus));

    // utopia::UVector R_minus_c = utopia::local_zeros(local_size(Q_plus));

   {    
        

    
        utopia::Read<utopia::UVector> w_qp(Q_plus), w_qm(Q_minus), w_pp(P_plus), w_pm(P_minus);

        utopia::Write<utopia::UVector> w_Rp(R_plus), w_Rm(R_minus);

        utopia::each_read(Q_plus, [&](const utopia::SizeType i, double value){
        
        if(std::abs(P_plus.get(i))>0){    
            double value_p=value/P_plus.get(i);
            double entry_p=std::min(1.0,value_p);
            R_plus.set(i, entry_p);
        }
    
        });

        utopia::each_read(Q_minus, [&](const utopia::SizeType i, double value){
        

        if(std::abs(P_minus.get(i))>0){   
            double value_m=value/P_minus.get(i);
            double entry_m=std::min(1.0,value_m);
            R_minus.set(i, entry_m);
        }
    
    
        });



    } 
     
    one_constraint_concentration_vec(rhs_m_c,  R_minus, _constraint_m);

    one_constraint_concentration_vec(rhs_m_c,  R_plus, _constraint_m);



    // utopia::disp("R_minus");

    // utopia::disp(R_minus_c);

    // utopia::disp("R_plus");

    // utopia::disp(R_plus_c);



    // utopia::disp(R_plus);

    // utopia::write("R_minus.m", R_minus);

    // utopia::write("R_plus.m", R_plus);


   

    utopia::USparseMatrix Alpha_m = _f;
    Alpha_m*=0;



    utopia::UVector sol_ghosted_rp = utopia::ghosted(dof_map_1.n_local_dofs(), dof_map_1.n_dofs(), send_list_v);
    sol_ghosted_rp = R_plus;


    utopia::UVector sol_ghosted_rm = utopia::ghosted(dof_map_1.n_local_dofs(), dof_map_1.n_dofs(), send_list_v);
    sol_ghosted_rm = R_minus;

    utopia::synchronize(sol_ghosted_rp);

    utopia::synchronize(sol_ghosted_rm);



    {    
        utopia::Write<utopia::USparseMatrix> w_Alpha(Alpha_m);

        utopia::Read<utopia::UVector> r_Rp(R_plus), r_Rm(R_minus);

        utopia::Read<utopia::USparseMatrix> r_f(_f);

        utopia::Range rr = utopia::row_range(_f);

        std::vector<int> index_r;

        utopia::Read<utopia::UVector> r_gp(sol_ghosted_rp);

        utopia::Read<utopia::UVector> r_gm(sol_ghosted_rm);

        for(auto i = rr.begin(); i != rr.end(); ++i) {

            utopia::RowView<const utopia::USparseMatrix> row_view(_M, i);
            
            decltype(i) n_values = row_view.n_values();

            index_r.clear();

            for(auto index = 0; index < n_values; ++index) {

                const decltype(i) j = row_view.col(index);

                const auto a_ij = row_view.get(index);

                //tirare fuori gli indici e loopare

                if(std::abs(a_ij) > 0) {

                    index_r.push_back(j);
                }
            }


            if (index_r.size() > 0){
            
                std::vector<double> values_p(index_r.size());

                std::vector<double> values_m(index_r.size());

                sol_ghosted_rp.get(index_r, values_p);

                sol_ghosted_rm.get(index_r, values_m);

                int n_values = values_p.size();

                for(int k=0; k<n_values; k++){

                    int j_ind = index_r.at(k);

                    auto f_value= _f.get(i,j_ind);

                    //std::cout<<"values_f "<< f_value <<std::endl;;

                    if(f_value > 0){
                        double alpha = std::min(R_plus.get(i), values_m.at(k));

                        //std::cout<<"values_m "<<values_m.at(k)<<std::endl;

                        auto _f_bar = 1.0 * alpha * f_value;
                        Alpha_m.set(i, j_ind, _f_bar);
                    }
                    else if(f_value < 0){
                        double alpha = std::min(R_minus.get(i), values_p.at(k));
                        auto _f_bar =  1.0 * alpha * f_value;
                        Alpha_m.set(i, j_ind, _f_bar);
                    }
                }
            } 
        }

    }
    


    // {
    //     utopia::Write<utopia::USparseMatrix> w_Alpha(Alpha_m);

    //     utopia::Read<utopia::UVector> r_Rp(R_plus), r_Rm(R_minus);

    //     utopia::Read<utopia::USparseMatrix> r_f(_f);

    //     utopia::each_read(_f, [&](const utopia::SizeType i, const utopia::SizeType j, double value){

    //         if(std::abs(value) > 0){
    //             double alpha = std::min(R_plus.get(i), R_minus.get(j));
    //             auto _f_bar = 1.0 * alpha * value;
    //             Alpha_m.set(i, j, _f_bar);
    //         }
    //         else if(value < 0){
    //             double alpha = std::min(R_minus.get(i), R_plus.get(j));
    //             auto _f_bar =  1.0 * alpha * value;
    //             Alpha_m.set(i, j, _f_bar);
    //         }
    //     });
     
    // }


    // utopia::disp(Alpha_m);

    _rhs = 1.0 * utopia::sum(Alpha_m, 1);


    // utopia::write("rhs.m", _rhs);



    unconstraint_concentration_vec(rhs_m_c,  _rhs, _constraint_m);


    _console << "Stabilize A matrix:: begin  "  << std::endl;
    
   
}
  


void 
FractureAppConforming::send_list(FEProblemBase & _problem, utopia::USparseMatrix &_M, std::vector<int> &a){

std::set<int> set_list;

_console << "create list::begin"  << std::endl;

auto &_sys = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

DofMap & dof_map_1 = _sys.get_dof_map();
 

utopia::Range rr = utopia::row_range(_M);

utopia::Read<utopia::USparseMatrix>  r_m(_M);


for(auto i = rr.begin(); i != rr.end(); ++i) {

    utopia::RowView<const utopia::USparseMatrix> row_view(_M, i);
    
    decltype(i) n_values = row_view.n_values();

    for(auto index = 0; index < n_values; ++index) {

        const decltype(i) j = row_view.col(index);

        const auto a_ij = row_view.get(index);

        set_list.insert(j);
    }
}

      a.clear();   
      a.insert(a.end(), set_list.begin(), set_list.end());

     _console << "create list::end"  << std::endl;
 }
    

void 
FractureAppConforming::unstabilize_coeffiecient(FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs)
{

    _console << "ustabilize_coeffiecient::begin"  << std::endl;

    utopia::USparseMatrix _f = _D + _M; 

    _f *=0;

    // {
    //     utopia::Read<utopia::UVector>  r_u(_u_dot), r_v(_u);

    //     utopia::Read<utopia::USparseMatrix>  r_s(_D), r_m(_M);

    //     utopia::Write<utopia::USparseMatrix> w_m(_f);

    //     utopia::each_read(_M, [&](const utopia::SizeType i, const utopia::SizeType j, double value){
    //       if (i!=j){
            
    //         double entry = value * (_u_dot.get(i) - _u_dot.get(j)) - _D.get(i,j) * (_u.get(i) - _u.get(j));

    //          _f.set(i,j,entry);
    //      }
    //     });
    // }


    auto &_sys = _problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    DofMap & dof_map_1 = _sys.get_dof_map();

    MooseVariable & _m_var = _fe_problem.getStandardVariable(0, _m_var_name);

    std::vector<int> send_list_v;

    FractureAppConforming::send_list(_fe_problem, _M, send_list_v);

    std::sort(send_list_v.begin(), send_list_v.end());

    send_list_v.erase( unique( send_list_v.begin(), send_list_v.end() ),  send_list_v.end());


    utopia::UVector sol_ghosted_u = utopia::ghosted(dof_map_1.n_local_dofs(), dof_map_1.n_dofs(), send_list_v);
    sol_ghosted_u = _u;


    utopia::UVector sol_ghosted_u_dot = utopia::ghosted(dof_map_1.n_local_dofs(), dof_map_1.n_dofs(), send_list_v);
    sol_ghosted_u_dot = _u_dot;

    utopia::synchronize(sol_ghosted_u);

    utopia::synchronize(sol_ghosted_u_dot);


    assert(sol_ghosted_u.implementation().has_ghosts());

    assert(sol_ghosted_u_dot.implementation().has_ghosts());

    
    {    
        utopia::Read<utopia::UVector>  r_u(_u_dot), r_v(_u), r_s_u(sol_ghosted_u), r_s_d(sol_ghosted_u_dot);

        utopia::Read<utopia::USparseMatrix>  r_s(_D), r_m(_M);

        utopia::Write<utopia::USparseMatrix> w_m(_f);

        utopia::Range rr = utopia::row_range(_M);

        std::vector<int> index_v;


        for(auto i = rr.begin(); i != rr.end(); ++i) {

            utopia::RowView<const utopia::USparseMatrix> row_view(_M, i);

            decltype(i) n_values = row_view.n_values();

            index_v.clear();

            for(auto index = 0; index < n_values; ++index) {

                const decltype(i) j = row_view.col(index);
                const auto a_ij = row_view.get(index);

                if(std::abs(a_ij) > 1.e-14) {

                    index_v.push_back(j);
                }
            }

            if (!index_v.empty()) {

                std::vector<double> values(index_v.size());
                std::vector<double> values_d(index_v.size());

                sol_ghosted_u.get(index_v, values);
                sol_ghosted_u_dot.get(index_v, values_d);

                auto v_i = _u.get(i);

                auto v_d_i = _u_dot.get(i);

                auto num_values = values.size();

                for(int k=0; k<num_values; k++){

                    int j_ind = index_v.at(k);

                    if(i != j_ind) {

                        double entry = _M.get(i,j_ind) * (v_d_i - values_d.at(k)) - _D.get(i,j_ind) * (v_i - values.at(k));

                        double check = ( v_i - values.at(k) ) * entry;

                        if(check==0) _f.set(i, j_ind,0.0);

                        else _f.set(i, j_ind,entry);
                    }                   
                }
            }
        }

    }

        _rhs = utopia::sum(_f,1);

        _console << "ustabilize_coeffiecient::end"  << std::endl;


}





/*    UVector P_minus = sum(transpose(P_m_minus),1);

    UVector P_max = sum(transpose(P_m_max),1);*/



    





      // USparseMatrix A_m_t_stab_c = A_m_t_stab;

      //   constraint_concentration_mat(rhs_m_c, A_m_t_stab_c, _constraint_m);

      //   constraint_concentration_mat(rhs_m_c, mass_m_c, _constraint_m);

      //   utopia::UVector _f_s = local_zeros(local_size(c_dot));

      //   // unstabilize_coeffiecient(_fe_problem, c_m, c_dot, A_m_t_stab, mass_m, _f_s);

      //   //unconstraint_concentration_vec(rhs_m_c,  _f_s, _constraint_m);




      //   // unstabilize_coeffiecient(_fe_problem, c_m, c_dot, A_m_t_stab, mass_m, _f_s);

      //   unconstraint_concentration_vec(rhs_m_c, _f_s, true);

      //   //utopia::disp(_f_s);

      //   // utopia::UVector _f_con = local_zeros(local_size(c_m));

      //   // _f_con = dt * f_s + rhs_m_t;

      //   // unconstraint_concentration_vec(rhs_m_c,  _f_s, _constraint_m);

      //   // _console << "transport_monolithic:: Correction"  << std::endl;

      //   utopia::UVector c_m_tot = c_m - 1.0 * inv_mass_lumped_mat * _f_s;

      //   //utopia::disp(corr);
        
      //   CopyMatrixSolution(c_m_tot);

      //   // rhs_m_t+= _f_s;

      //   //utopia::disp(_f_con);

      //   // c_m_tot = c_m - 1.0 * _f_con * dt;

      //   // op->update(make_ref(A_m_tot));
        
      //   // op->apply(_f_con, c_m);

      //   // CopyMatrixSolution(c_m);
  



