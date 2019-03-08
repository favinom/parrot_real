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



#include "utopia.hpp"
#include "utopia_Socket.hpp"
#include "utopia_FractureFlowUtils.hpp"
#include "utopia_TransferUtils.hpp"
#include "utopia_TransferAssembler.hpp"
#include "FractureNetworkUserObject.h"
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


#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"

typedef utopia::DSMatrixd SparseMatT;
typedef utopia::DVectord VecT;


registerMooseObject("parrot_realApp",FractureNetworkUserObject);

template <>
InputParameters
validParams<FractureNetworkUserObject>()
{
    //  MooseEnum proj_type("L2 pseudo-L2");
    
    InputParameters params = validParams<GeneralUserObject>();
    // params.addRequiredParam<AuxVariableName>("lagrange_variable",
    //                                          "The auxiliary variable to store the transferred values in.");
    params.addRequiredParam<VariableName>("matrix_variable",
                                          "The variable to transfer from.");
    params.addRequiredParam<VariableName>("fracture_variable",
                                          "The variable to transfer from.");
    params.addRequiredParam<UserObjectName>("operator_userobject","The userobject that stores our operators");
    params.addParam<std::string>("operator_type","opeartor_type" /*INTERPOLATION| L2_PROJECTION| PSEUDO_L2_PROJECTION | APPROX_L2_PROJECTION*/);
    params.addParam<bool>("solve_cg","transpose option");
    params.addParam<double>("porosity_m","posrosity matrix");
    params.addParam<double>("porosity_f","posrosity fracture");
    params.addParam<bool>("pressure","false","put true if you solve pressure system");
    params.addParam<bool>("transport","false","put true if you transport system");
    params.addParam<bool>("constraint_m","false","put true if matrix has Dirichlet BC");
    params.addParam<bool>("constraint_f","false","put true if fibres has Dirichlet BC");
    params.addRequiredParam<MultiAppName>("multi_app", "The MultiApp's name in your input file!");
    params.addParam<int>("id_slave","The boundary ID associated with the master side");
    //params.addParam<std::vector<BoundaryID>>("id_master","The boundary ID associated with the slave  side");
    // params.addParam< std::string >( "meshes",
    //                                "The name of the mesh files for the coarser levels"
    //                                " (must be an exodusII file). Listed from finest to coarsest.");
    return params;
    
}

FractureNetworkUserObject::FractureNetworkUserObject(const InputParameters & parameters):
GeneralUserObject(parameters),
// _name_mesh(parameters.get< std::string >("meshes")),
_slave_id(parameters.get<int>("id_slave")),
// _master_id(parameters.get<std::vector<BoundaryID>>("contact_master")),
_f_var_name(getParam<VariableName>("fracture_variable")),
_m_var_name(getParam<VariableName>("matrix_variable")),
_operator_storage(getUserObject<StoreTransferOperators>("operator_userobject")),
_operator_type(getParam<std::string>("operator_type")),
_solve_cg(getParam<bool>("solve_cg")),
_pressure(getParam<bool>("pressure")),
_transport(getParam<bool>("transport")),
_constraint_m(getParam<bool>("constraint_m")),
_constraint_f(getParam<bool>("constraint_f")),
_porosity_m(getParam<double>("porosity_m")),
_porosity_f(getParam<double>("porosity_f")),
_multiapp_name(getParam<MultiAppName>("multi_app"))

{
//    std::vector<std::string> tmp = utopia::split_string(parameters.get<std::string>("dc_boundaries"), ' ');
//    for(auto str_tmp=tmp.begin(); str_tmp != tmp.end(); str_tmp++)
//    {
//        _dc_boundary_id.push_back(atoi(str_tmp->c_str()));
//    }
//    // reading BC variables corresponding to BC ids
//    FractureNetworkApp::determine_dc_bnd_var_id(utopia::split_string(parameters.get<std::string>("dc_variables"), ' '),_n_master_variables);
}




void
FractureNetworkUserObject::initialize()
{
    _console << "Initial Setup of Fracture App " << std::endl;
    
    using namespace utopia;
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
    
    FEProblemBase & from_problem =_multi_app.problemBase();
    
    FEProblemBase & to_problem = _multi_app.appProblemBase(0);
    
    MeshBase *_m_mesh  = &to_problem.mesh().getMesh();
    
    MeshBase *_f_mesh  = &from_problem.mesh().getMesh();
    
    MooseVariable & _f_var = from_problem.getStandardVariable(0,_f_var_name);
    
    MooseVariable & _m_var = to_problem.getStandardVariable(0, _m_var_name);
    
    BoundaryMesh _b_mesh(_f_mesh->comm(), cast_int<unsigned char>(_f_mesh->mesh_dimension()-1));

    _f_mesh->get_boundary_info().sync(_b_mesh);

    EquationSystems _b_es (_b_mesh);

    _b_es.add_system<LinearImplicitSystem> ("boundary_sys");

    auto _l_var_num = _b_es.get_system("boundary_sys").add_variable("lambda", FIRST); 

    _b_es.init();

    libMesh::Variable _l_var = _b_es.get_system("boundary_sys").variable(_l_var_num);

    //bool is_shell = _b_mesh.mesh_dimension() < _b_mesh.spatial_dimension();

    auto V_m = utopia::LibMeshFunctionSpace(utopia::make_ref(_m_var.sys().system().get_equation_systems()),_m_var.sys().system().number(), _m_var.number());
    
    auto V_f = utopia::LibMeshFunctionSpace(utopia::make_ref(_f_var.sys().system().get_equation_systems()),_f_var.sys().system().number(), _f_var.number());

    auto V_l = utopia::LibMeshFunctionSpace(utopia::make_ref(_b_es),_b_es.get_system<LinearImplicitSystem>("boundary_sys").number(), _l_var.number());

    _P = std::make_shared<SparseMatT>();

    bundary_volume_permulation_matrix(_b_es, _l_var_num, _f_var, _P);
    
    assemble_projection(V_m, V_l, B_temp, D_temp);

    D_temp *=-1.;

    D = transpose(*_P) * D_temp * (*_P);

    //D+= -1.0 * local_identity(local_size(D).get(0),local_size(D).get(1));

    B = transpose(*_P) * B_temp;
    
    //B+= 1.0 * identity(B.size().get(0),B.size().get(1));

    D_t = transpose(D);
    
    B_t = transpose(B);

   
    utopia::disp(B.size());

    utopia::disp(D.size());
 
    // DVectord v = local_values(local_size(B_temp).get(1), 1.);

    // utopia::USparseMatrix D = transpose(*_P) * D_temp * (*_P);

    // utopia::disp(B_temp.size());

    // utopia::disp((*_P).size());

    // utopia::USparseMatrix B = transpose(*_P) * B_temp;

    // utopia::disp(B_MF.size());

    // DVectord u, sol;

    // u = B * v;

    // local_zeros(local_size(u));

    // auto lsolver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>();
    
    // lsolver->solve(D, u, sol);

    // disp(sol);

    // disp(sol.size());

    // CopyFractureSolution(sol);

    // utopia::disp(B.size().get(0));

    // utopia::disp(B.size().get(1));

    // utopia::disp(D.size().get(0));

    // utopia::disp(D.size().get(1));
}

void FractureNetworkUserObject::
bundary_volume_permulation_matrix(EquationSystems &_b_es, unsigned int _b_var_num, MooseVariable & _v_var, std::shared_ptr<SparseMatT> P_matrix)
{

  _console << "Permutation Matrix " << std::endl;

  // Get references to the Systems from the Variables
  libMesh::Variable _b_var_lib = _b_es.get_system("boundary_sys").variable(_b_var_num);
  
  libMesh::Variable _v_var_lib = _v_var.sys().system().variable(_v_var.number());


  // Get system number, variable numbers and number of components
  const unsigned int v_sys_number = _v_var.sys().system().number();
  const unsigned int v_var_number = _v_var_lib.number();

  const unsigned int b_sys_number = _b_es.get_system<LinearImplicitSystem>("boundary_sys").number();
  const unsigned int b_var_number = _b_var_lib.number();

  const unsigned int v_n_comp = _v_var_lib.n_components();
  const unsigned int b_n_comp = _b_var_lib.n_components();

  typedef std::map<numeric_index_type, numeric_index_type> DofMapping;

  DofMap & v_dofmap = _v_var.sys().system().get_dof_map();

  DofMap & b_dofmap = _b_es.get_system<LinearImplicitSystem>("boundary_sys").get_dof_map(); 

  (*P_matrix)=utopia::local_sparse(b_dofmap.n_local_dofs(),v_dofmap.n_local_dofs(),1.0);
  
  DofMapping dof_mapping;

  std::vector<dof_id_type> dof_boundary;

  {
        utopia::Write<SparseMatT> d(*(P_matrix));

        auto RowRange=utopia::row_range(*P_matrix);

        auto ColRange=utopia::col_range(*P_matrix);
       
        // loop through all boundary elements.
        for (const auto & b_elem : _b_es.get_mesh().active_local_element_ptr_range())
        {

            const Elem * v_elem = b_elem->interior_parent();


            // loop through all nodes in each boundary element.
            for (unsigned int node=0; node < b_elem->n_nodes(); node++)
            {

                // Node in boundary element.
                const Node * b_node = b_elem->node_ptr(node);

                for (unsigned int node_id=0; node_id < v_elem->n_nodes(); node_id++)
                {
                    // Nodes in interior_parent element.
                    const Node * v_node = v_elem->node_ptr(node_id);

                    const dof_id_type v_dof = v_node->dof_number(v_sys_number,
                                                                 v_var_number,
                                                                 v_n_comp - 1);
     
                    // See if we've already encountered this DOF in the loop
                    // over boundary elements.
                    DofMapping::iterator it = dof_mapping.find(v_dof);

                    // If we've already mapped this dof, we don't need to map
                    if (it == dof_mapping.end() &&
                        v_node->absolute_fuzzy_equals(*b_node, TOLERANCE))
                    {
                        // Global dof_index for node in BoundaryMesh
                        const dof_id_type b_dof = b_node->dof_number(b_sys_number,
                                                                     b_var_number,
                                                                     b_n_comp - 1);
                        dof_mapping[v_dof] = b_dof;

                        if (RowRange.inside(b_dof) && ColRange.inside(v_dof)){
                                (*P_matrix).set(b_dof,v_dof,1.0);
                        }
                    }
                }
            }
        }
    }
      
    // utopia::disp(B);
    //utopia::disp((*P_matrix).size());
    // utopia::disp(*P_matrix);

}


bool
FractureNetworkUserObject::solve(){
    
    if (_solve_cg){
        ok = solve_cg_dual();
    }
    else{
        ok = solve_monolithic();
    }
    //solve_cg_dual();
    
      return ok;
    
}


void
FractureNetworkUserObject::execute()
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

        solve_transport_monolithic();
    }
}


bool FractureNetworkUserObject::solve_monolithic()
{
    _console << "Solve_monolithic()"  << std::endl;
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
    
    FEProblemBase & _f_problem =_multi_app.problemBase();
    
    FEProblemBase & _m_problem = _multi_app.appProblemBase(0);
    
    MeshBase *_m_mesh = &_m_problem.mesh().getMesh();
    
    MeshBase *_f_mesh = &_f_problem.mesh().getMesh();
    
    MooseVariable & _m_var = _m_problem.getStandardVariable(0, _m_var_name);
    
   // MooseVariable & _l_var = _f_problem.getStandardVariable(0, _f_var_name);
    
    MooseVariable & _f_var   = _f_problem.getStandardVariable(0, _f_var_name);
    
    auto V_m = utopia::LibMeshFunctionSpace(utopia::make_ref(_m_var.sys().system().get_equation_systems()),_m_var.sys().system().number(), _m_var.number());
    
    auto V_f = utopia::LibMeshFunctionSpace(utopia::make_ref(_f_var.sys().system().get_equation_systems()),_f_var.sys().system().number(), _f_var.number());
    
//    to_problem.es().print_info();
    
    auto &_f_sys = _f_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_f = _f_problem.getNonlinearSystemBase();
    // Pointer to underlying PetscMatrix type
    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_f = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_f_sys.matrix);
    
    _f_problem.computeResidualSys(_f_sys, *_nl_f.currentSolution(), _nl_f.RHS());
    
    _f_problem.computeJacobianSys(_f_sys, *_nl_f.currentSolution(), *petsc_mat_f);
    
//    from_problem.es().print_info();
    
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

    _console << "Solve_monolithic():: START SOLVING"  << std::endl;


    utopia::disp(B.size());

    utopia::disp(D.size());

    utopia::disp(B_t.size());

    utopia::disp(D_t.size());

    utopia::disp(A_m.size());

    utopia::disp(A_f.size());

    utopia::USparseMatrix A = utopia::Blocks<utopia::USparseMatrix>(3, 3,
                                            {
                                                utopia::make_ref(A_m), nullptr, utopia::make_ref(B_t),
                                                nullptr, utopia::make_ref(A_f), utopia::make_ref(D_t),
                                                utopia::make_ref(B), utopia::make_ref(D), nullptr
                                            });

    utopia::write("mat.m", A);
    
    utopia::UVector z = utopia::local_zeros(local_size(B).get(0));

    utopia::UVector rhs = utopia::blocks(rhs_m, rhs_f, z);
    
    x_m  = utopia::local_zeros(local_size(rhs_m));

    x_f  = utopia::local_zeros(local_size(rhs_f));

    lagr = utopia::local_zeros(local_size(z));

    //lagr = utopia::local_zeros(local_size(D_temp).get(0));
    
    utopia::UVector sol = blocks(x_m, x_f, lagr);
    
    utopia::Factorization<utopia::USparseMatrix, utopia::UVector> op;
   
    op.update(make_ref(A));
    
    bool ok = op.apply(rhs, sol);

    _console << "Solve_monolithic():: FINISH SOLVING"  << std::endl;
    
    utopia::undo_blocks(sol, x_m, x_f, lagr);
    
    return ok;
}
bool FractureNetworkUserObject::solve_cg_dual(){
    
    _console << "Solve_cg_dual()  "  << std::endl;
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
    
    FEProblemBase & from_problem =_multi_app.problemBase();
    
    FEProblemBase & to_problem = _multi_app.appProblemBase(0);
    
    MeshBase *from_mesh = &from_problem.mesh().getMesh();
    
    MeshBase * to_mesh = &to_problem.mesh().getMesh();
    
    MooseVariable & _m_var = from_problem.getStandardVariable(0, _m_var_name);
    
    MooseVariable & _l_var = to_problem.getStandardVariable(0, _f_var_name);
    
    MooseVariable & _f_var   = to_problem.getStandardVariable(0, _f_var_name);
    
    auto V_m = utopia::LibMeshFunctionSpace(utopia::make_ref(_m_var.sys().system().get_equation_systems()),_m_var.sys().system().number(), _m_var.number());
    
    auto V_f = utopia::LibMeshFunctionSpace(utopia::make_ref(_f_var.sys().system().get_equation_systems()),_f_var.sys().system().number(), _f_var.number());
    
    utopia::SPBlockConjugateGradient<utopia::USparseMatrix, utopia::UVector> solver;
    
    solver.verbose(true);
    
    solver.max_it(2000);
    
    solver.atol(1e-14);
    
    solver.use_simple_preconditioner();
    
    bool use_mg;
    
    use_mg = false;
    int mg_sweeps = 1;
    int mg_levels = 5;
    
    if(use_mg) {
        auto mg = make_mg_solver(V_m, mg_levels);
        solver.set_master_solver(mg);
        
        solver.set_master_sweeps(mg_sweeps);
        solver.set_master_max_it(mg->max_it());
    }
    
    //to_problem.es().print_info();
    
    auto &_to_sys = to_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_to = to_problem.getNonlinearSystemBase();
    // Pointer to underlying PetscMatrix type
    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_to = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_to_sys.matrix);
    
    to_problem.computeResidualSys(_to_sys, *_nl_to.currentSolution(), _nl_to.RHS());
    
    to_problem.computeJacobianSys(_to_sys, *_nl_to.currentSolution(), *petsc_mat_to);
    
    //from_problem.es().print_info();
    
    auto &_from_sys = from_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");
    
    NonlinearSystemBase & _nl_from = from_problem.getNonlinearSystemBase();
    
    // Pointer to underlying PetscMatrix type
    libMesh::PetscMatrix<libMesh::Number> *petsc_mat_from = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_from_sys.matrix);
    
    from_problem.computeResidualSys(_from_sys, *_nl_from.currentSolution(), _nl_from.RHS());
    
    from_problem.computeJacobianSys(_from_sys, *_nl_from.currentSolution(), *petsc_mat_from);
    
    utopia::UVector rhs_f;
    utopia::USparseMatrix A_f;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_to.RHS()), rhs_f);
    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_to).mat(), A_f);
    //disp(mat_to);
    
    utopia::UVector rhs_m;
    utopia::USparseMatrix A_m;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_from.RHS()), rhs_m);
    utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_from).mat(), A_m);
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
  bool ok = solver.apply(rhs_m, rhs_f, x_m, x_f, lagr);
    
  return ok;
    
    
}

void
FractureNetworkUserObject::CopyMatrixSolution(utopia::UVector _sol_m)
{
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);

    FEProblemBase & _problem_m = _multi_app.appProblemBase(0);
    
    MooseVariable & _var_m = _problem_m.getStandardVariable(0, _m_var_name);
    
    System & _sys_m = _var_m.sys().system();
    
    //FEProblemBase & to_problem = _multi_app.appProblemBase(0);
    
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

    if(_pressure)
       ExodusII_IO (*_mesh_m).write_equation_systems("matrix_p.e", _var_m.sys().system().get_equation_systems());
    else
       ExodusII_IO (*_mesh_m).write_equation_systems("matrix_c.e", _var_m.sys().system().get_equation_systems());
    
    
}

void
FractureNetworkUserObject::CopyFractureSolution(utopia::UVector _sol_f)
{
    
    
    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);

    FEProblemBase & _problem_f =_multi_app.problemBase();
    
//    std::cout<<" _f_var_name" << _f_var_name<<std::endl;
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
FractureNetworkUserObject::constraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
{
//    auto &V = space->space().last_subspace();

    using namespace utopia;
    
    {
        Write<utopia::UVector> w_v(vec);

        Read<utopia::UVector> r_v(boundary);

        if(has_constaints) {

            Range r = range(vec);

            for(SizeType i = r.begin(); i < r.end(); ++i) {

                //std::cout<<"value "<< boundary.get(i) <<std::endl;

                if(boundary.get(i)!=0) {

                    //std::cout<<"i "<< i <<"valpos->second "<<boundary.get(i)<<std::endl;

                    // if(valpos != rhs_values.end()) {
                    vec.set(i, boundary.get(i));
                    // }
                }
            }
        }
    }

    synchronize(vec);
    
}

void
FractureNetworkUserObject::constraint_concentration_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints)
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
FractureNetworkUserObject::assemble_mass_matrix(double porosity, FEProblemBase & _problem, utopia::USparseMatrix &mass_matrix){
    
    _console << "Assemble_Mass_matrix() begin "  << std::endl;
    
    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _problem.es().get_mesh();
    
    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();
    
    _problem.es().add_system<LinearImplicitSystem>("aux").add_variable("var",FIRST);
    
    _problem.es().reinit();
    
    // Get a reference to our system.
    LinearImplicitSystem & _system = _problem.es().get_system<LinearImplicitSystem>("aux");
    
    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _system.get_dof_map().variable_type(0);
    
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    
    SparseMatrix<Number> & matrix_M = *_system.matrix;
    
    //_console << "is matrix closed: " << matrix_A.closed() << std::endl;
    
    // The element mass matrix.
    DenseMatrix<Number> Me;
    
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
    matrix_M.zero();
    
    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;
        Elem * ele = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        fe->reinit (elem);
        
        Me.resize (dof_indices.size(), dof_indices.size());
        
        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
            for (unsigned int i=0; i<phi.size(); i++)
                for (unsigned int j=0; j<phi.size(); j++)
                        Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
        
        //dof_map.constrain_element_matrix(Me, dof_indices, false);
        //        Me.print();
        
        matrix_M.add_matrix (Me, dof_indices);
        
    }
    
    matrix_M.close();
    
    utopia::USparseMatrix mass;

    utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M), mass);

    mass_matrix = porosity * mass;
    
    //out = mass_matrix * in;
    
    _console << "Assemble_Mass_matrix() end "  << std::endl;

    
}


void
FractureNetworkUserObject::solve_transport_monolithic(){
    
    using namespace utopia;
    
    _console << "solve_transport_monolithic()  "  << std::endl;

    MultiApp &  _multi_app = * _fe_problem.getMultiApp(_multiapp_name);
        
    FEProblemBase & from_problem =_multi_app.problemBase();
    
    FEProblemBase & to_problem = _multi_app.appProblemBase(0);

    auto &_to_sys = to_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    NonlinearSystemBase & _nl_to = to_problem.getNonlinearSystemBase();
    
    auto &_from_sys = from_problem.es().get_system<TransientNonlinearImplicitSystem>("nl0");

    NonlinearSystemBase & _nl_from = from_problem.getNonlinearSystemBase();



    
 
    if  (from_problem.timeStep()==1)
     {
              
        // Pointer to underlying PetscMatrix type
        libMesh::PetscMatrix<libMesh::Number> *petsc_mat_to = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_to_sys.matrix);
        
        to_problem.computeJacobianSys(_to_sys, *_nl_to.currentSolution(), *petsc_mat_to);
        
        // Pointer to underlying PetscMatrix type
        libMesh::PetscMatrix<libMesh::Number> *petsc_mat_from = dynamic_cast<libMesh::PetscMatrix<libMesh::Number>* >(_from_sys.matrix);

        
        from_problem.computeJacobianSys(_from_sys, *_nl_from.currentSolution(), *petsc_mat_from);
        
        
        utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_to).mat(), A_f_t);

        utopia::convert(const_cast<libMesh::PetscMatrix<libMesh::Number> &>(*petsc_mat_from).mat(), A_m_t);

        from_problem.computeResidualSys(_from_sys, *_nl_from.currentSolution(), _nl_from.RHS());
  
        to_problem.computeResidualSys(_to_sys, *_nl_to.currentSolution(), _nl_to.RHS());

        utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_to.RHS()), rhs_f_t);

        utopia::convert(const_cast<NumericVector<libMesh::Number> &>(_nl_from.RHS()), rhs_m_t);


        rhs_m_c = rhs_m_t; 

        rhs_f_c = rhs_f_t; 


        Real dt = static_cast<Transient*>(from_problem.getMooseApp().getExecutioner())->getDT();

        double inv_dt = 1.0/dt;

        USparseMatrix A_m_tot = A_m_t + mass_m * inv_dt;

        constraint_concentration_mat(rhs_m_c, A_m_tot, _constraint_m);

        USparseMatrix A_f_tot = A_f_t + mass_f * inv_dt;

        constraint_concentration_mat(rhs_f_c, A_f_tot, _constraint_f);


        USparseMatrix A_t = Blocks<USparseMatrix>(3, 3,
                                                {
                                                    make_ref(A_m_tot),              nullptr,                  make_ref(B_t),
                                                    nullptr,                   make_ref(A_f_tot),             make_ref(D_t),
                                                    make_ref(B),               make_ref(D),               nullptr
                                                });
        


        _A_t_store   = std::make_shared<USparseMatrix>(A_t);

        const_cast<StoreTransferOperators&>(_operator_storage).setTransferOperator() = _A_t_store;

        auto op = std::make_shared<utopia::Factorization<utopia::USparseMatrix, utopia::UVector> >();

        const_cast<StoreTransferOperators&>(_operator_storage).getVoidPointer() = op;

        dlagr_t = utopia::local_zeros(local_size(B).get(0));   

        lagr_t = utopia::local_zeros(local_size(dlagr_t));

        c_m  = utopia::local_zeros(local_size(rhs_m_t));

        c_f  = utopia::local_zeros(local_size(rhs_f_t));


        CopyMatrixSolution(c_m);
        
        CopyFractureSolution(c_f);
    
    }


    if (from_problem.timeStep()>1)
    {



        auto sol_m_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_from_sys.old_local_solution.get())->vec();

        auto sol_f_old = static_cast<libMesh::PetscVector<libMesh::Number> *>(_to_sys.old_local_solution.get())->vec();

        // UVector c_m_old, c_f_old;

        utopia::convert(sol_m_old, c_m_old);

        utopia::convert(sol_f_old, c_f_old);

        // disp("c_m_old");
        
        // disp(c_m_old);

        // disp("rhs_m_t");
        
        // disp(rhs_m_t);

        Real dt = static_cast<Transient*>(from_problem.getMooseApp().getExecutioner())->getDT();

        double inv_dt = 1.0/dt;

        UVector  mass_c_m_old,  mass_c_f_old;

        mass_c_m_old = mass_m * c_m_old;

        mass_c_f_old = mass_f * c_f_old;

        mass_c_f_old_dot = mass_c_f_old * inv_dt;

        mass_c_m_old_dot = mass_c_m_old * inv_dt;

        rhs_m_t = -1.0 * mass_c_m_old_dot;

        rhs_f_t = -1.0 * mass_c_f_old_dot;

   

        constraint_concentration_vec(rhs_m_c,  rhs_m_t, _constraint_m);
    
        constraint_concentration_vec(rhs_f_c,  rhs_f_t, _constraint_f);

        
        auto op = std::static_pointer_cast< Factorization<USparseMatrix, UVector> >(const_cast<StoreTransferOperators&>(_operator_storage).getVoidPointer());

        utopia::UVector rhs_t = utopia::blocks(rhs_m_t, rhs_f_t, dlagr_t);
        
        utopia::UVector sol_t = blocks(c_m, c_f, lagr_t);

         _A_t_store = const_cast<StoreTransferOperators&>(_operator_storage).getTransferOperator();
        
        op->update(_A_t_store);
        
        op->apply(rhs_t, sol_t);

        utopia::undo_blocks(sol_t, c_m, c_f, lagr_t);

        c_m *=1;

        c_f*=1;

        CopyMatrixSolution(c_m);
        
        CopyFractureSolution(c_f);
    }

}

        //UVector sol;

        // Real dt = static_cast<Transient*>(from_problem.getMooseApp().getExecutioner())->getDT();

        // double inv_dt = 1.0/dt;

        //A_m_t*=dt;

        // USparseMatrix A_m_tot = A_m_t + mass_m * inv_dt;

        // constrain_concentration_mat(rhs_m_c, A_m_tot, true);

        // USparseMatrix A_f_tot = A_f_t + mass_f * inv_dt;

        // constrain_concentration_mat(rhs_f_c, A_f_tot, false);


        // USparseMatrix A_tot = Blocks<USparseMatrix>(3, 3,
        //                                         {
        //                                             make_ref(A_m_tot),              nullptr,                  make_ref(B_t),
        //                                             nullptr,                   make_ref(A_f_tot),             make_ref(D_t),
        //                                             make_ref(B),               make_ref(D),               nullptr
        //                                         });
        


        //op->update(make_ref(A_m_t));
        
        //op->apply(rhs_m_t, sol);

        //op->update(make_ref(A_tot));

        // write("A_m_tot.m", A_m_tot);

        // write("rhs_m.m", rhs_m_t);
        
        

     


        // c_m_old = c_m;

        // c_f_old = c_f;

        // disp("sol");

        // disp(sol);


// void
// L2ProjectionLibMeshTransferBidirectional::Permutation(FEProblemBase *problem, std::shared_ptr<SparseMatT> B_Mm, std::shared_ptr<SparseMatT> B_Ss, std::string _complete_var_name, std::string problem_name, unsigned int number){


//     std::cout<<"Permutation"<<std::endl;

//     MeshBase & mesh = problem->mesh();

//     MeshBase::const_element_iterator it = mesh.local_elements_begin();

//     const MeshBase::const_element_iterator end_it = mesh.local_elements_end();

//     std::vector<dof_id_type> temp_main, temp_aux;

//     MooseVariable & main_var = problem->getStandardVariable(0, main_var_name);

//     DofMap & complete_dof = main_var.sys().system().get_dof_map(); //solid_problem

//     DofMap & _main_dof = problem->es().get_system(problem_name).get_dof_map(); //auxi_var.sys().system().get_dof_map();  //fluid_problem

//     std::cout<<"Permutation Loop"<<std::endl;

//     (*P_matrix)=utopia::local_sparse(auxi_dof.n_local_dofs(),main_dof.n_local_dofs(),1.0);

//     std::cout<<"Permutation After allocation"<<std::endl;
//    {
//      utopia::Write<SparseMatT> d(*(P_matrix));
//      auto RowRange=utopia::row_range(*P_matrix);
//      auto ColRange=utopia::col_range(*P_matrix);

//      for ( ; it != end_it; ++it)
//      {
//          const Elem * elem = *it;

//          main_dof.dof_indices(elem, temp_main, main_var.number());

//          auxi_dof.dof_indices(elem, temp_aux, number);

//          for(int i=0; i<temp_aux.size(); i++){

//             const long dof_I = temp_aux[i];

//             const long dof_J = temp_main[i];

//             if (RowRange.inside(dof_I) && ColRange.inside(dof_J)){

//                 for(utopia::SizeType d = 0; d < _n_var; ++d){
//                    (*P_matrix).set(dof_I+d,dof_J+d,1.0);
//                 }
//             }
//         }
//      }
//    }
// }



    



