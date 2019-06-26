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
    params.addRequiredParam<std::vector<int>>("block_id",
                                                     "The name of the nodeset to create");
    params.addRequiredParam<std::vector<Real>>("value_p",
                                                     "The name of the nodeset to create");
    return params;
    
}

FractureAppConforming::FractureAppConforming(const InputParameters & parameters):
GeneralUserObject(parameters),
_m_var_name(getParam<VariableName>("matrix_variable")),
_vector_p(getParam<std::vector<int>>("block_id")),
_vector_value(getParam<std::vector<Real>>("value_p")),
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

         // Get a reference to our system.
   if (_fe_problem.timeStep()==0) {
    _fe_problem.es().add_system<LinearImplicitSystem>("aux1").add_variable("var",FIRST);
   
    _fe_problem.es().add_system<LinearImplicitSystem>("aux2").add_variable("var",FIRST);
    
    _fe_problem.es().reinit();
  }
    

    
    
    
}


void
FractureAppConforming::execute()
{

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
    LinearImplicitSystem & _system = _fe_problem.es().get_system<LinearImplicitSystem>("aux1");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = _system.get_dof_map().variable_type(0);

    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

    SparseMatrix<Number> & matrix_M = *_system.matrix;

 
    DenseMatrix<Number> Me;

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


        //dof_map.constrain_element_matrix(Me,dof_indices,false);

        matrix_M.add_matrix (Me, dof_indices);



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
FractureAppConforming::assemble_poro_mass_matrix(FEProblemBase & _problem, utopia::USparseMatrix &mass_matrixp, utopia::USparseMatrix &lumped_mass_matrixp){
    
    _console << "Assemble_Poro_Mass_matrix() begin "  << std::endl;

    Real dt = static_cast<Transient*>(_fe_problem.getMooseApp().getExecutioner())->getDT();

    // Get a constant reference to the mesh object.
    const MeshBase & mesh = _problem.es().get_mesh();

    // The dimension that we are running.
    const unsigned int dim = mesh.mesh_dimension();

    std::cout<<"dim"<<dim<<std::endl;
    
    LinearImplicitSystem & _system = _fe_problem.es().get_system<LinearImplicitSystem>("aux2");

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    
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
        
        fe->reinit (elem);

        dof_map.dof_indices(elem, dof_indices);
        
        Me_p.resize (dof_indices.size(), dof_indices.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

            for (unsigned int i=0; i<phi.size(); i++){

                for (unsigned int j=0; j<phi.size(); j++){

                        Me_p(i,j) +=   ComputeMaterialProprties(elem) * JxW[qp] * phi[i][qp] * phi[j][qp];
          
            }
        }
    }

        
        dof_map.constrain_element_matrix(Me_p,dof_indices,false);

        matrix_M_p.add_matrix (Me_p, dof_indices);

        
    }


    
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
FractureAppConforming::stabilize_A_matrix(FEProblemBase & _problem, utopia::USparseMatrix &A_0, utopia::USparseMatrix &S_matrix)
{
    
    
    _console << "Stabilize A matrix:: begin  "  << std::endl;

    utopia::USparseMatrix A_0_t;
    
    A_0_t = utopia::transpose(A_0);

    S_matrix = A_0_t;

    S_matrix*=0;
    
    
    
    
    {
        utopia::Read<utopia::USparseMatrix>  r_s(A_0), r_s_t(A_0_t);
        utopia::Write<utopia::USparseMatrix> w_s(S_matrix);
        utopia::each_read(A_0, [&](const utopia::SizeType i, const utopia::SizeType j, double value){
            if(i!=j)
            {
                double value_1 = 1.0 * value;
                
                //utopia::disp(value_1);
                
                double value_2 = 1.0 * A_0_t.get(i,j);

                Real max=std::max(value_1,value_2);

                if (max>0.0){

                    max*=-1.0;

                    S_matrix.set(i,j,max);
                }
            }
                
           
        });
    }
    
    
    
    
    utopia::UVector diag_elem = -1.0 * sum(S_matrix,1);
    
    utopia::USparseMatrix S_diag=diag(diag_elem);
    
    // S_matrix += transpose(S_matrix);
    
    // S_matrix *=0.5;
    
    S_matrix+=S_diag;
    
    _console << "Stabilize A matrix:: end  "  << std::endl;
        
    
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



Real
FractureAppConforming::ComputeMaterialProprties(const Elem *elem){
    
   // _console << "_vector_p.size()"  << _vector_p.size() <<std::endl;

Real poro=0.0;

    for(int ll=0; ll<_vector_p.size(); ll++){
        if (elem->subdomain_id()==_vector_p[ll]) {

            poro = _vector_value[ll]; 
        }
    }
  
    return poro;   
}


void 
FractureAppConforming::stabilize_coeffiecient(utopia::UVector &rhs_m_c, FEProblemBase & _problem, utopia::UVector &_u, utopia::UVector &_u_dot, utopia::USparseMatrix &_D, utopia::USparseMatrix &_M, utopia::UVector &_rhs)
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
     
    boundary_constraint_vec(rhs_m_c,  R_minus, true);

    boundary_constraint_vec(rhs_m_c,  R_plus, true);



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


    _console << "stabilize_coeffiecient::end"  << std::endl;
    
   
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
