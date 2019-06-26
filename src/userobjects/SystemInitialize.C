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


#include <algorithm>    // std::max


#include "SystemInitialize.h"
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


registerMooseObject("parrot_realApp",SystemInitialize);

template <>
InputParameters
validParams<SystemInitialize>()
{
    
    InputParameters params = validParams<GeneralUserObject>();
    params.addRequiredParam<VariableName>("matrix_variable",
                                          "The variable to transfer from.");

    return params;
    
}

SystemInitialize::SystemInitialize(const InputParameters & parameters):
GeneralUserObject(parameters),
_m_var_name(getParam<VariableName>("matrix_variable"))


{
    
}




void
SystemInitialize::initialize()
{

  


         //Get a reference to our system.
    // _fe_problem.es().add_system<TransientImplicitSystem>("stiffness").add_variable("var_new",FIRST);

    // _fe_problem.es().reinit();

    // _fe_problem.es().print_info();

  

    

    
    
}


void
SystemInitialize::execute()
{
//    solve_stabilize_monolithic();
}

void
SystemInitialize::stabilize_A_matrix(FEProblemBase & _problem, utopia::USparseMatrix &A_0, utopia::USparseMatrix &S_matrix)
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

    utopia::write("New_S.m", S_matrix);
        
    
}


void
SystemInitialize::unconstraint_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints)
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
        
        set_zero_rows(mat, rows, 0.);
    }
}

void
SystemInitialize::constraint_mat(utopia::UVector &boundary, utopia::USparseMatrix &mat, bool has_constaints)
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
SystemInitialize::boundary_constraint_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
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
SystemInitialize::CopyMatrixSolution(utopia::UVector _sol_m)
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
    //ExodusII_IO (*_mesh_m).write_equation_systems("matrix_c.e", _var_m.sys().system().get_equation_systems());
    

}


void
SystemInitialize::unconstraint_concentration_vec(utopia::UVector &boundary, utopia::UVector &vec, bool has_constaints)
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
SystemInitialize::assemble_mass_matrix(FEProblemBase & _problem, utopia::USparseMatrix &mass_matrix, utopia::USparseMatrix &lumped_mass_matrix){
    
    _console << "Assemble_Mass_matrix() begin "  << std::endl;

    Real dt = static_cast<Transient*>(_fe_problem.getMooseApp().getExecutioner())->getDT();

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


        // dof_map.constrain_element_matrix(Me,dof_indices,false);

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






