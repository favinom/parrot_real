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
#include "utopia_MeshTransferOperator.hpp"




typedef utopia::DSMatrixd SparseMatT;
typedef utopia::DVectord VecT;


registerMooseObject("parrot_realApp",Test_Assembly);

template <>
InputParameters
validParams<Test_Assembly>()
{

    InputParameters params = validParams<GeneralUserObject>();
    return params;
    
}

Test_Assembly::Test_Assembly(const InputParameters & parameters):
GeneralUserObject(parameters)
{

}

void
Test_Assembly::initialize()
{
    _console << "Initial Setup of Fracture App " << std::endl;

    assemble_mass_matrix(_fe_problem,  mass_m, mass_lumped);
}



void
Test_Assembly::execute()
{

}



void
Test_Assembly::assemble_mass_matrix(FEProblemBase & _problem, utopia::USparseMatrix &mass_matrix, utopia::USparseMatrix &lumped_mass_matrix){
    
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

    DenseVector<Number> Fe;


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

        //Elem * ele = *el;
        
        fe->reinit (elem);

        dof_map.dof_indices(elem, dof_indices);
        
        Me.resize(dof_indices.size(), dof_indices.size());

        Fe.resize(dof_indices.size());

        // Me_p.resize (dof_indices.size(), dof_indices.size());
        
        for (unsigned int qp=0; qp<qrule.n_points(); qp++){

            for (unsigned int i=0; i<phi.size(); i++){

                for (unsigned int j=0; j<phi.size(); j++){

                    if(i==j) Me(i,j) = 1.0;
            }
        }
    }


    for (unsigned int i = 0; i < phi.size(); i++){
        // RHS
        Fe(i) = 1.0;
    }


    // if(elem->id()==5){


         std::cout<<"Me"<<Me<<std::endl;

    //     //std::cout<<"Fe"<<Me<<std::endl;



    //     //std::cout<<"elem_id "<<elem->id()<<""<<*elem<<std::endl;


    //     for(int jj=0; jj<dof_indices.size(); jj++){
    //             std::cout<<"dof_indices "<<dof_indices[jj]<<std::endl;
    //             //std::cout<<"node "<<elem->node_ref(jj)<<std::endl;
    //     }
          
            
    //         //dof_map.constrain_element_matrix(Me,dof_indices,true);

    //     dof_map.constrain_element_vector(Fe,dof_indices,true);

        
    //     for(int jj=0; jj<dof_indices.size(); jj++){
    //             std::cout<<"dof_indices_after "<<dof_indices[jj]<<std::endl;
    //             //std::cout<<"node "<<elem->node_ref(jj)<<std::endl;
    //     }


       



    //     // exit(1);

    // }
        dof_map.constrain_element_matrix(Me,dof_indices,true);

        std::cout<<"Me_after"<<Me<<std::endl;

        matrix_M.add_matrix (Me, dof_indices);


        _system.rhs->add_vector(Fe, dof_indices);
        //matrix_M_p.add_matrix (Me_p, dof_indices);

        
    }


    
    matrix_M.close();
    _system.rhs->close();


    //matrix_M_p.close();

    utopia::SizeType nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()); //dof_map.get_n_nz();
        
    utopia::USparseMatrix mass_temp = utopia::local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

    utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M), mass_temp);

    utopia::UVector rhs_f;
    
    utopia::convert(const_cast<NumericVector<libMesh::Number> &>(*_system.rhs), rhs_f);

    disp(mass_temp);

    disp(rhs_f);

    // utopia::convert(const_cast<SparseMatrix<Number> &>(matrix_M_p), mass_mp);

    mass_matrix =  mass_temp;

    lumped_mass_matrix = diag(sum(mass_matrix,1));


    _console << "Assemble_Mass_matrix() end "  << std::endl;

    
}











/*    UVector P_minus = sum(transpose(P_m_minus),1);

    UVector P_max = sum(transpose(P_m_max),1);*/



    






  



