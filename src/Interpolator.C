#include "Interpolator.h"

// MOOSE includes
#include "FEProblem.h"
// MOOSEMESH not needed, strange

// libmesh includes
#include "libmesh/equation_systems.h"

//mesh_base not needed, strange
#include "libmesh/mesh_tools.h"
#include "libmesh/dof_map.h"
#include "libmesh/face_quad4.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/quadrature_gauss.h"

#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

//#include "libmesh/libmesh.h"
#include "libmesh/fe.h"

registerMooseObject("parrot_realApp", Interpolator);

// typedef utopia::DSMatrixd SparseMatT;
// typedef utopia::DVectord VecT;

template <>
InputParameters validParams<Interpolator>()
{
    // Get the input parameters from the parent class
    InputParameters params = validParams<GeneralUserObject>();
    // params.addParam<std::string>("dc_boundaries", "-1", "Dirichlet Boundary ID");
    // params.addParam<std::string>("dc_variables" , "-1", "Variable to which given BC_id applies");
    // params.addRequiredParam<std::vector< VariableName >>("primal_variables" , "Primal Variables names");
    // params.addParam<std::vector< VariableName >>("dual_variables" , "Dual Variables names");
    // params.addParam< std::string >("primal_order","Primal Variables Order");
    // params.addParam< std::string >("dual_order", "Dual Variables Order");
    return params;
}


Interpolator::Interpolator(const InputParameters & parameters) : GeneralUserObject(parameters),
_fe_problem(parameters.get<FEProblem *>("_fe_problem")), // il problema
_mooseMesh( _fe_problem[0].mesh() ),
_equationSystem(_fe_problem[0].es()),
_mesh( static_cast<Mesh *>(&_equationSystem.get_mesh()) ),
_lis(_equationSystem.get_system<LinearImplicitSystem> (0) ),
_dof_map(_lis.get_dof_map())
// _nl(_fe_problem->getNonlinearSystem()),
// _nVars(_nl.nVariables()),
// _primal_variable_order(parameters.get< std::string>("primal_order")),
// _dual_variable_order(parameters.get< std::string >("dual_order")),
// _primal_variable_name(getParam<std::vector< VariableName> > ("primal_variables")),
// _dual_variable_name(getParam<std::vector< VariableName> > ("dual_variables")),
{
    //    MeshBase const * meshBase = &_mooseMesh -> getMesh();
    //    int const levels0=libMesh::MeshTools::n_levels(meshBase[0]);
    //    std::cout<<"meshBase "<<levels0<<std::endl;
    
    _numOfLevels=libMesh::MeshTools::n_levels(_mesh[0]);
    std::cout<<"The read mesh has "<<_numOfLevels<<" levels!"<<std::endl;
    
    _meshes = new Mesh *[_numOfLevels];
    
    _meshDimension=_mesh[0].mesh_dimension ();
    _spatialDimension=_mesh[0].spatial_dimension ();
    
//    _local2global_nodes.resize(_numOfLevels);
//    _global2local_nodes.resize(_numOfLevels);
    _fine2coarse.resize(_numOfLevels-1);
    _coarse2fine.resize(_numOfLevels-1);
    
    // if(_primal_variable_order=="SECOND")
    //   FEType_primal=SECOND;
    // else
    //   FEType_primal=FIRST;
    
    // if(_dual_variable_order=="SECOND")
    //   FEType_dual=SECOND;
    // else
    //   FEType_dual=FIRST;
    
    // if(_primal_variable_name.size()< 1 && _dual_variable_name.size() < 1 )
    // mooseError("Please set _primal_variable_name and/or _primal_variable_name to have at least one element");
    //
    //
    // _mesh_dimension = static_cast<int>(_fe_problem->es().get_mesh().mesh_dimension());
    
    
    // reading BC id markers
    // std::vector<std::string> tmp = Interpolator::split(parameters.get<std::string>("dc_boundaries"), ' ');
    // for(auto str_tmp=tmp.begin(); str_tmp != tmp.end(); str_tmp++)
    // {
    // _dc_boundary_id.push_back(atoi(str_tmp->c_str()));
    // }
    
    // reading BC variables corresponding to BC ids
    // Interpolator::determine_dc_bnd_var_id(Interpolator::split(parameters.get<std::string>("dc_variables"), ' '));
}

void
Interpolator::execute()
{
    _meshes[_numOfLevels-1]=new Mesh(_mesh[0]);
    {
    std::map<dof_id_type, std::vector<dof_id_type> > hanging_nodes;
    libMesh::MeshTools::find_hanging_nodes_and_parents(_meshes[_numOfLevels-1][0],hanging_nodes);
    
    std::cout<<"Number of hanging nodes on level "<<_numOfLevels-1<<" is "<<hanging_nodes.size()<<std::endl;
    }
    
    
    createMeshOnLevel(2);
    createMeshOnLevel(1);
    createMeshOnLevel(0);


    
    std::string filename0("ciao0.e");
    _meshes[0][0].write(filename0);
    
    std::string filename1("ciao1.e");
    _meshes[1][0].write(filename1);
    
    std::string filename2("ciao2.e");
    _meshes[2][0].write(filename2);
    
    std::string filename3("ciao3.e");
    _meshes[3][0].write(filename3);
    
    //verifyMaps();
    
    std::string systemName("Poisson");
    _equationSystems = new EquationSystems * [_numOfLevels];
    for (int level=0; level<_numOfLevels; ++level)
    {
        _equationSystems[level]=new EquationSystems(_meshes[level][0]);
        //std::string systemName("system");
        //systemName=systemName+std::to_string(level);
        _equationSystems[level][0].add_system<LinearImplicitSystem> (systemName.c_str());
        _equationSystems[level][0].get_system(systemName.c_str()).add_variable("u", FIRST);
        _equationSystems[level][0].init();
    }
    
    LinearImplicitSystem& system_c = _equationSystems[0][0].get_system<LinearImplicitSystem> (systemName.c_str());
    LinearImplicitSystem& system_f = _equationSystems[1][0].get_system<LinearImplicitSystem> (systemName.c_str());
    
    NumericVector<Number> & solution_c =  *(system_c.solution.get());
    
    std::cout<<solution_c.size()<<std::endl;
    for (int i=0; i<solution_c.size(); ++i)
    {
        solution_c.set(i,1.0*i);
    }
    
    solution_c.close();
    
    ExodusII_IO (_meshes[0][0]).write_equation_systems ("out.e", _equationSystems[0][0]);
    
    _equationSystems[0][0].print_info();
    exit(1);
    
    const DofMap& dof_map_c = system_c.get_dof_map();
    const DofMap& dof_map_f = system_f.get_dof_map();
    
    FEType fe_type_c = dof_map_c.variable_type(0);
    FEType fe_type_f = dof_map_f.variable_type(0);
    
    std::shared_ptr<FEBase> fe_c (FEBase::build(2, fe_type_c));
    std::shared_ptr<FEBase> fe_f (FEBase::build(2, fe_type_f));
    
    QGauss qrule_c (2, FIFTH);
    QGauss qrule_f (2, FIFTH);
    
    fe_c->attach_quadrature_rule (&qrule_c);
    fe_f->attach_quadrature_rule (&qrule_f);
    
    const std::vector<Point>& q_point_c = fe_c->get_xyz();
    const std::vector<Point>& q_point_f = fe_f->get_xyz();
    
    const std::vector<std::vector<Real> >& phi_c = fe_c->get_phi();
    const std::vector<std::vector<Real> >& phi_f = fe_f->get_phi();

    FE<2, XYZ> myfe_c(fe_type_c);
    FE<2, XYZ> myfe_f(fe_type_f);
    
    MeshBase::const_element_iterator       el     = _meshes[0][0].active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = _meshes[0][0].active_local_elements_end();
    
    std::map<dof_id_type, Point > ciao;
    
    int count(0);
    
    PetscMatrix<Number> proj_matrix(_mesh[0].comm());
    proj_matrix.init(dof_map_f.n_dofs(),dof_map_c.n_dofs(),
                     dof_map_f.n_local_dofs (),dof_map_c.n_local_dofs ());
    
    std::cout<<dof_map_c.n_dofs()<<std::endl;
    std::cout<<dof_map_f.n_dofs()<<std::endl;
    std::cout<<dof_map_c.n_local_dofs ()<<std::endl;
    std::cout<<dof_map_f.n_local_dofs ()<<std::endl;
    
    for ( ; el != end_el ; ++el)
    {
        //std::cout<<"count="<<count++<<std::endl;
        
        Elem * elem_c = *el;
        Elem * elem_f = _coarse2fine.at(0).find(elem_c)->second;
        
        //std::cout<<"elem_c id = "<<elem_c[0].id()<<std::endl;
        //std::cout<<"elem_f id = "<<elem_f[0].id()<<std::endl;

        //std::cout<<"elem_c hc = "<<elem_c[0].has_children()<<std::endl;
        //std::cout<<"elem_f hc = "<<elem_f[0].has_children()<<std::endl;

        //std::cout<<"elem_c nc = "<<elem_c[0].n_children()<<std::endl;
        //std::cout<<"elem_f nc = "<<elem_f[0].n_children()<<std::endl;

        //std::cout<<"elem_c ia = "<<elem_c[0].active()<<std::endl;
        //std::cout<<"elem_f ia = "<<elem_f[0].active()<<std::endl;

        //std::cout<<elem_c[0]<<std::endl;
        //std::cout<<elem_f[0]<<std::endl;
        

        std::vector<Point> slave_points_reference;
        ciao.clear();
        if ( elem_f[0].has_children() )
        {
            for (int i=0; i<elem_f[0].n_children(); ++i)
            {
                Elem const * childOfFine=elem_f[0].child_ptr(i);
                //std::cout<<childOfFine[0]<<std::endl;
                for (int fine_node=0; fine_node<childOfFine[0].n_nodes(); ++fine_node)
                {
                    Node const fineNode=childOfFine[0].node_ref(fine_node);
                    dof_id_type local_dof_id=fineNode.dof_number(0, 0, 0);
                    Point const finePoint(fineNode);
                    ciao.insert(std::pair<dof_id_type,Point>(local_dof_id,finePoint) );
                }
            }
        }
        else
        {
        }
        
        std::map<dof_id_type,Point>::iterator it = ciao.begin();
        std::map<dof_id_type,Point>::iterator const endIt = ciao.end();
        
        std::vector<dof_id_type> dof_indices_f;
        
        for ( ; it!=endIt; ++it)
        {
            slave_points_reference.push_back( myfe_c.inverse_map (elem_c,it->second) );
            dof_indices_f.push_back(it->first);
            //std::cout<<"verifica"<<slave_points_reference.at(id)<<std::endl;
        }
        
        //for (int i=0; i<slave_points_reference.size(); ++i)
        //{
        //    std::cout<<i<<" "<<slave_points_reference.at(i)<<std::endl<<std::endl<<std::endl;
        //}
        
        std::vector<dof_id_type> dof_indices_c;
        dof_map_c.dof_indices (elem_c, dof_indices_c);
        
        fe_c->reinit(elem_c,&slave_points_reference);
        std::cout<<"phi_c.size() = "<<phi_c.size()<<std::endl;
                        for (int i=0; i<phi_c.size(); ++i)
                        {
                            std::cout<<"phi_c.at("<<i<<").size()="<<phi_c.at(i).size()<<std::endl;
                            for (int j=0; j<phi_c.at(i).size(); ++j)
                            {
                                std::cout<<phi_c.at(i).at(j)<<" ";
                                
                                proj_matrix.set( dof_indices_f.at(j) , dof_indices_c.at(i) ,phi_c.at(i).at(j));
                                
                            }
                            std::cout<<std::endl;
                        }
        
    }
    proj_matrix.close();
    proj_matrix.print_matlab("projectionMatrix.m");
    
//    MeshBase::const_node_iterator it = coarse_mesh.local_nodes_begin();
//    //MeshBase::const_node_iterator it = coarse_mesh.active_nodes_begin();
//    const MeshBase::const_node_iterator it_end = coarse_mesh.local_nodes_end();
//    //const MeshBase::const_node_iterator it_end = coarse_mesh.active_nodes_end();
//    for ( ; it != it_end; ++it)
//    {
//        const Node * node = *it;
//        dof_id_type const id=node[0].id();
//
//        if ( hanging_nodes.find(id)==hanging_nodes.end() )
//        {
//            reorganization_matrix.set(id,id,1.0);
//        }
//        else
//        {
//            std::map<dof_id_type, std::vector<dof_id_type> >::iterator hang_it=hanging_nodes.find(id);
//            //            std::cout<<hang_it->second.size()<<std::endl;
//            //            std::cout<<id<<" "<<hang_it->first<<std::endl;
//
//            for (int i=0;i<hang_it->second.size(); ++i)
//            {
//                dof_id_type hang_id=hang_it->second.at(i);
//                reorganization_matrix.set(hang_id,id,0.5);
//            }
//        }
//
//    }
    
    
//    std::vector<dof_id_type>   dof_indices;
//
//    MeshBase::const_element_iterator el=_mesh[0].level_elements_begin(0);
//    MeshBase::const_element_iterator end_el=_mesh[0].level_elements_end(0);
//
//    int count=0;
//    for ( ; el != end_el ; ++el)
//    {
//        ++count;
//        std::cout<<"count="<<count<<std::endl;
//        const Elem * const elem=*el;
//
//        _dof_map.dof_indices (elem, dof_indices);
//
//        for (int i=0; i<dof_indices.size(); ++i)
//        {
//            std::cout<<i<<" "<<dof_indices.at(i)<<std::endl;
//        }
//        std::cout<<std::endl<<std::endl;
//
//        for (int i=0; i< elem[0].n_nodes(); ++i)
//        {
//            Node const & node=elem[0].node_ref(i);
//            //Point coarsePoint=static_cast<const Point>(coarseNode);
//            _dof_map.dof_indices (&node, dof_indices);
//
//            std::cout<<"i="<<i<<" id=" <<node.id()<<" dof="<<dof_indices.at(0)<<std::endl;
//            std::cout<<node<<std::endl;
//        }
//
//
//
//    }
//
//    el=_mesh[0].level_elements_begin(1);
//    end_el=_mesh[0].level_elements_end(1);
//
//    count=0;
//    for ( ; el != end_el ; ++el)
//    {
//        ++count;
//        std::cout<<"count="<<count<<std::endl;
//
//        const Elem * const elem=*el;
//        _dof_map.dof_indices (elem, dof_indices);
//
//        for (int i=0; i<dof_indices.size(); ++i)
//        {
//            std::cout<<i<<" "<<dof_indices.at(i)<<std::endl;
//        }
//        std::cout<<std::endl<<std::endl;
//    }
    

    //
    //    el=mesh2.level_elements_begin(0);
    //    end_el=mesh2.level_elements_end(0);
    //    for ( ; el != end_el ; ++el)
    //    {
    //        const Elem * const elem=*el;
    //        std::cout<<elem[0]<<std::endl;
    //    }
    //mesh2.update_parallel_id_counts();
    
    
    
    //Mesh *slave_mesh=static_cast<Mesh *>(_mesh);
    
    //DistributedMesh distMesh(*static_cast<DistributedMesh *>(_mesh);
    //__mesh[0].print_info();
    //    EquationSystems equation_systems (_mesh);
    
    
    // if(_meshes.size() == 0)
    //   mooseError("You need to provide meshes or use uniform refinement.");
    //
    
    // problem.setNumberOfLevels(_mesh_files.size()+1);
    // std::cout<<_mesh_files.size()<<"_mesh_files.size()"<<std::endl;
    //
    // DistributedMesh master_mesh(*static_cast<DistributedMesh *>(&_fe_problem->es().get_mesh()));
    // DistributedMesh slave_mesh(*static_cast<DistributedMesh *>(&_fe_problem->es().get_mesh()));
    //
    //   // loop over the levels creating operators from fine to coarse
    //   for(unsigned int i =_mesh_files.size(); i>0 ; --i){
    //
    //     if  (i==_mesh_files.size()){
    //
    //         master_mesh.read(_mesh_files[0]);
    //
    //         if (_primal_variable_order == "SECOND" || _dual_variable_order == "SECOND")
    //         {
    //           master_mesh.all_second_order(true);
    //           slave_mesh.all_second_order(true);
    //         }
    //
    //
    //         master_mesh.print_info();
    //         slave_mesh.print_info();
    //
    //         EquationSystems master_es (master_mesh);
    //         EquationSystems slave_es (slave_mesh);
    //         buildTransferOperators(/*slave_mesh*/slave_es,/*master_mesh*/master_es, 1, problem, _mesh_files.size() - i);
    //
    //
    //
    //       }
    //     else{
    //
    //         master_mesh.read(_mesh_files[i]);
    //         slave_mesh.read(_mesh_files[i-1]);
    //
    //        if (_primal_variable_order == "SECOND" || _dual_variable_order == "SECOND"){
    //             master_mesh.all_second_order(true);
    //             slave_mesh.all_second_order(true);
    //           }
    //         master_mesh.print_info();
    //         slave_mesh.print_info();
    //
    //         EquationSystems master_es (master_mesh);
    //         EquationSystems slave_es (slave_mesh);
    //         buildTransferOperators(/*slave_mesh*/slave_es,/*master_mesh*/master_es, 0, problem,_mesh_files.size() - i);
    //
    //     }
    //
    //
    //   }
    std::cout<<"fine\n";
    exit(1);
}

void Interpolator::verifyMaps()
{
    for (int i=0; i<_coarse2fine.size(); ++i)
    {
        std::cout<<"The number of element on level "<<i<<" is "<<_coarse2fine.at(i).size()<<std::endl;
        std::map<Elem *,Elem *>::iterator it = _coarse2fine.at(i).begin();
        std::map<Elem *,Elem *>::iterator endIt = _coarse2fine.at(i).end();
        for (; it!=endIt; ++it)
        {
            Elem * elem_c=it->first;
            Elem * elem_f=it->second;
            std::cout<<"Id c = "<<elem_c[0].id()<<" Id f = "<<elem_f[0].id()<<std::endl;
            if (elem_c==elem_f)
            {
                std::cout<<"the pointers are the same! error\n";
                exit(1);
            }
        }

    }
    
}

void Interpolator::createMeshOnLevel(int level)
{
    _meshes[level]=new Mesh(_meshes[level+1][0]);
    
    // Here we try to create a local to glocal map
    MeshBase::const_element_iterator  el_fine     = _meshes[level+1][0].active_elements_begin();
    MeshBase::const_element_iterator  end_el_fine = _meshes[level+1][0].active_elements_end();

    MeshBase::const_element_iterator  el_coarse     = _meshes[level][0].active_elements_begin();
    MeshBase::const_element_iterator  end_el_coarse = _meshes[level][0].active_elements_end();

    for ( ; el_coarse != end_el_coarse ; ++el_coarse)
    {
        Elem * const elem_c=*el_coarse;
        Elem * const elem_f=*el_fine;
        
        // check if pointers are different
        if (elem_c==elem_f)
        {
            std::cout<<"error\n";
            exit(1);
        }

        // check if centers are different
        RealVectorValue c_c=elem_c[0].centroid();
        RealVectorValue c_f=elem_f[0].centroid();
        RealVectorValue diff=c_c-c_f;
        if (diff.norm()>1e-15)
        {
            std::cout<<"the centers are not the same "<<diff.norm()<<"\n";
            exit(1);
        }
        
        if (elem_c[0].level()==level+1)
        {
            if (elem_f[0].level()!=level+1)
            {
                std::cout<<"the element is the same but the level is different => impossible\n";
                exit(1);
            }
            Elem * elem_c_p=elem_c[0].parent();
            Elem * elem_f_p=elem_f[0].parent();
            // check if pointers are different
            if (elem_c_p==elem_f_p)
            {
                std::cout<<"error\n";
                exit(1);
            }
            _coarse2fine.at(level).insert ( std::pair<Elem *,Elem *>(elem_c_p,elem_f_p) );
            _fine2coarse.at(level).insert ( std::pair<Elem *,Elem *>(elem_f_p,elem_c_p) );
        }
        else
        {
            _coarse2fine.at(level).insert ( std::pair<Elem *,Elem *>(elem_c,elem_f) );
            _fine2coarse.at(level).insert ( std::pair<Elem *,Elem *>(elem_f,elem_c) );
        }
        
        ++el_fine;
    }
    
    MeshRefinement mesh_refinement(_meshes[level][0]);
    
    MeshBase::const_element_iterator el=_meshes[level][0].level_elements_begin(level+1);
    MeshBase::const_element_iterator end_el=_meshes[level][0].level_elements_end(level+1);
    
    for ( ; el != end_el ; ++el)
    {
        Elem * const elem=*el;
        elem[0].set_refinement_flag(Elem::COARSEN);
        //Elem::RefinementState a=elem[0].refinement_flag();
    }
    mesh_refinement.refine_and_coarsen_elements();
    _mesh[0].contract();

    el=_meshes[level][0].level_elements_begin(level+1);
    for ( ; el != end_el ; ++el)
    {
        Elem * const elem=*el;
        _meshes[level][0].delete_elem(elem);
    }
    _meshes[level][0].prepare_for_use();
    
    std::map<dof_id_type, std::vector<dof_id_type> > hanging_nodes;
    libMesh::MeshTools::find_hanging_nodes_and_parents(_meshes[level][0],hanging_nodes);
    
    std::cout<<"Number of hanging nodes on level "<<level<<" is "<<hanging_nodes.size()<<std::endl;
    
}

// #undef __FUNCT__
// #define __FUNCT__ "buildTransferOperators"
// void
// Interpolator::buildTransferOperators(EquationSystems &slave_es, EquationSystems &master_es, bool _apply_bc, HartFEProblem & problem, int level)
// {
//
//     PetscFunctionBegin;
//     using namespace libMesh;
//     utopia::DSMatrixd _B1 ;
//
//     Mat _B, Inter;
//
//     MeshBase & slave_mesh = slave_es.get_mesh();
//     MeshBase & master_mesh = master_es.get_mesh();
//
//     std::string pv_string;
//     std::string dv_string;
//
//     if (!slave_es.has_system("TMP_slave"))
//     {
//
//       slave_es.add_system<NonlinearImplicitSystem> ("TMP_slave");
//
//
//       for (int iName=0; iName<_primal_variable_name.size(); ++iName)
//       {
//           int temp=slave_es.get_system("TMP_slave").add_variable(_primal_variable_name.at(iName), FEType_primal);
//           variableNumberSlave.push_back(temp);
//       }
//
//       for (int iName=0; iName<_dual_variable_name.size(); ++iName)
//       {
//         int temp=slave_es.get_system("TMP_slave").add_variable(_dual_variable_name.at(iName), FEType_dual);
//         variableNumberSlave.push_back(temp);
//       }
//       slave_es.init();
//     }
//
//
//     // coarse mesh
//     if (!master_es.has_system("TMP_master")){
//
//       master_es.add_system<NonlinearImplicitSystem> ("TMP_master");
//
//       for (int iName=0; iName<_primal_variable_name.size(); ++iName)
//       {
//         int temp=master_es.get_system("TMP_master").add_variable(_primal_variable_name.at(iName), FEType_primal);
//         variableNumberMaster.push_back(temp);
//
//         if (pv_string.size()!=0)
//           pv_string+=" ";
//         pv_string+=_primal_variable_name.at(iName);
//
//       }
//
//       for (int iName=0; iName<_dual_variable_name.size(); ++iName)
//       {
//         int temp=master_es.get_system("TMP_master").add_variable(_dual_variable_name.at(iName), FEType_dual);
//         variableNumberMaster.push_back(temp);
//
//         if (dv_string.size()!=0)
//           dv_string+=" ";
//         dv_string+=_dual_variable_name.at(iName);
//       }
//       master_es.init();
//
//   }
//
//
//     DofMap &slave_dof  = slave_es.get_system("TMP_slave").get_dof_map();
//     DofMap &master_dof = master_es.get_system("TMP_master").get_dof_map();
//
//     moonolith::Communicator comm(_communicator.get());
//     PetscPrintf(PETSC_COMM_WORLD, "InterpolatorUniform: Info master mesh: \n");
//     master_mesh.print_info();
//     PetscPrintf(PETSC_COMM_WORLD, "InterpolatorUniform: Info slave mesh: \n");
//     slave_mesh.print_info();
//     std::chrono::high_resolution_clock::time_point _t_start_interpol = std::chrono::high_resolution_clock::now();
//
//
//    //   IS  idx_disp;
//    // is_from_variables(&idx_disp, "disp_x disp_y disp_z", & slave_es.get_system("TMP_slave"));
//    //       std::cout<<" ======================================"<<std::endl;
//    //  std::cout<<" L2 idx_disp"<<std::endl;
//
//    //  vi(idx_disp);
//
//    //   std::cout<<" ======================================"<<std::endl;
//
//     for (int kk=0; kk<variableNumberSlave.size(); ++kk)
//     {
//
//     Mat temp;
//
//
//     assemble_volume_transfer(comm,
//     utopia::make_ref(master_mesh),
//     utopia::make_ref(slave_mesh),
//     utopia::make_ref(master_dof),
//     utopia::make_ref(slave_dof),
//     variableNumberMaster.at(kk),
//     variableNumberSlave.at(kk),
//     false,
//     1,
//     _B1,
//     {},
//     true);
//
//
//
//     utopia::Interpolator(make_ref(_B1)).describe(std::cout);
//
//
//     MatDuplicate(raw_type(_B1), MAT_COPY_VALUES, &temp);
//   if(kk==0)
//   {
//     MatDuplicate(temp, MAT_COPY_VALUES, &_B);
//     MatDestroy(&temp);
//
//     // MatSetOption(_B, MAT_NEW_NONZERO_LOCATIONS,   PETSC_TRUE) ;
//     // MatSetOption(_B, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE) ;
//     //  MatSetOption(_B, MAT_NO_OFF_PROC_ENTRIES,     PETSC_FALSE) ;
//
//
//   }else
//   {
//     MatAXPY(_B, 1.0, temp,  DIFFERENT_NONZERO_PATTERN);
//     MatDestroy(&temp);
//    }
//
//
//     }
//
//     MatAssemblyBegin(_B,MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(_B,MAT_FINAL_ASSEMBLY);
//
//
//     MatDuplicate(_B, MAT_COPY_VALUES,&Inter);
//
//     Vec onesLarge;
//     PetscInt r,c,R,C;
//
//     MatGetLocalSize(Inter,&r,&c);
//     MatGetSize(Inter,&R,&C);
//
//     VecCreate(PETSC_COMM_WORLD,&onesLarge);
//     VecSetType(onesLarge, VECMPI);
//     VecSetSizes(onesLarge, r,R);
//     MatGetRowSum(Inter, onesLarge);
//
//     PetscScalar onesLargeMax;
//     PetscScalar onesLargeMin;
//
//     VecMax(onesLarge,NULL,&onesLargeMax);
//     VecMin(onesLarge,NULL,&onesLargeMin);
//
//     std::cout<<"max ="<<onesLargeMax<<std::endl;
//     std::cout<<"min ="<<onesLargeMin<<std::endl;
//
//     VecReciprocal(onesLarge);
//     MatDiagonalScale(Inter,onesLarge,NULL);
//     VecDestroy(&onesLarge);
//
//     MatAssemblyBegin(Inter,MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(Inter,MAT_FINAL_ASSEMBLY);
//
//     MatDestroy(&_B);
//
//     if(_apply_bc == true)
//      apply_BC_to_interpolation(Inter, _dc_boundary_id);
//
//    MatDuplicate(Inter,MAT_COPY_VALUES, &(problem.interpol_store[level]));
//    PetscPrintf(PETSC_COMM_WORLD, "\n Level  =%i \n", level);
//
//
//   MatGetSize(Inter, &r, &c);
//   PetscPrintf(PETSC_COMM_WORLD, "\nGlobal dimensions of extracted second order interpolator N=%i, M=%i. \n", c,r);
//
//
//   //sm(Inter,"matInter.m","I");
//
//  if(_primal_variable_name.size() > 0)
//   {
//
//
//     IS primal_master, primal_slave;
//
//     is_from_variables(&primal_slave , pv_string,  &slave_es.get_system("TMP_slave") );
//     is_from_variables(&primal_master, pv_string, &master_es.get_system("TMP_master"));
//
//     Mat A_primal;
//     MatGetSubMatrix(Inter,primal_slave,primal_master,MAT_INITIAL_MATRIX,&A_primal);
//     MatDuplicate(A_primal,MAT_COPY_VALUES,  &(problem.interpol_store_primal[level])  );
//     MatGetSize(A_primal, &r, &c);
//     PetscPrintf(PETSC_COMM_WORLD, "\nGlobal dimensions of extracted second order interpolator A_primal N=%i, M=%i. \n", c,r);
//     MatDestroy(&A_primal);
//     ISDestroy(&primal_slave);
//     ISDestroy(&primal_master);
//
//   }
//
//   if(_dual_variable_name.size() > 0)
//   {
//     IS dual_master,dual_slave;
//
//     is_from_variables(&dual_slave , dv_string,  &slave_es.get_system("TMP_slave") );
//     is_from_variables(&dual_master, dv_string, &master_es.get_system("TMP_master"));
//
//     Mat A_dual;
//     MatGetSubMatrix(Inter,dual_slave,dual_master,MAT_INITIAL_MATRIX,&A_dual);
//     MatDuplicate(A_dual,MAT_COPY_VALUES,  &(problem.interpol_store_dual[level])  );
//     MatGetSize(A_dual, &r, &c);
//     PetscPrintf(PETSC_COMM_WORLD, "\nGlobal dimensions of extracted second order interpolator A_dual N=%i, M=%i. \n", c,r);
//     MatDestroy(&A_dual);
//     ISDestroy(&dual_slave);
//     ISDestroy(&dual_master);
//
//
//   }
//
//   MatDestroy(&Inter);
//
//   variableNumberSlave.erase(variableNumberSlave.begin(),variableNumberSlave.end());
//   variableNumberMaster.erase(variableNumberMaster.begin(),variableNumberMaster.end());
//    // std::chrono::high_resolution_clock::time_point _t_end = std::chrono::high_resolution_clock::now();
//    // std::chrono::duration<double, std::milli> fp_ms = _t_end - _t_start_interpol ;
//    // PetscPrintf(PETSC_COMM_WORLD, "\nInterpolator computed in %f ms. \n\n", fp_ms);
//
//     PetscFunctionReturnVoid();
//}


/**
 * @brief      Function to zero rows in interpolation operator related to D-BC
 *
 * @param      I                Interpolation operator
 * @param[in]  _dc_boundary_id  IDs of active boundaries
 */
// #undef __FUNCT__
// #define __FUNCT__ "apply_BC_to_interpolation"
//  void Interpolator::apply_BC_to_interpolation(Mat & I, const std::vector<int> _dc_boundary_id)
//  {
//   std::chrono::high_resolution_clock::time_point _t_start = std::chrono::high_resolution_clock::now();
//   PetscFunctionBegin;
//
//   PetscErrorCode ierr;
//   std::vector<int> zero_rows;
//   MooseMesh * _mesh = &(_fe_problem->mesh());
//   ConstBndNodeRange & bnd_nodes = *_mesh->getBoundaryNodeRange();
//
//   unsigned int i = 0;
//   for(auto boundary = _dc_boundary_id.begin(); boundary != _dc_boundary_id.end(); ++boundary, i++)
//   {
//
//
//     // iterate just over boundary nodes
//     for(const auto & bnode : bnd_nodes)
//     {
//       libMesh::Node * current_node = bnode->_node;
//
//       // check if node is in active boundary list
//       if(_mesh->isBoundaryNode(current_node->id(), *boundary))
//       {
//         // loop over all variables at this node
//         for(auto v = 0; v < _primal_variable_name.size(); v++) //change
//         {
//           const Variable & var  = _nl.sys().variable(v);
//           unsigned int var_num  = var.number();
//
//           // see if this variable has any dofs at this node
//           if(current_node->n_dofs(_nl.number(), var_num) > 0)
//           {
//             // check if given variable has BC on node
//             if(std::find(_dc_variables_id[i].begin(), _dc_variables_id[i].end(), var_num) != _dc_variables_id[i].end())
//             {
//               // different components are not supported by moose at the moment...
//               zero_rows.push_back(current_node->dof_number(_nl.number(), var_num, 0));
//
//               std::cout<<"current_node->dof_number(_nl.number(), var_num, 0): "<< current_node->dof_number(_nl.number(), var_num, 0) << "  \n";
//             }
//           }
//         }
//       }
//     }
//   }
//
//   IS i_dc_nodes;
//   ierr = ISCreateGeneral(PETSC_COMM_WORLD, zero_rows.size(), &zero_rows[0], PETSC_USE_POINTER, &i_dc_nodes); CHKERRV(ierr);
//
//   ierr = MatZeroRows(I, zero_rows.size(), &zero_rows[0], 0, NULL, NULL); CHKERRV(ierr);
//   zero_rows.clear();
//   ierr = ISDestroy(&i_dc_nodes); CHKERRV(ierr);
//
//   std::chrono::high_resolution_clock::time_point _t_end = std::chrono::high_resolution_clock::now();
//   std::chrono::duration<double, std::milli> fp_ms = _t_end - _t_start ;
//
//   std::cout<<"  Walltime of BC_TO interpolation: " << std::to_string(fp_ms.count()) << "  ms. \n";
//
//
//   PetscFunctionReturnVoid();
// }

// #undef __FUNCT__
// #define __FUNCT__ "determine_dc_bnd_var_id"
// void Interpolator::determine_dc_bnd_var_id(const std::vector<std::string> & BC_var)
// {
//   // automatic fill-in
//   std::vector<int> vec(_nVars);
//   std::iota(vec.begin(), vec.end(), 0);
//
//   unsigned int i;
//   auto str_tmp = BC_var.begin();
//
//   PetscFunctionBegin;
//   // going over all BC_ids
//   for(i = 0, str_tmp; str_tmp != BC_var.end(); i++, str_tmp++)
//   {
//     std::vector<std::string> tmp = Interpolator::split(*str_tmp, '-');
//
//     // check if variable assigned in the input file exists for given simulation
//     bool var_flg = 1;
//     for(auto t = tmp.begin(); t != tmp.end(); ++t)
//     {
//       if(atoi(t->c_str()) >= _nVars)
//         var_flg = 0;
//     }
//
//     // in case u havent put anything into input file, or u put too much
//     if(*str_tmp == "-1" || var_flg == 0)
//     {
//       _dc_variables_id.push_back(vec);
//     }
//     else
//     {
//
//       unsigned int j;
//       std::vector<int > one_BC_id;
//       auto str_in = tmp.begin();
//       for(j = 0, str_in; str_in != tmp.end(); j ++, str_in++)
//       {
//         one_BC_id.push_back(atoi(str_in->c_str()));
//       }
//       _dc_variables_id.push_back(one_BC_id);
//     }
//   }
//
//   // check if u have same number of BC_ids in both parameters
//   if(_dc_variables_id.size() != _dc_boundary_id.size())
//   {
//     _dc_variables_id.clear();
//     for(auto i = 0; i != _dc_boundary_id.size(); i++)
//     {
//       _dc_variables_id.push_back(vec);
//     }
//   }
//
//   // print out what is considered for zero-ing
//   std::cout<<" ------ BC CONDITIONS  ------ \n";
//   unsigned int t = 0;
//   for(auto i = _dc_variables_id.begin(); i != _dc_variables_id.end();  t++, i++)
//   {
//     std::cout<<"\n BC_id:  "<< _dc_boundary_id[t] << "   var_ids:  ";
//     std::for_each(i->begin(), i->end(), [](int i){ std::cout << i << "  " ; });
//     std::cout<<"\n"<<std::endl;
//   }
//   PetscFunctionReturnVoid();
// }

Interpolator::~Interpolator()
{
    for (int i=0; i<_numOfLevels; ++i)
        delete _meshes[i];
    //delete [] _meshes;
}


/////////////////////////////////////////////////////////////////////////////////
////////////////////////
////////////////////////    KEEP THIS
////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// This is a useful code to create a mesh but it is not useful for our purposes
//void Interpolator::createMeshOnLevel(int level)
//{
//    _meshes[level]=new Mesh(_mesh[0].comm(),_meshDimension);
//    _meshes[level][0].clear();
//    _meshes[level][0].set_mesh_dimension(_meshDimension);
//    _meshes[level][0].set_spatial_dimension(_spatialDimension);
//
//    _global2local_nodes.at(level).resize(_mesh[0].n_nodes(),-999);
//    _local2global_nodes.at(level).resize(0);
//    _global2local_elem.at(level).resize(_mesh[0].n_elem(),-999);
//    _local2global_elem.at(level).resize(0);
//
//    int elemCount=0;
//    int nodeCount=0;
//
//    for (int l=level; l>=0; --l)
//    {
//        MeshBase::const_element_iterator el=_mesh[0].level_elements_begin(l);
//        MeshBase::const_element_iterator end_el=_mesh[0].level_elements_end(l);
//
//        for ( ; el != end_el ; ++el)
//        {
//            const Elem * const elem=*el;
//
//            if (elem[0].level()==level || ( elem[0].level()<level && !elem[0].has_children() ) ) //
//            {
//                //std::cout<<elem[0].level()<<" "<<elem[0].id()<<std::endl;
//                if ( _global2local_elem.at(level).at(elem->id())==-999 )
//                {
//                    _global2local_elem.at(level).at(elem->id())=elemCount;
//                    _local2global_elem.at(level).push_back(elem->id());
//                }
//                else
//                {
//                    std::cout<<"you should not be here\n";
//                    exit(1);
//                }
//
//                for (int i=0; i< elem[0].n_nodes(); ++i)
//                {
//                    Node const & node=elem[0].node_ref(i);
//                    if (_global2local_nodes.at(level).at(node.id())==-999)
//                    {
//                        // This not has not been added
//                        _global2local_nodes.at(level).at(node.id())=nodeCount;
//                        _local2global_nodes.at(level).push_back(node.id());
//                        Point point=static_cast<const Point>(node);
//                        _meshes[level][0].add_point(point,nodeCount);
//                        ++nodeCount;
//                    }
//
//                }
//            }
//
//            Elem * elemToAdd = new Quad4;
//            elemToAdd->set_id(elemCount);
//            elemToAdd = _meshes[level][0].add_elem (elemToAdd);
//            for (int i=0; i< elem[0].n_nodes(); ++i)
//            {
//                Node const & node=elem[0].node_ref(i);
//                elemToAdd->set_node(i) = _meshes[level][0].node_ptr( _global2local_nodes.at(level).at(node.id()) );
//            }
//            ++elemCount;
//        }
//    }
//
//    _meshes[level][0].prepare_for_use (/*skip_renumber =*/ false);
//    _meshes[level][0].update_parallel_id_counts();
//    _meshes[level][0].prepare_for_use (/*skip_renumber =*/ false);
//}
