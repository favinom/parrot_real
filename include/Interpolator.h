#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H

// MOOSE includes
#include "GeneralUserObject.h"
// #include "HartFEProblem.h"
// #include "ksp_hart_impl.h"
// #include "NonlinearSystem.h"

// #include "utilities.h"
// #include <utopia_fe.hpp>
// #include <utopia.hpp>
// #include <utopia_LibMeshBackend.hpp>

#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h"
// #include <libmesh/petsc_matrix.h>
// #include <libmesh/parallel_mesh.h>
// #include <LibmeshTransferForMoose.hpp>

// #include "MooseError.h"
// #include "MooseUtils.h"
// #include "MooseMesh.h"
// #include "petsc_libmesh_utilities.h"
// #include "projection_utilities.h"

// Forward declarations
class Interpolator;
// class HartFEProblem;

template<>
InputParameters validParams<Interpolator>();

/**
 * User object that reads takes an existing mesh coarsens it a specified
 * number of times, finds the transfer operator between meshes, and then sends
 * it onto passo.
 */
class Interpolator : public GeneralUserObject
{
public:
    Interpolator(const InputParameters & parameters);

    virtual ~Interpolator();
    virtual void initialize() override {std::cout<<"initialize\n";} ;
    virtual void finalize()   override {std::cout<<"finalize\n";} ;
    virtual void execute()    override    ;
    
    void verifyMaps();
    void createMeshOnLevel(int level);

protected:
    
    FEProblem       * _fe_problem;
    MooseMesh const & _mooseMesh;
    //MeshBase const * meshBase
    EquationSystems & _equationSystem;
    // Why does not it work if I put const Mesh ?
    Mesh            * _mesh;
    LinearImplicitSystem const & _lis;
    DofMap const & _dof_map;
    
    int _numOfLevels;
    int _meshDimension;
    int _spatialDimension;
    
    Mesh            ** _meshes;
    
    // test unordered_map to improve speed?
    std::vector< std::map<Elem *,Elem *> > _fine2coarse;
    std::vector< std::map<Elem *,Elem *> > _coarse2fine;
    
    EquationSystems ** _equationSystems;
    

    // KSP_HART * _ksp_p;
    
    // NonlinearSystem & _nl;
    // unsigned int _nVars;
    //int _mesh_dimension;
    // std::vector<int> _dc_boundary_id;
    // std::vector<std::vector<int> > _dc_variables_id;
    

    // std::string _primal_variable_order;
    // std::string _dual_variable_order;
    //
    // std::vector<int> variableNumberSlave;
    // std::vector<int> variableNumberMaster;

    // virtual void apply_BC_to_interpolation(Mat & , const std::vector<int>);
    // virtual void determine_dc_bnd_var_id(const std::vector<std::string> & );
    // virtual void buildTransferOperators(EquationSystems &slave_es,EquationSystems &master_es, bool _apply_bc, HartFEProblem &problem, int level);

//private:



};

#endif
