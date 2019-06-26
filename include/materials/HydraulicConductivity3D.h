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

#ifndef HydraulicConductivity3D_H
#define HydraulicConductivity3D_H

#include "Material.h"

//Forward Declarations
class HydraulicConductivity3D;

template<>
InputParameters validParams<HydraulicConductivity3D>();

/**
 * Example material class that defines a few properties.
 */
class HydraulicConductivity3D : public Material
{
public:
  HydraulicConductivity3D(const InputParameters & parameters);
    
    void ComputeNormalsFromAngles(
                                  RealVectorValue const & angles,
                                  RealVectorValue & n1,
                                  RealVectorValue & n2,
                                  RealVectorValue & n3);

    int is_inside(RealVectorValue const & point);
    
  ~HydraulicConductivity3D()
    {
//        delete [] _fx;
//        delete [] _fy;
//        delete [] _ft;
//        delete [] _fl;
//        delete [] _fa;
//        delete [] _a;
//        delete [] _b;
//        delete [] _c;
//        delete [] _ao;
//        delete [] _bo;
//        delete [] _co;
//
        std::cout<<"Called desctructor\n";
    };

protected:
  virtual void computeQpProperties();
//  virtual void initQpStatefulProperties();

    void outerProduct
    (RealVectorValue const & in1, RealVectorValue const &  in2, RealTensorValue & out);
    int findRegion(RealVectorValue const & point,std::vector<int> & in);
    
    int _fn;
    
    int _prec;
    
    Real _min_dimension;
    
    std::string _fx_string;
    std::string _fy_string;
    std::string _fz_string;
    std::string _fa1_string;
    std::string _fa2_string;
    std::string _fa3_string;
    std::string _fd1_string;
    std::string _fd2_string;
    std::string _fd3_string;
    
    RealVectorValue * _center;
    RealVectorValue * _rotation;
    RealVectorValue * _dimension;
    
    RealVectorValue ** _n;
    
    RealVectorValue *_d;
    
    RealVectorValue * _max_f;
    RealVectorValue * _min_f;
    
    bool * shall_check;

    
    Real _pi;
    
    RealTensorValue _identity;
    

    MaterialProperty<RealTensorValue> &_K_filettata;
    MaterialProperty<Real> &_phi;
    Real _phiFracture,_phiMatrix, _K_matrix;
    bool _cond0, _cond1;
    const VariableGradient &_gradP;
    MaterialProperty<RealVectorValue> &_U;

    MaterialProperty<Real> &_numOfFrac;
    
    std::vector<int> _whichFrac;
    
    std::vector<RealVectorValue> _regionMin;
    std::vector<RealVectorValue> _regionMax;
    MaterialProperty<int> &_regionID;
    MaterialProperty<Real> &_regionIDReal;
};

#endif //FreqNavierStokesFracture2D_H
