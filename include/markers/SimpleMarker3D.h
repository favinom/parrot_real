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

#ifndef SIMPLEMARKER3D_H
#define SIMPLEMARKER3D_H

#include"Marker.h"

class SimpleMarker3D;

template<>
InputParameters validParams<SimpleMarker3D>();

class SimpleMarker3D : public Marker
{
public:
    SimpleMarker3D(const InputParameters & parameters);
    virtual ~SimpleMarker3D();
    
    virtual void markerSetup();
    
protected:
    virtual MarkerValue computeElementMarker();
    
    void ComputeNormalsFromAngles(RealVectorValue const & angles,
                                  RealVectorValue & n1,
                                  RealVectorValue & n2,
                                  RealVectorValue & n3);

    bool is_inside(RealVectorValue const & point);
    
private:
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
    
};

#endif /* SIMPLEMARKER3D_H */
