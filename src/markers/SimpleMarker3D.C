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

#include "SimpleMarker3D.h"

// libMesh includes
#include "libmesh/error_vector.h"

registerMooseObject("parrot_realApp", SimpleMarker3D);

//#include <algorithm>

template<>
InputParameters validParams<SimpleMarker3D>()
{
    InputParameters params = validParams<Marker>();
    params.addRequiredParam<int>("fn", "number of fractures");
    params.addRequiredParam<std::string>("fx_string", "x-coordinates of center of fractures");
    params.addRequiredParam<std::string>("fy_string", "y-coordinates of center of fractures");
    params.addRequiredParam<std::string>("fz_string", "z-coordinates of center of fractures");
    params.addRequiredParam<std::string>("fa1_string", "rotation along z axis");
    params.addRequiredParam<std::string>("fa2_string", "rotation along y axis");
    params.addRequiredParam<std::string>("fa3_string", "rotation along x axis");
    params.addRequiredParam<std::string>("fd1_string", "fracture dimension 1");
    params.addRequiredParam<std::string>("fd2_string", "fracture dimension 2");
    params.addRequiredParam<std::string>("fd3_string", "fracture dimension 3");
    
    return params;
}


SimpleMarker3D::SimpleMarker3D(const InputParameters & parameters) :
Marker(parameters),
_fn(getParam<int>("fn")),
_fx_string(getParam<std::string>("fx_string")),
_fy_string(getParam<std::string>("fy_string")),
_fz_string(getParam<std::string>("fz_string")),
_fa1_string(getParam<std::string>("fa1_string")),
_fa2_string(getParam<std::string>("fa2_string")),
_fa3_string(getParam<std::string>("fa3_string")),
_fd1_string(getParam<std::string>("fd1_string")),
_fd2_string(getParam<std::string>("fd2_string")),
_fd3_string(getParam<std::string>("fd3_string"))
{
    _center   =new RealVectorValue [_fn];
    _rotation =new RealVectorValue [_fn];
    _dimension=new RealVectorValue [_fn];
    _d        =new RealVectorValue [_fn];
    
    shall_check= new bool [_fn];
    
    _n = new RealVectorValue * [_fn];
    for (int i=0; i<_fn; ++i)
    {
        _n[i]=new RealVectorValue [3];
    }
    
    _max_f=new RealVectorValue[_fn];
    _min_f=new RealVectorValue[_fn];
    
    std::istringstream  fx_ss( _fx_string);
    std::istringstream  fy_ss( _fy_string);
    std::istringstream  fz_ss( _fz_string);
    std::istringstream fa1_ss(_fa1_string);
    std::istringstream fa2_ss(_fa2_string);
    std::istringstream fa3_ss(_fa3_string);
    std::istringstream fd1_ss(_fd1_string);
    std::istringstream fd2_ss(_fd2_string);
    std::istringstream fd3_ss(_fd3_string);
    
    std::string token;
    
    for (int i=0; i<_fn; ++i)
    {
        std::getline(fx_ss, token, ',');
        _center[i](0)=std::atof(token.c_str());
        std::getline(fy_ss, token, ',');
        _center[i](1)=std::atof(token.c_str());
        std::getline(fz_ss, token, ',');
        _center[i](2)=std::atof(token.c_str());
        
        std::getline(fa1_ss, token, ',');
        _rotation[i](0)=atof(token.c_str());
        std::getline(fa2_ss, token, ',');
        _rotation[i](1)=atof(token.c_str());
        std::getline(fa3_ss, token, ',');
        _rotation[i](2)=atof(token.c_str());
        _rotation[i]=_rotation[i]/180.0*pi;
        
        std::getline(fd1_ss, token, ',');
        _dimension[i](0)=atof(token.c_str());
        std::getline(fd2_ss, token, ',');
        _dimension[i](1)=atof(token.c_str());
        std::getline(fd3_ss, token, ',');
        _dimension[i](2)=atof(token.c_str());
        
        ComputeNormalsFromAngles(_rotation[i],_n[i][0],_n[i][1],_n[i][2]);
        
//        std::cout<<_n[0][0]<<std::endl;
//        std::cout<<_n[0][1]<<std::endl;
//        std::cout<<_n[0][2]<<std::endl;
//        exit(1);
        
        for (int j=0; j<3; ++j)
        {
            _d[i](j)=_n[i][j]*_center[i];
        }
    }
    
    _min_dimension=1000000000;
    for (int i=0; i<_fn; ++i)
        for (int j=0; j<3; ++j)
            _min_dimension=std::min(_min_dimension,_dimension[i](j));
    
    
    _fx_string.clear();
    _fy_string.clear();
    _fz_string.clear();
    _fa1_string.clear();
    _fa2_string.clear();
    _fa3_string.clear();
    _fd1_string.clear();
    _fd2_string.clear();
    _fd3_string.clear();
    
    _prec = 10000000;
    srand (time(NULL));
    
    
    RealVectorValue * vertices;
    vertices=new RealVectorValue[8];
    
    for (int l=0; l<_fn; ++l)
    {
        
        
        int counter =0;
        for (int i=-1; i<=1; ++++i)
            for (int j=-1; j<=1; ++++j)
                for (int k=-1; k<=1; ++++k)
                {
                    vertices[counter]=_center[l]+
                    0.5*_dimension[l](0)*i*_n[l][0]+
                    0.5*_dimension[l](1)*j*_n[l][1]+
                    0.5*_dimension[l](2)*k*_n[l][2];
                    counter+=1;
                }
        
        _max_f[l]=vertices[0];
        _min_f[l]=vertices[0];
        for (int i=1; i<8; ++i)
        {
            for (int j=0; j<3; ++j)
            {
                _max_f[l](j)=std::max(_max_f[l](j), vertices[i](j));
                _min_f[l](j)=std::min(_min_f[l](j), vertices[i](j));
            }
        }
        
    }
    delete [] vertices;
    
    delete [] _center;
    delete [] _rotation;
    
}

void
SimpleMarker3D::markerSetup()
{
}

Marker::MarkerValue
SimpleMarker3D::computeElementMarker()
{
    RealVectorValue max_e;
    RealVectorValue min_e;
    
    max_e=(*_current_elem).point(0);
    min_e=(*_current_elem).point(0);
    for (int i = 1; i < (*_current_elem).n_nodes(); i++)
    {
        for (int j=0; j<3; ++j)
        {
            max_e(j)=std::max(max_e(j),(*_current_elem).point(i)(j));
            min_e(j)=std::min(min_e(j),(*_current_elem).point(i)(j));
        }
    }
    
    int counter=_fn;
    for (int l=0; l<_fn; ++l)
    {
        shall_check[l]=true;
        if (max_e(0)<_min_f[l](0) || _max_f[l](0)<min_e(0) )
        {
            shall_check[l]=false;
            counter-=1;
            continue;
        }
        else
        {
            if (max_e(1)<_min_f[l](1) || _max_f[l](1)<min_e(1) )
            {
                shall_check[l]=false;
                counter-=1;
                continue;
            }
            else
            {
                if (max_e(2)<_min_f[l](2) || _max_f[l](2)<min_e(2) )
                {
                    shall_check[l]=false;
                    counter-=1;
                    continue;
                }
                
            }
        }
    }

    
//    if (counter==0)
//    {
//        return DO_NOTHING;
//    }
    
    RealVectorValue _needed_points;
    _needed_points(0)=1.0;//(max_e(0)-min_e(0))/_min_dimension;
    _needed_points(1)=1.0;//(max_e(1)-min_e(1))/_min_dimension;
    _needed_points(2)=1.0;//(max_e(2)-min_e(2))/_min_dimension;
    
    std::cout<<"start   "<<std::flush;
    RealVectorValue point;
    for (int i=0; i<std::ceil( _needed_points(0) ); ++i )
        for (int j=0; j<std::ceil( _needed_points(0) ); ++j )
            for (int k=0; k<std::ceil( _needed_points(0) ); ++k )
            {
                for (int l=0;l<3;++l)
                    point(l)=(max_e(l)+min_e(l))/2.0;
                
                if ( is_inside(point))
                    return REFINE;
                
            }
    
    std::cout<<"stop"<<std::endl<<std::flush;
    
    return DO_NOTHING;
    
}


void SimpleMarker3D::ComputeNormalsFromAngles(RealVectorValue const & angles,
                                              RealVectorValue & n1,
                                              RealVectorValue & n2,
                                              RealVectorValue & n3)
{
    RealTensorValue R1;
    RealTensorValue R2;
    RealTensorValue R3;
    
    R1(0,0)=std::cos(angles(0));
    R1(0,1)=-std::sin(angles(0));
    R1(0,2)=0.0;
    R1(1,0)=std::sin(angles(0));
    R1(1,1)=std::cos(angles(0));
    R1(1,2)=0.0;
    R1(2,0)=0.0;
    R1(2,1)=0.0;
    R1(2,2)=1.0;
    
    R2(0,0)=std::cos(angles(1));
    R2(0,1)=0.0;
    R2(0,2)=-std::sin(angles(1));
    R2(1,0)=0.0;
    R2(1,1)=1.0;
    R2(1,2)=0.0;
    R2(2,0)=std::sin(angles(1));
    R2(2,1)=0.0;
    R2(2,2)=std::cos(angles(1));
    
    R3(0,0)=1.0;
    R3(0,1)=0.0;
    R3(0,2)=0.0;
    R3(1,0)=0.0;
    R3(1,1)=std::cos(angles(2));
    R3(1,2)=-std::sin(angles(2));
    R3(2,0)=0.0;
    R3(2,1)=std::sin(angles(2));
    R3(2,2)=std::cos(angles(2));
    
    //std::cout<<R1<<std::endl;
    //std::cout<<R2<<std::endl;
    //std::cout<<R3<<std::endl;
    //exit(1);
    RealTensorValue R=R1*R2*R3;
    std::cout<<R<<std::endl;
    
    for (int i=0; i<3; ++i)
    {
        n1(i)=R(i,0);
        n2(i)=R(i,1);
        n3(i)=R(i,2);
    }
//    std::cout<<n1<<std::endl;
//    std::cout<<n2<<std::endl;
//    std::cout<<n3<<std::endl;
    
}

bool SimpleMarker3D::is_inside(RealVectorValue const & point)
{
    std::cout<<"sono qui?\n";
    if (0.49 < point(0) && point(0)<0.51 )
    {
        std::cout<<point<<std::endl;
                exit(1);
    }

    for (int i=0; i<_fn; ++i)
    {
        //if (shall_check[i]==true)
        {
            Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
            if (temp1<=_dimension[i](0)/2.0)
            {
                Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
                if (temp2<=_dimension[i](1)/2.0)
                {
                    Real temp3=std::fabs( _n[i][2]*point-_d[i](2) );
                    if (temp3<=_dimension[i](2)/2.0)
                    {
                        return true;
                    }
                    
                }
            }
        }
    }
    return false;
}

SimpleMarker3D::~SimpleMarker3D()
{
    delete [] _max_f;
    delete [] _min_f;
    delete [] _dimension;
    for (int i=0; i<_fn; ++i)
    {
        delete [] _n[i];
    }
    delete [] _n;
    delete [] _d;
    
    delete [] shall_check;
    
}
