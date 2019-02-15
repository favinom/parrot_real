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

#include "HydraulicConductivity3D.h"

#include <sstream>
#include "MooseMesh.h"

registerMooseObject("parrot_realApp", HydraulicConductivity3D);

template<>
InputParameters validParams<HydraulicConductivity3D>()
{
    InputParameters params = validParams<Material>();
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
    params.addRequiredParam<std::vector<Real>>("conductivityMatrix","conductivity of the Matrix");
    params.addRequiredParam<std::vector<Real>>("conductivityFracture0","conductivity of intersection 0");
    params.addRequiredParam<std::vector<Real>>("conductivityFracture1","conductivity of intersection 1");
    params.addRequiredParam<std::vector<Real>>("conductivityFracture2","conductivity of the fracture");
    return params;
}

HydraulicConductivity3D::HydraulicConductivity3D(const InputParameters & parameters) :
Material(parameters),
_fn(getParam<int>("fn")),
_fx_string(getParam<std::string>("fx_string")),
_fy_string(getParam<std::string>("fy_string")),
_fz_string(getParam<std::string>("fz_string")),
_fa1_string(getParam<std::string>("fa1_string")),
_fa2_string(getParam<std::string>("fa2_string")),
_fa3_string(getParam<std::string>("fa3_string")),
_fd1_string(getParam<std::string>("fd1_string")),
_fd2_string(getParam<std::string>("fd2_string")),
_fd3_string(getParam<std::string>("fd3_string")),
_K(declareProperty<RealTensorValue>("conductivityTensor")),
_cond0(getParam<std::vector<Real>>("conductivityFracture0")),
_cond1(getParam<std::vector<Real>>("conductivityFracture1")),
_cond2(getParam<std::vector<Real>>("conductivityFracture2")),
_cond3(getParam<std::vector<Real>>("conductivityMatrix")),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"): _grad_zero),
_U(declareProperty<RealVectorValue>("VelocityVector")),
_numOfFrac(declareProperty<Real>("numero"))
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
    
    _identity=RealTensorValue(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    
}

void
HydraulicConductivity3D::computeQpProperties()
{
    RealVectorValue point=_q_point[_qp];
    
    //        for (int i=0; i<_fn; ++i)
    //        {
    //            std::cout<<i<<std::endl;
    //            std::cout<<_n[i][0]<<std::endl;
    //            std::cout<<_n[i][1]<<std::endl;
    //            std::cout<<_n[i][2]<<std::endl;
    //            std::cout<<_d[i](0)<<std::endl;
    //            std::cout<<_d[i](1)<<std::endl;
    //            std::cout<<_d[i](2)<<std::endl;
    //            std::cout<<_dimension[i](0)<<std::endl;
    //            std::cout<<_dimension[i](1)<<std::endl;
    //            std::cout<<_dimension[i](2)<<std::endl;
    //            std::cout<<point<<std::endl;
    //        }
    
    
    int count = is_inside(point);
    _numOfFrac[_qp]=count;
    
    if (count>0)
    {
        //We are inside at least a fracture
        _K[_qp]=1.0*_identity;
        
        if (count == 1)
        {
            // we are inside a single fracture.
            int q=_whichFrac.at(0);
            //std::cout<<"we are in Frac "<<q<<"con normale "<<_n[q][2]<<std::endl;
           _K[q]=_cond0[q]*_identity;
        }
        if (count == 2)
        {
            // we are inside on a line.
            int q0=_whichFrac.at(0);
            int q1=_whichFrac.at(1);
            //std::cout<<"we are on a line intesacting surface "<<q0<<" and "<<q1<<std::endl;
            //std::cout<<"normals are"<<_n[q0][2]<<" and "<<_n[q1][2]<<std::endl;
             _K[q0]=_cond1[q0]*_identity;
             _K[q1]=_cond1[q1]*_identity;
         }
        if (count == 3)
        {
            // we are inside on a dot.
            int q0=_whichFrac.at(0);
            int q1=_whichFrac.at(1);
            int q2=_whichFrac.at(2);
            //std::cout<<"we are on a dot intesacting surfaces "<<q0<<" and "<<q1<<" and "<<q2<<std::endl;
            //std::cout<<"normals are"<<_n[q0][2]<<", "<<_n[q1][2]<<" and " <<_n[q2][2]<<std::endl;
             _K[q0]=_cond2[q0]*_identity; 
             _K[q1]=_cond2[q1]*_identity;
             _K[q2]=_cond2[q2]*_identity;
       }
        if (count > 3)
        {
            std::cout<<"error\n";
            exit(1);
        }
    }
    else
    {
        //we are in background
        _K[_qp]=_cond3[_qp]*_identity;
    }

    
    //    Real x_coord = _q_point[_qp](0);
    //    Real y_coord = _q_point[_qp](1);
    //
    //
    //    for (int i = 0; i < _fn; i++)
    //    {
    //        Real temp =_a[i]*x_coord+_b[i]*y_coord + _c[i];
    //        if ( std::fabs(temp) <=  _ft[i]/2.0 )
    //        {
    //            Real tempo=_ao[i]*x_coord+_bo[i]*y_coord + _co[i];
    //            if (std::fabs(tempo) <=  _fl[i]/2.0 )
    //            {
    //                std::cout<<"dentro\n";
    //
    //            }
    //        }
    //    }
    //
    //    // Kinematics
    
}

void HydraulicConductivity3D::ComputeNormalsFromAngles(RealVectorValue const & angles,
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
    //std::cout<<R<<std::endl;
    
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

int HydraulicConductivity3D::is_inside(RealVectorValue const & point)
{
    int numOfFrac=0;
    _whichFrac.clear();
    for (int i=0; i<_fn; ++i)
    {
        
        Real temp1=std::fabs( _n[i][0]*point-_d[i](0) );
        if (temp1<_dimension[i](0)/2.0)
        {
            Real temp2=std::fabs( _n[i][1]*point-_d[i](1) );
            if (temp2<_dimension[i](1)/2.0)
            {
                Real temp3=std::fabs( _n[i][2]*point-_d[i](2) );
                if (temp3<_dimension[i](2)/2.0)
                {
                    numOfFrac+=1;
                    _whichFrac.push_back(i);
                }
                
            }
        }
        
    }
    return numOfFrac;
}
