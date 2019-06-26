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

#include "HydraulicConductivity2D.h"

#include <sstream>
#include "MooseMesh.h"

registerMooseObject("parrot_realApp", HydraulicConductivity2D);

template<>
InputParameters validParams<HydraulicConductivity2D>()
{
    InputParameters params = validParams<Material>();
    params.addRequiredParam<int>("fn", "number of fractures");
    params.addRequiredParam<std::string>("fx_string", "x-coordinates of center of fractures");
    params.addRequiredParam<std::string>("fy_string", "y-coordinates of center of fractures");
    params.addRequiredParam<std::string>("fa1_string", "rotation along z axis");
    params.addRequiredParam<std::string>("fd1_string", "fracture dimension 1");
    params.addRequiredParam<std::string>("fd2_string", "fracture dimension 2");
    params.addCoupledVar("pressure",
                         "The gradient of this variable will be used as "
                         "the velocity vector.");
    params.addRequiredParam<Real>("K_matrix","K_matrix");
    params.addRequiredParam<bool>("cond0","condition0");
    params.addParam<bool>("cond1","false","condition1");
    params.addRequiredParam<Real>("phi_f","phi fracture");
    params.addRequiredParam<Real>("phi_m","phi matrix");
    return params;
}

HydraulicConductivity2D::HydraulicConductivity2D(const InputParameters & parameters) :
Material(parameters),
_fn(getParam<int>("fn")),
_fx_string(getParam<std::string>("fx_string")),
_fy_string(getParam<std::string>("fy_string")),
_fa1_string(getParam<std::string>("fa1_string")),
_fd1_string(getParam<std::string>("fd1_string")),
_fd2_string(getParam<std::string>("fd2_string")),
_K_filettata(declareProperty<RealTensorValue>("conductivityTensor")),
_phi(declareProperty<Real>("Porosity")),
_phiFracture(getParam<Real>("phi_f")),
_phiMatrix(getParam<Real>("phi_m")),
_K_matrix(getParam<Real>("K_matrix")),
_cond0(getParam<bool>("cond0")),
_cond1(getParam<bool>("cond1")),
_gradP(parameters.isParamValid("pressure") ? coupledGradient("pressure"): _grad_zero),
_U(declareProperty<RealVectorValue>("VelocityVector")),
_numOfFrac(declareProperty<Real>("numero")),
_regionID(declareProperty<int>("RegionID"))
{
    _center   =new RealVectorValue [_fn];
    _rotation =new RealVectorValue [_fn];
    _dimension=new RealVectorValue [_fn];
    _d        =new RealVectorValue [_fn];
    
    shall_check= new bool [_fn];
    
    _n = new RealVectorValue * [_fn];
    for (int i=0; i<_fn; ++i)
    {
        _n[i]=new RealVectorValue [2];
    }
    
    _max_f=new RealVectorValue[_fn];
    _min_f=new RealVectorValue[_fn];
    
    std::istringstream  fx_ss( _fx_string);
    std::istringstream  fy_ss( _fy_string);
    std::istringstream fa1_ss(_fa1_string);
    std::istringstream fd1_ss(_fd1_string);
    std::istringstream fd2_ss(_fd2_string);
    
    std::string token;
    
    for (int i=0; i<_fn; ++i)
    {
        std::getline(fx_ss, token, ',');
        _center[i](0)=std::atof(token.c_str());
        std::getline(fy_ss, token, ',');
        _center[i](1)=std::atof(token.c_str());
        
        std::getline(fa1_ss, token, ',');
        _rotation[i](0)=atof(token.c_str());
        _rotation[i]=_rotation[i]/180.0*pi;
        
        std::getline(fd1_ss, token, ',');
        _dimension[i](0)=atof(token.c_str());
        std::getline(fd2_ss, token, ',');
        _dimension[i](1)=atof(token.c_str());

        
        ComputeNormalsFromAngles(_rotation[i],_n[i][0],_n[i][1]);
        
        //        std::cout<<_n[0][0]<<std::endl;
        //        std::cout<<_n[0][1]<<std::endl;
        //        std::cout<<_n[0][2]<<std::endl;
        //        exit(1);
        
        for (int j=0; j<2; ++j)
        {
            _d[i](j)=_n[i][j]*_center[i];
        }
    }
    
    _min_dimension=1000000000;
    for (int i=0; i<_fn; ++i)
        for (int j=0; j<2; ++j)
            _min_dimension=std::min(_min_dimension,_dimension[i](j));
    
    
    _fx_string.clear();
    _fy_string.clear();
    _fa1_string.clear();
    _fd1_string.clear();
    _fd2_string.clear();

    
    _prec = 10000000;
    srand (time(NULL));
    
    
    RealVectorValue * vertices;
    vertices=new RealVectorValue[4];
    
    for (int l=0; l<_fn; ++l)
    {
        
        
        int counter =0;
        for (int i=-1; i<=1; ++++i)
            for (int j=-1; j<=1; ++++j)
                {
                    vertices[counter] = _center[l]+
                                        0.5*_dimension[l](0)*i*_n[l][0]+
                                        0.5*_dimension[l](1)*j*_n[l][1];
                    counter+=1;
                }
        
        _max_f[l]=vertices[0];
        _min_f[l]=vertices[0];
        
        for (int i=1; i<4; ++i)
        {
            for (int j=0; j<2; ++j)
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
HydraulicConductivity2D::computeQpProperties()
{
    RealVectorValue point=_q_point[_qp];
    
    _regionID[_qp]=-1;
    
    std::vector<int> ciao;
    
    _regionID[_qp]=findRegion(point,ciao);
    
    if (ciao.size()>1)
    {
        std::cout<<"maggiore\n";
        exit(1);
    }
    
    Real permFrac;
    Real permMatrix=_K_matrix;
    Real permMatrixReduced=0.1;
    
    if (_cond0)
    {
        permFrac=1.0e4;
    }
    else
    {
        permFrac=1.0e-4;
    }
    
    int count = is_inside(point);
    _numOfFrac[_qp]=count;
    
    if (count>0)
    {
        //We are inside at least a fracture
        _phi[_qp]=_phiFracture;
        _K_filettata[_qp]= permFrac * _identity;
        if (count > 2)
        {
            std::cout<<"error\n";
            exit(1);
        }
    }
    else
    {
        _K_filettata[_qp]= permMatrix * _identity;
        _phi[_qp]=_phiMatrix;

    }
    
    Real a=_JxW[_qp];
    
    _U[_qp] =  -1.0 *  _K_filettata[_qp] * _gradP[_qp];
    
}

void HydraulicConductivity2D::ComputeNormalsFromAngles(RealVectorValue const & angles,
                                                       RealVectorValue & n1,
                                                       RealVectorValue & n2)
{
    RealTensorValue R1;

    
    R1(0,0)=std::cos(angles(0));
    R1(0,1)=-std::sin(angles(0));
    R1(1,0)=std::sin(angles(0));
    R1(1,1)=std::cos(angles(0));

    
    RealTensorValue R=R1;
    
    for (int i=0; i<2; ++i)
    {
        n1(i)=R(i,0);
        n2(i)=R(i,1);
    }
}

int HydraulicConductivity2D::is_inside(RealVectorValue const & point)
{
    int numOfFrac=0;
    _whichFrac.clear();
    for (int i=0; i < _fn; ++i)
    {
        
        Real temp1=std::fabs( _n[i][0]*point - 1.0 * _d[i](0) );
        if (temp1<_dimension[i](0)/2.0)
        {
            Real temp2=std::fabs( _n[i][1]*point - 1.0 * _d[i](1) );
            if (temp2<_dimension[i](1)/2.0)
            {
                numOfFrac+=1;
                _whichFrac.push_back(i);
                
            }
        }
        
    }
    return numOfFrac;
}

void HydraulicConductivity2D::outerProduct
(RealVectorValue const & in1, RealVectorValue const &  in2, RealTensorValue & out)
{
    
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
        {
            out(i, j) = in1(i)*in2(j);
        }
}

int HydraulicConductivity2D::findRegion(RealVectorValue const & point,std::vector<int> & in)
{
    in.clear();
    int returnValue = -1;
    
    for ( int i=0; i<_regionMin.size(); ++i )
    {
        bool isIn=1;
        for (int dim=0; dim<2; ++dim)
        {
            if ( _regionMin[i](dim) < point(dim) && point(dim) < _regionMax[i](dim) )
            {
                // do nothing
            }
            else
            {
                isIn = 0;
            }
        }
        if (isIn)
        {
            if (returnValue!= -1)
            {
                std::cout<<"c'e' un problema\n";
                exit(1);
            }
            returnValue = i;
            in.push_back(i);
        }
    }
    
    return returnValue;
}
