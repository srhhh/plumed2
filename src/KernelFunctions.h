/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.
  
   This file is part of plumed, version 2.0.
  
   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
  
   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_KernelFunctions_h
#define __PLUMED_KernelFunctions_h

#include "PlumedException.h"
#include "Tools.h"
#include "Matrix.h"
#include "Value.h"
#include <vector>
#include <math.h>

namespace PLMD {

class KernelOptions {
friend class Kernel;
public:
  std::vector<double> pos;
  std::vector<double> width;
  double height;
  bool normalize;
  KernelOptions( const std::vector<double>& , const std::vector<double>& , const double& , const bool& );
};

class Kernel;

class KernelRegister {
public:
  static Kernel* create( const std::string& type, const KernelOptions& ko, const bool& );
};

class Kernel {
private:
/// Is the metric matrix diagonal
  bool diagonal;
/// The center of the kernel function
  std::vector<double> center;
/// The width of the kernel
  std::vector<double> width;
/// Convert the width into matrix form
  Matrix<double> getMatrix() const;
protected:
/// Get the determinant of the metric
  double getDeterminant() const;
public:
/// Does the kernel have derivatives
  bool hasderivatives;
/// You can specify kernels as a function of r2 and thus avoid square roots
  bool is_function_of_r2;
  Kernel( const KernelOptions& ko ); 
/// Return the dimensionality of the kernel
  unsigned ndim() const ;
/// Get the position of the center of the kernel function
  std::vector<double> getCenter() const;
/// Get the numbers of neighbors of the center required in each direction
  std::vector<unsigned> getSupport( const std::vector<double>& dx );
/// Get how far out we need to go from the center
  virtual double getCutoff( double& width )=0;
/// Evaluate the kernel function  
  double evaluate( const std::vector<Value*>& pos, std::vector<double>& derivatives, bool usederiv=true );
/// Get the value of the kernel at this point (note we pass here (p-c)/b)
  virtual double getValue( const double& x, double& dx )=0;
/// Print the header for the kernel function to a file
  std::string fieldNames( const std::vector<std::string>& arg_names );
/// Print the header for the parameters of the kernel
  virtual std::string parameterNames()=0; 
/// Print the kernel function to a file
  void print( FILE* ofile );
/// Print extra parameters of the kernel
  virtual void printParameters( FILE* ofile )=0;
};

inline
Matrix<double> Kernel::getMatrix() const {
  unsigned k=0, ncv=ndim(); Matrix<double> mymatrix(ncv,ncv); 
  for(unsigned i=0;i<ncv;i++){
    for(unsigned j=i;j<ncv;j++){
        mymatrix(i,j)=mymatrix(j,i)=width[k]; // recompose the full inverse matrix
        k++;
    }
  }
  return mymatrix;
}

inline
unsigned Kernel::ndim() const {
  return center.size();
}

inline
std::vector<double> Kernel::getCenter() const {
  return center;
}

class UniformKernel : public Kernel {
private:
  double height;
public:
  UniformKernel( const KernelOptions& ko ); 
  double getCutoff( double& width );
  double getValue( const double& x, double& dx );
  std::string parameterNames();
  void printParameters( FILE* ofile );
};

inline
double UniformKernel::getCutoff( double& width ){
  return width;
}

inline
double UniformKernel::getValue( const double& x, double& dx ){
  if( x>1.0 ) return 0.;
  dx=0;
  return height;
}

class GaussianKernel : public Kernel {
private:
  double height;
  double DP2CUTOFF;
public:
  GaussianKernel( const KernelOptions& ko );
  double getCutoff( double& width );
  double getValue( const double& x, double& dx );
  std::string parameterNames();
  void printParameters( FILE* ofile );
};

inline
double GaussianKernel::getCutoff( double& width ){
  return sqrt(2.0*DP2CUTOFF)*width;
}

inline
double GaussianKernel::getValue( const double& x, double& dx ){
  if( x<DP2CUTOFF ){
      double val=height*exp(-0.5*x); // N.B. x here is x^2 so we can avoid expensive square roots.
      dx=-0.5*val;
      return val; 
  }
  dx=0;
  return 0.0;
}

class TriangularKernel : public Kernel {
private:
  double height;
public:
  TriangularKernel( const KernelOptions& ko );
  double getCutoff( double& width );
  double getValue( const double& x, double& dx );
  std::string parameterNames();
  void printParameters( FILE* ofile );
};

inline
double TriangularKernel::getCutoff( double& width ){
  return width;
}

inline
double TriangularKernel::getValue( const double& x, double& dx ){
  if( x<1.0 ){
     if(x==0) dx=0;
     else if(x>0) dx=-height; 
     else dx=height;
     return height*( 1. - fabs(x) ); 
  }
  dx=0.0;
  return 0;
}

}
#endif
