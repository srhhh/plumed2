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

#include "NonLinearFunctions.h"
#include "Tools.h"
#include "Matrix.h"
#include "Value.h"
#include <vector>
#include <math.h>

namespace PLMD {

/// \ingroup TOOLBOX
/// Small class to compute kernel functions.  Within plumed a 
/// kernel function has the following form:
/// \f[
///    f( \mathbf{x-x_0}^T \mathbf{\Sigma}^-1 \mathbf{x-x_0} )
/// \f]
/// Within this equation \f$\mathbf{\Sigma}\f$ is a square matrix that
/// describes the covariance in the region of space around \f$\mathbf{x_0}\f$.
/// This matrix can be diagonal. Within plumed you can use various different 
/// function for \f$f()\f$ for more details see \ref PLUMED::NonLinearFunction. 

class Kernel {
private:
/// Is the metric matrix diagonal
  bool diagonal;
/// The center of the kernel function
  std::vector<double> center;
/// The width of the kernel
  std::vector<double> width;
/// The height of the kernel
  double height;
/// The non-linear function we are using
  NonLinearFunction nlfunc;
/// Convert the width into matrix form
  Matrix<double> getMatrix() const;
protected:
/// Get the determinant of the metric
  double getDeterminant() const;
public:
/// Does the kernel have derivatives
  bool hasderivatives;
  Kernel( const std::vector<double>& at, const std::vector<double>& sig, const std::string& type, const double& w, const bool& norm );
/// Get the dimensionality of the kernel
  unsigned ndim() const;
/// Get the position of the center 
  std::vector<double> getCenter() const;
/// Get the support
  std::vector<unsigned> getSupport( const std::vector<double>& dx ) const; 
/// Evaluate the kernel function  
  double evaluate( const std::vector<Value*>& pos, std::vector<double>& derivatives, bool usederiv=true ) const;
/// Print the kernel function to a file
  void print( const std::vector<std::string>& cv_names, PlumedOFile& ofile ) const ;
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

}
#endif
