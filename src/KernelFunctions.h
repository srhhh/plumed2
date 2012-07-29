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
#include <vector>
#include <math.h>

namespace PLMD {

class KernelOptions {
friend class Kernel;
private:
  std::vector<double> pos;
public:
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
/// The center of the kernel function
  std::vector<double> center;
/// The width of the kernel
  std::vector<double> width;
/// The height at the center of the kernel
  double height;
protected:
/// Set the value of the height
  void setHeight( const double& h );
/// Get the value of the height
  double getHeight() const;
public:
/// Does the kernel have derivatives
  bool hasderivatives;
  Kernel( const KernelOptions& ko ); 
/// Return the dimensionality of the kernel
  unsigned ndim() const ;
/// Get the position of the center of the kernel function
  std::vector<double> getCenter() const;
/// Get the numbers of neighbors of the center required in each direction
  std::vector<unsigned> getSupport( const std::vector<double>& dx );
/// Get how far out we need to go from the center
  virtual double getCutoff( double& width )=0;
/// Evaluate the kernel function (this bit looks after PBCs) - it would need to be overwritten to do non-digaonal Normal distributions
  virtual double evaluate( const std::vector<bool>& pbc, const std::vector<double>& range, const std::vector<double>& pos );
/// Get the value of the kernel at this point (note we pass here (p-c)/b)
  virtual double getValue( const std::vector<double>& diff )=0;
};

inline
unsigned Kernel::ndim() const {
  return width.size();
}

inline
void Kernel::setHeight(const double& h){
  height=h;
}

inline
double Kernel::getHeight() const {
  return height;
}

inline
std::vector<double> Kernel::getCenter() const {
  return center;
}

class UniformKernel : public Kernel {
public:
  UniformKernel( const KernelOptions& ko ); 
  double getCutoff( double& width );
  double getValue( const std::vector<double>& diff );
};

inline
double UniformKernel::getCutoff( double& width ){
  return width;
}

inline
double UniformKernel::getValue( const std::vector<double>& diff ){
  for(unsigned i=0;i<diff.size();++i){
      if( diff[i]>1.0 ) return 0.;
  }
  return getHeight();
}

}
#endif
