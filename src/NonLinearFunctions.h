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
#ifndef __PLUMED_NonLinearFunction_h
#define __PLUMED_NonLinearFunction_h

#include <string>
#include "Tools.h"

namespace PLMD {

/// This class contains various functions that map the real numbers 
/// to values between 0 and 1. For all these functions f(0)>f(1).
/// The following functions are available:
///
/// rational    : \f$ f(x) = \frac{ 1 - x^n }{ 1 - x^m } \f$
/// triangular  : \f$ f(x) = 1 - |x| \qquad |x| < 1 \qquad f(x)=0 \qquad |x|>1 \f$
/// step        : \f$ f(x) = 1       \qquad |x| < 1 \qquad f(x)=0 \qquad |x|>1 \f$
/// exponential : \f$ f(x) = exp(-x) \f$
/// gaussian    : \f$ f(x) = exp(-x^2/2) \f$
///
/// This class is used within SwitchingFunctions, reversedSwitchingFunctions and 
/// Kernels. In essence the classes that contain instances of the NonLinearFunction
/// class do some form of linear algebra on their input. The result of this linear
/// algebra is a single real number.  Obviously, this real number
/// can be transformed by a variety by a variety of non-linear functions. This class
/// thus provides a number of non-linear functions.   

class NonLinearFunction {
private:
/// Has setup been called
  bool setup;
/// The parameters for rational functions
  unsigned nn, mm;
/// The type of non linear function we are using
  enum {rational,triangular,step,exponential,gaussian} type;
public:
  static std::string printLatexFunctions( const std::string& x, const bool needderiv );
  static std::string printInput( const std::string& aparams, const bool needderiv );
/// The constructor
  NonLinearFunction();
/// This is the function we call to setup the non linear function
  void set( const std::string& type, std::vector<std::string>& data, bool needderiv=true );
/// Does this function take the value of x squared as input
  bool inputXSquared() const ;
/// Get the name of the particular function we are using
  std::string getName() const ;
/// Calculate the value of the function
  double calculate( const double& x, double& dx ) const ;
/// Get the volume of the function
  double getVolume( const unsigned& dim, const double& determinant ) const ;
/// Write out the parameters of the function
  std::string writeParameters() const;
/// Get the value of the cutoff
  double getCutoff( const double& width ) const ;
};

inline
bool NonLinearFunction::inputXSquared() const {
  return (type==gaussian);
}

}
#endif

