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
#ifndef __PLUMED_SwitchingFunction_h
#define __PLUMED_SwitchingFunction_h

#include <string>
#include "NonLinearFunctions.h"

namespace PLMD {

class Log;

/// \ingroup TOOLBOX
/// Small class to compure switching functions. Within plumed a
/// switching function has the following form:
/// \f{eqnarray*}{
///   x' &=& \frac{ x - d_0 }{ r_0 } \\
///  \sigma(x') &=& 1 \\qquad x'<0 \\
///  \sigma(x') &=& s(x') \qquad x'>0 \land \qquad x< d_{\textrm{max}} \\
///  \sigma(x') &=& 0 \\qquad x> d_{\textrm{max}}
/// \f} 
/// Within plumed you can use various different functions for \f$s(x')\f$
/// for more details as to your options see \ref PLUMED::NonLinearFunction.
class SwitchingFunction{
  bool init;
  double invr0,d0,dmax;
  NonLinearFunction nlfunc;
public:
  static std::string documentation();
  SwitchingFunction();
  void set(const std::string& definition);
  void set(int nn,int mm,double r0,double d0); 
  std::string description() const ;
  double calculate(double x,double&df)const;
  double get_r0() const;
};

}

#endif

