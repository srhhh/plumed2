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
#ifndef __PLUMED_vesselbase_FunctionAndDerivativesOnGrid_h
#define __PLUMED_vesselbase_FunctionAndDerivativesOnGrid_h

#include "GridVesselBase.h"

namespace PLMD {
namespace vesselbase{

class FunctionAndDerivativesOnGrid : public GridVesselBase {
private:
/// The derivatives wrt to the low-dimensional properties
  std::vector<double> derlow;
/// Were forces applied on this object
  bool wasforced;
/// The forces that are acting on each of the derivatives in this object
  std::vector<double> forces;
protected:
/// Add value to the field of values, add low-dimensional derivatives and high-dimensional derivatives
  void accumulate( const double& , const double& , const double& , const double& , const unsigned& ); 
public:
/// Create the keywords
  static void registerKeywords( Keywords& keys );
/// The constructor
  FunctionAndDerivativesOnGrid( const VesselOptions& );
/// Resize the field
  void resize();
/// Apply some forces to the field
  bool applyForce(std::vector<double>&);
/// Set the forces on the quantities underlying the fields
  void setForces( const std::vector<double>& );
};

inline
void FunctionAndDerivativesOnGrid::setForces( const std::vector<double>& ff ){
  plumed_dbg_assert( ff.size()==forces.size() );
  wasforced=true;
  for(unsigned i=0;i<ff.size();++i) forces[i]=ff[i];
}

}
}
#endif

