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
#ifndef __PLUMED_RMSDMap_h
#define __PLUMED_RMSDMap_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include "ActionWithDistribution.h"
#include "RMSD.h"
#include <vector>

namespace PLMD {

class RMSDMap :
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithDistribution
  {
private:
  double lambda;
  std::vector<RMSD> frames; 
  Value thevalue;
  std::vector< std::vector<double> > low_dims;
  std::vector<Vector> derivs;
  Tensor vir;
public:
  RMSDMap(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  void calculate();
/// Apply the forces on the values
  void apply();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Return the number of Colvars this is calculating
  unsigned getNumberOfFunctionsInAction();
/// Return the number of derivatives for a given colvar
  unsigned getNumberOfDerivatives( const unsigned& j );
/// Retrieve the value that was calculated last
  void retreiveLastCalculatedValue( Value& myvalue );
/// Activate the value
  void activateValue( const unsigned j ){}
/// Ensure that nothing gets done for your deactivated colvars
  void deactivateValue( const unsigned j ){}
/// Merge the derivatives
  void mergeDerivatives( const unsigned j, const Value& value_in, const double& df, Value& value_out );
/// Are the base quantities periodic
  bool isPeriodic(); 
/// Calculate one of the functions in the distribution
  bool calculateThisFunction( const unsigned& j );
/// Get the dimensionality that we are mapping into
  unsigned getLowDim() const ;
/// Get a point from the projection
  double getProjectionPoint(const unsigned& jp, const unsigned& kp ) const ;
};

inline
unsigned RMSDMap::getLowDim() const {
  return low_dims[0].size();
}

inline
double RMSDMap::getProjectionPoint(const unsigned& jp, const unsigned& kp ) const {
  plumed_assert( jp<frames.size() && kp<low_dims[jp].size() );
  return low_dims[jp][kp];
}

inline
unsigned RMSDMap::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms()+9;
}

inline
unsigned RMSDMap::getNumberOfFunctionsInAction(){
  return frames.size();
}

inline
unsigned RMSDMap::getNumberOfDerivatives( const unsigned& j ){
  plumed_assert( j<frames.size() );
  return 3*getNumberOfAtoms()+9;
}

inline
void RMSDMap::retreiveLastCalculatedValue( Value& myvalue ){
  copy( thevalue, myvalue );
}

inline
void RMSDMap::mergeDerivatives( const unsigned j, const Value& value_in, const double& df, Value& value_out ){
  copy( value_in, value_out ); value_out.chainRule(df);
}

inline
bool RMSDMap::isPeriodic(){
  return false;
}

}
#endif
