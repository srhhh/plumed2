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
#ifndef __PLUMED_Analysis_h
#define __PLUMED_Colvar_h

#include "Matrix.h"
#include "ActionPilot.h"
#include "ActionWithArguments.h"

#define PLUMED_ANALYSIS_INIT(ao) Action(ao),Analysis(ao)

namespace PLMD {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new methods for analyzing the trajectory, within it there 
is information as to how to go about implementing a new analysis method.

*/

class Analysis :
  public ActionPilot,
  public ActionWithArguments
  {
private:
/// The frequency with which we are performing analysis
  unsigned freq;
/// The temperature at which we are running the calculation
  double simtemp;
/// The temperature at which we want the histogram
  double rtemp;
/// Do we need the energy (are we reweighting at a different temperature)
  bool needeng;
/// Is there a bias acting on the system
  bool hasbias;
/// The piece of data we are inserting
  unsigned idata;
/// The data we are going to analyze
  Matrix<double> data;
/// The weights of all the data points
  std::vector<double> weights;
protected:
/// Retrieve the ith point
  void getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const ;
public:
  static void registerKeywords( Keywords& keys );
  Analysis(const ActionOptions&);
  void prepare();
  void calculate();
  void update();
  virtual void performAnalysis()=0;
  void apply(){}
};

inline
void Analysis::getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const {
  plumed_assert( idata<weights.size() &&  point.size()==data.ncols() );
  for(unsigned i=0;i<point.size();++i) point[i]=data(idata,i);
  weight=weights[idata];
}

}

#endif
