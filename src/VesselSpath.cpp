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
#include "FunctionVessel.h"
#include "ActionWithArguments.h"
#include "ActionAtomistic.h"
#include "RMSDMap.h"

namespace PLMD {

class spath : public NormedSumVessel {
private:
  unsigned ncoords;
  RMSDMap* myrmsdmap;
public:
  static void reserveKeyword( Keywords& keys );
  spath( const VesselOptions& da );
  void getWeight( const unsigned& i, Value& weight );
  void compute( const unsigned& i, const unsigned& j, Value& theval );
};

PLUMED_REGISTER_VESSEL(spath,"SPATH")

void spath::reserveKeyword( Keywords& keys ){
  keys.reserveFlag("SPATH",false,"calculated the position on the manifold using as a weighted average. "
                                 "The position on the manifold is calculated as "
                                 "\\f$s(r)=\\frac{\\sum_{i=1}^N \\mathbf{v}_i \\exp(-(\\mathbf{R}-\\mathbf{R}_i)/\\lambda) }{\\sum_{i=1}^N\\exp(-(\\mathbf{R}-\\mathbf{R}_i)/\\lambda) } \\f$." );
}

spath::spath( const VesselOptions& da ) :
NormedSumVessel(da)
{
  useNorm();

  ActionWithArguments* aargs=dynamic_cast<ActionWithArguments*>( getAction() );
  if(aargs){
     myrmsdmap=NULL;
     plumed_massert(0,"not implemented yet");
  } else {
     ActionAtomistic* aatoms=dynamic_cast<ActionAtomistic*>( getAction() );
     plumed_massert(aatoms, "spath is used to calculate weighted averages for mappings what are you trying to do");
     //mycvmap=NULL;
     myrmsdmap=dynamic_cast<RMSDMap*>( getAction() );
     plumed_massert(myrmsdmap, "spath is used to calculate weighted averages for mappings what are you trying to do");
     ncoords=myrmsdmap->getLowDim();
  }
  if(ncoords==1){
     addOutput("s");
     log.printf("  value %s.s contains the position along the spath\n",(getAction()->getLabel()).c_str());
  } else {
     for(unsigned i=0;i<ncoords;++i){
       std::string num; Tools::convert(i+1,num); addOutput("s"+num);
       log.printf("  value %s.s%d contains the %d th component of the position on the manifold\n",(getAction()->getLabel()).c_str(),i+1,i+1); 
     }
  }
}

void spath::compute( const unsigned& i, const unsigned& j, Value& theval ){
  plumed_assert( j<ncoords );

  double lowd;
  if( myrmsdmap ){
     myrmsdmap->retreiveLastCalculatedValue( theval );
     lowd=myrmsdmap->getProjectionPoint( i, j );
  } else { 
     getAction()->retreiveLastCalculatedValue( theval );
  }
  theval.chainRule( lowd ); theval.set( lowd*theval.get() );
}

void spath::getWeight( const unsigned& i, Value& weight ){
  getAction()->retreiveLastCalculatedValue( weight );
}

}
