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

#include "Analysis.h"
#include "PlumedMain.h"
#include "Atoms.h"

namespace PLMD {

void Analysis::registerKeywords( Keywords& keys ){
  ActionPilot::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.use("ARG");
  keys.add("compulsory","STRIDE","the frequency with which data should be stored for analysis");
  keys.reserve("compulsory","RUN","the frequency with which to run the analysis algorithm");
  keys.add("optional","BIAS","the bias acting on the system that we are using to reweight the data");
  keys.add("compulsory","TEMP","the system temperature");
  keys.add("optional","REWEIGHT_TEMP","reweight the data at this temperature");
}

Analysis::Analysis(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
idata(0)
{
  if( keywords.exists("RUN") ){
      unsigned ndata; parse("RUN",freq );
      if( freq%getStride()!= 0 ) error("Frequncy of running is not a multiple of the stride");
      ndata=std::floor(freq/getStride() );
      data.resize( ndata, getNumberOfArguments() ); weights.resize( ndata );
  } else {
      freq=getStride();
      data.resize( 1, getNumberOfArguments() );
      weights.resize(1);
  }

  // Can probably add a faculty for restarts here

  parse("TEMP",simtemp);
  rtemp=0; parse("REWEIGHT_TEMP",rtemp);
  if( rtemp!=0 ) needeng=true; 
  else needeng=false;

  // Need to do stuff to get the bias
  parseArgumentList("BIAS",biases);
  // Check everything is a bias
  std::string thename;
  for(unsigned i=0;i<biases.size();++i){
      thename=biases[i]->getName(); 
      std::size_t dot=thename.find_first_of('.');
      if(thename.substr(dot+1)!="bias") error("value " + thename + " is not a bias"); 
  }
}

void Analysis::prepare(){
  if(needeng) plumed.getAtoms().setCollectEnergy(true);
}

void Analysis::calculate(){
  // Get the arguments and store them in a matrix
  for(unsigned i=0;i<getNumberOfArguments();++i) data(idata,i)=getArgument(i);
  // We reweight according to the temperature and eventually according to the bias
  if(needeng){
     double energy=plumed.getAtoms().getEnergy();
     weights[idata]=exp( -( (1.0/rtemp) - (1.0/simtemp) )*energy / plumed.getAtoms().getKBoltzmann() );
  } else {
     weights[idata]=1.0;
  }
  if( biases.size()>0 ){
      double bias=0.0; 
      for(unsigned i=0;i<biases.size();++i) bias+=biases[i]->get();
      weights[idata]*=exp( bias/( plumed.getAtoms().getKBoltzmann()*simtemp ) );
  }
  idata++;
}

void Analysis::update(){
  if( getStride()%freq==0 ){
      plumed_assert( idata==weights.size() );
      performAnalysis(); idata=0;
  }
}

}
