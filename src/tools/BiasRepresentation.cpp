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
#include "core/Value.h"
#include "Communicator.h"
#include "BiasRepresentation.h"
#include <iostream>

namespace PLMD {

//+PLUMEDOC INTERNAL kernelfunctions
/*
*/
//+ENDPLUMEDOC

/// the constructor here
BiasRepresentation::BiasRepresentation(vector<Value*> tmpvalues, Communicator &cc ):hasgrid(false),mycomm(cc){
    ndim=tmpvalues.size();
    for(int i=0;i<ndim;i++){
         values.push_back(tmpvalues[i]);
         names.push_back(values[i]->getName());
         //cerr<<"NN "<<names.back()<<endl;
         //cerr<<"ptr "<<&(values[i])<<endl;
    } 
}; 
/// overload the constructor: add the grid at constructor time 
BiasRepresentation::BiasRepresentation(vector<Value*> tmpvalues, Communicator &cc , vector<string> gmin, vector<string> gmax, vector<unsigned> nbin ):hasgrid(true), mycomm(cc){
    ndim=tmpvalues.size();
    for(int  i=0;i<ndim;i++){
         values.push_back(tmpvalues[i]);
         names.push_back(values[i]->getName());
    } 
    // initialize the grid 
    hasgrid=true;	
    vector<Value*> vv;for(unsigned i=0;i<values.size();i++)vv.push_back(values[i]);
    string ss; ss="file.bias"; 
    cerr<<" initializing grid "<<endl;
    BiasGrid_=new Grid(ss,vv,gmin,gmax,nbin,false,true);
}; 


unsigned BiasRepresentation::getNumberOfDimensions(){
    return values.size();
}; 
vector<string> BiasRepresentation::getNames(){
    return names;
}; 
const string & BiasRepresentation::getName(unsigned i){
    return names[i];
}; 

const vector<Value*>& BiasRepresentation::getPtrToValues(){
    return values;
}; 
Value*  BiasRepresentation::getPtrToValue(unsigned i){
    return values[i];
}; 
void BiasRepresentation::pushGaussian( const vector <double> &sigma, const double &height ){
        vector<double> center;
        for(int i=0;i<ndim ;i++)center.push_back(values[i]->get()); 
 	KernelFunctions *kk=new  KernelFunctions( center ,  sigma, "gaussian",  height, false );
 	hills.push_back(kk);		
        // if grid is defined then it should be added on the grid    
 	cerr<<"now with "<<hills.size()<<endl;
        if(hasgrid){
                 vector<unsigned> nneighb=kk->getSupport(BiasGrid_->getDx());
                 vector<unsigned> neighbors=BiasGrid_->getNeighbors(center,nneighb);
                 vector<double> der(ndim);
                 vector<double> xx(ndim);
                 if(mycomm.Get_size()==1){
                  for(int i=0;i<neighbors.size();++i){
                   unsigned ineigh=neighbors[i];
                   for(int j=0;j<ndim;++j){der[j]=0.0;}
                   BiasGrid_->getPoint(ineigh,xx);   
                   // assign xx to a new vector of values
                   for(int j=0;j<ndim;++j){values[j]->set(xx[j]);}	 
                   double bias=kk->evaluate(values,der,true);
                   BiasGrid_->addValueAndDerivatives(ineigh,bias,der);
                  } 
                 } else {
                  unsigned stride=mycomm.Get_size();
                  unsigned rank=mycomm.Get_rank();
                  vector<double> allder(ndim*neighbors.size(),0.0);
                  vector<double> allbias(neighbors.size(),0.0);
	          vector<double> tmpder(ndim); 
                  for(unsigned i=rank;i<neighbors.size();i+=stride){
                   unsigned ineigh=neighbors[i];
                   BiasGrid_->getPoint(ineigh,xx);
                   for(int j=0;j<ndim;++j){values[j]->set(xx[j]);}	 
                   allbias[i]=kk->evaluate(values,tmpder,true);
 	           // this solution with the temporary vector is rather bad, probably better to take 		
		   // a pointer of double as it was in old gaussian 
                   for(int j=0;j<ndim;++j){ allder[ndim*i+j]=tmpder[j];tmpder[j]=0.;}
                  }
                  mycomm.Sum(&allbias[0],allbias.size());
                  mycomm.Sum(&allder[0],allder.size());
                  for(unsigned i=0;i<neighbors.size();++i){
                   unsigned ineigh=neighbors[i];
                   for(unsigned j=0;j<ndim;++j){der[j]=allder[ndim*i+j];}
                   BiasGrid_->addValueAndDerivatives(ineigh,allbias[i],der);
                  }
                }
        };
};










}
