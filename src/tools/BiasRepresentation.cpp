/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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

//+PLUMEDOC INTERNAL BiasRepresentation 
/*
An internal tool in plumed that is used to represent a bias

*/
//+ENDPLUMEDOC

/// the constructor here
BiasRepresentation::BiasRepresentation(vector<Value*> tmpvalues, Communicator &cc ):hasgrid(false),mycomm(cc){
    ndim=tmpvalues.size();
    for(int i=0;i<ndim;i++){
         values.push_back(tmpvalues[i]);
         names.push_back(values[i]->getName());
    } 
}
/// overload the constructor: add the sigma  at constructor time 
BiasRepresentation::BiasRepresentation(vector<Value*> tmpvalues, Communicator &cc,  vector<double> sigma ):hasgrid(false),histosigma(sigma),mycomm(cc){
    ndim=tmpvalues.size();
    for(int i=0;i<ndim;i++){
         values.push_back(tmpvalues[i]);
         names.push_back(values[i]->getName());
    } 
} 
/// overload the constructor: add the grid at constructor time 
BiasRepresentation::BiasRepresentation(vector<Value*> tmpvalues, Communicator &cc , vector<string> gmin, vector<string> gmax, vector<unsigned> nbin ):hasgrid(false), rescaledToBias(false), mycomm(cc){
    ndim=tmpvalues.size();
    for(int  i=0;i<ndim;i++){
         values.push_back(tmpvalues[i]);
         names.push_back(values[i]->getName());
    } 
    // initialize the grid 
    addGrid(gmin,gmax,nbin);
    // test the grids
	//int pp=Grid2::gridTester();
} 
/// overload the constructor with some external sigmas: needed for histogram
BiasRepresentation::BiasRepresentation(vector<Value*> tmpvalues, Communicator &cc , vector<string> gmin, vector<string> gmax, vector<unsigned> nbin , vector<double> sigma):hasgrid(false), rescaledToBias(false),histosigma(sigma),mycomm(cc){
    ndim=tmpvalues.size();
    for(int  i=0;i<ndim;i++){
         values.push_back(tmpvalues[i]);
         names.push_back(values[i]->getName());
    }
    // initialize the grid 
    addGrid(gmin,gmax,nbin);
}


void  BiasRepresentation::addGrid( vector<string> gmin, vector<string> gmax, vector<unsigned> nbin ){
    plumed_massert(hills.size()==0,"you can set the grid before loading the hills");
    plumed_massert(hasgrid==false,"to build the grid you should not having the grid in this bias representation");
    string ss; ss="file.free"; 
    vector<Value*> vv;for(unsigned i=0;i<values.size();i++)vv.push_back(values[i]);
    //cerr<<" initializing grid "<<endl;
    //BiasGrid_=new Grid(ss,vv,gmin,gmax,nbin,false,true);
    BiasGrid2_=new Grid2(ss,vv,gmin,gmax,nbin,false,true);

    hasgrid=true;	

}
bool BiasRepresentation::hasSigmaInInput(){
   if(histosigma.size()==0){return false;}else{return true;} 
}
void BiasRepresentation::setRescaledToBias(bool rescaled){
	plumed_massert(hills.size()==0,"you can set the rescaling function only before loading hills");
        rescaledToBias=rescaled;
}
const bool & BiasRepresentation::isRescaledToBias(){
	return rescaledToBias;
}

unsigned BiasRepresentation::getNumberOfDimensions(){
    return values.size();
} 
vector<string> BiasRepresentation::getNames(){
    return names;
} 
const string & BiasRepresentation::getName(unsigned i){
    return names[i];
} 

const vector<Value*>& BiasRepresentation::getPtrToValues(){
    return values;
} 
Value*  BiasRepresentation::getPtrToValue(unsigned i){
    return values[i];
} 

KernelFunctions* BiasRepresentation::readFromPoint(IFile *ifile){
	vector<double> cc( names.size() );
	for(unsigned i=0;i<names.size();++i){
         ifile->scanField(names[i],cc[i]);
        }
        double h=1.0; 
	return new KernelFunctions(cc,histosigma,"gaussian",false,h,false);	
}
void BiasRepresentation::pushKernel( IFile *ifile ){
	KernelFunctions *kk;
	// here below the reading of the kernel is completely hidden
	if(histosigma.size()==0){
		ifile->allowIgnoredFields();
		kk=KernelFunctions::read(ifile,names)   ;
	}else{
		// when doing histogram assume gaussian with a given diagonal sigma 
		// and neglect all the rest
		kk=readFromPoint(ifile)   ;
	}
	hills.push_back(kk);
	// the bias factor is not something about the kernels but 
	// must be stored to keep the  bias/free energy duality
	string dummy; double dummyd;
	if(ifile->FieldExist("biasf")){
		ifile->scanField("biasf",dummy);
		Tools::convert(dummy,dummyd);
	}else{dummyd=1.0;}
	biasf.push_back(dummyd);
	// the domain does not pertain to the kernel but to the values here defined
	string	mins,maxs,minv,maxv,mini,maxi;mins="min_";maxs="max_";
	for(unsigned i=0 ; i<unsigned(ndim); i++){
		if(values[i]->isPeriodic()){
			ifile->scanField(mins+names[i],minv);
			ifile->scanField(maxs+names[i],maxv);
			// verify that the domain is correct
			values[i]->getDomain(mini,maxi);
			plumed_massert(mini==minv,"the input periodicity in hills and in value definition does not match"  );
			plumed_massert(maxi==maxv,"the input periodicity in hills and in value definition does not match"  );
		}
	}
	// if grid is defined then it should be added on the grid
	if(hasgrid){

		//vector<unsigned> nneighb=kk->getSupport(BiasGrid_->getDx());

		vector<unsigned> nneighb2=kk->getSupport(BiasGrid2_->getDx());


		//vector<unsigned> neighbors=BiasGrid_->getNeighbors(kk->getCenter(),nneighb);

		vector<unsigned> neighbors2=BiasGrid2_->getNeighbors(kk->getCenter(),nneighb2);

		vector<double> der(ndim);
		vector<double> xx(ndim);
		if(mycomm.Get_size()==1){
			for(unsigned i=0;i<neighbors2.size();++i){
				unsigned ineigh=neighbors2[i];
				for(unsigned j=0;j<unsigned(ndim);++j){der[j]=0.0;}

				//cerr<<endl;
				BiasGrid2_->getPoint(ineigh,xx);


				// assign xx to a new vector of values
				for(int j=0;j<ndim;++j){values[j]->set(xx[j]);}
				double bias=kk->evaluate(values,der,true);
				if(rescaledToBias){
					double f=(biasf.back()-1.)/(biasf.back());
					bias*=f;
					for(int j=0;j<ndim;++j){der[j]*=f;}
				}

				BiasGrid2_->addValueAndDerivatives(ineigh,bias,der);


			}

		} else {
			unsigned stride=mycomm.Get_size();
			unsigned rank=mycomm.Get_rank();
			vector<double> allder(ndim*neighbors2.size(),0.0);
			vector<double> allbias(neighbors2.size(),0.0);
			vector<double> tmpder(ndim);
			for(unsigned i=rank;i<neighbors2.size();i+=stride){
				unsigned ineigh=neighbors2[i];
				BiasGrid2_->getPoint(ineigh,xx);
				for(int j=0;j<ndim;++j){values[j]->set(xx[j]);}
				allbias[i]=kk->evaluate(values,tmpder,true);
				if(rescaledToBias){
					double f=(biasf.back()-1.)/(biasf.back());
					allbias[i]*=f;
					for(int j=0;j<ndim;++j){tmpder[j]*=f;}
				}
				// this solution with the temporary vector is rather bad, probably better to take
				// a pointer of double as it was in old gaussian
				for(int j=0;j<ndim;++j){ allder[ndim*i+j]=tmpder[j];tmpder[j]=0.;}
			}
			mycomm.Sum(&allbias[0],allbias.size());
			mycomm.Sum(&allder[0],allder.size());
			for(unsigned i=0;i<neighbors2.size();++i){
				unsigned ineigh=neighbors2[i];
				for(unsigned j=0;j<unsigned(ndim);++j){der[j]=allder[ndim*i+j];}
				BiasGrid2_->addValueAndDerivatives(ineigh,allbias[i],der);
			}
		}
	}
}

	int BiasRepresentation::getNumberOfKernels(){
	return hills.size();
}
Grid2* BiasRepresentation::getGridPtr(){
        plumed_massert(hasgrid,"if you want the grid pointer then you should have defined a grid before"); 
	//return BiasGrid_;
	return BiasGrid2_;

}
void BiasRepresentation::getMinMaxBin(vector<double> &vmin, vector<double> &vmax, vector<unsigned> &vbin){
	vector<double> ss,cc,binsize; 
	vmin.clear();vmin.resize(ndim,10.e20);
	vmax.clear();vmax.resize(ndim,-10.e20);
	vbin.clear();vbin.resize(ndim);
	binsize.clear();binsize.resize(ndim,10.e20);
	int ndiv=10; // adjustable parameter: division per support
	for(unsigned i=0;i<hills.size();i++){
		if(histosigma.size()!=0){
			ss=histosigma;
		}else{
			ss=hills[i]->getContinuousSupport();	
		}
		cc=hills[i]->getCenter();
		for(unsigned j=0;j<unsigned(ndim);j++){
			double dmin=cc[j]-ss[j];
			double dmax=cc[j]+ss[j];
			double ddiv=ss[j]/double(ndiv);
			if(dmin<vmin[j])vmin[j]=dmin; 
			if(dmax>vmax[j])vmax[j]=dmax; 
			if(ddiv<binsize[j])binsize[j]=ddiv;		
		}
	}
	for(unsigned j=0;j<unsigned(ndim);j++){
		// reset to periodicity
		if(values[j]->isPeriodic()){
//			double minv,maxv;
			// old pbc treatment does not allow cross pbc
//			values[j]->getDomain(minv,maxv);
//			if(minv>vmin[j])vmin[j]=minv;
//			if(maxv<vmax[j])vmax[j]=maxv;
			// allow crossing of periodicity:
			// first bring the periodicity in with three steps
			// 1-rescale from -0.5 to 0.5
			double pmin,pmax;
			values[j]->getDomain(pmin,pmax);
			double v1min=(vmin[j]-pmin)/(pmax-pmin)-0.5;
			double v1max=(vmax[j]-pmin)/(pmax-pmin)-0.5;
			// register if the boundary was crossed from one way or the other
			int deltamin=int(v1min+copysign(0.5,v1min));
			int deltamax=int(v1max+copysign(0.5,v1max));
			bool crossedperiod=false; if(deltamin!=0 || deltamax!=0)crossedperiod=true;
			// 2-bring back into pbc
			v1min=v1min-deltamin;
			v1max=v1max-deltamax;
			// 3-rescale back to the old periodicity
			vmin[j]=(v1min+0.5)*(pmax-pmin)+pmin;
			vmax[j]=(v1max+0.5)*(pmax-pmin)+pmin;
			//
			// if you crossed periodicity in one direction
			// you have 2 cases:  min>max ok
			//                    min<max : take the whole period
			if(crossedperiod && vmin[j]<vmax[j] ){
				vmin[j]=pmin;vmax[j]=pmax;
				vbin[j]=static_cast<unsigned>(ceil((vmax[j]-vmin[j])/binsize[j]) );
			}else{
				vbin[j]=static_cast<unsigned>(ceil((vmax[j]-pmin)/binsize[j])+ceil((pmax-vmin[j])/binsize[j]) );
			}
		}else{
			vbin[j]=static_cast<unsigned>(ceil((vmax[j]-vmin[j])/binsize[j]) );
		}
	}
}
void BiasRepresentation::clear(){
        // clear the hills
	for(vector<KernelFunctions*>::const_iterator it = hills.begin(); it != hills.end(); ++it)
        {
	    delete *it;
        } 
        hills.clear(); 
        // clear the grid
        if(hasgrid){
              BiasGrid2_->clear();
        }
}


}
