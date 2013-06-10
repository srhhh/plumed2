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
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cfloat>

// these two below can be commented if the code for the tester is commented
#include <time.h>
#include <Grid.h>


#include "Grid2.h"
#include "Tools.h"
#include "core/Value.h"
#include "File.h"
#include "Exception.h"
#include "KernelFunctions.h"

using namespace std;
namespace PLMD{

Grid2::Grid2(const std::string& funcl, std::vector<Value*> args, const vector<std::string> & gmin,
		const vector<std::string> & gmax, const vector<unsigned> & nbin, bool dospline, bool usederiv, bool doclear, bool docommensurate){
	// various checks
	plumed_massert(args.size()==gmin.size(),"grid dimensions in input do not match number of arguments");
	plumed_massert(args.size()==nbin.size(),"grid dimensions in input do not match number of arguments");
	plumed_massert(args.size()==gmax.size(),"grid dimensions in input do not match number of arguments");
	unsigned dim=gmax.size();
	std::vector<std::string> names;
	std::vector<bool> isperiodic;
	std::vector<string> pmin,pmax;
	names.resize( dim );
	isperiodic.resize( dim );
	pmin.resize( dim );
	pmax.resize( dim );
	for(unsigned int i=0;i<dim;++i){
		names[i]=args[i]->getName();
		if( args[i]->isPeriodic() ){
			isperiodic[i]=true;
			args[i]->getDomain( pmin[i], pmax[i] );
		} else {
			isperiodic[i]=false;
			pmin[i]="0.";
			pmax[i]="0.";
		}
	}
	// this is a value-independent initializator
	Init(funcl,names,gmin,gmax,nbin,dospline,usederiv,doclear,isperiodic,pmin,pmax,docommensurate);
}

Grid2::Grid2(const std::string& funcl, const std::vector<string> &names, const std::vector<std::string> & gmin,
		const vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline, bool usederiv, bool doclear, const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin, const std::vector<std::string> &pmax, bool docommensurate ){
	// this calls the initializator
	Init(funcl,names,gmin,gmax,nbin,dospline,usederiv,doclear,isperiodic,pmin,pmax,docommensurate);
}

void Grid2::Init(const std::string& funcl, const std::vector<std::string> &names, const vector<std::string> & gmin,
		const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline, bool usederiv, bool doclear,
		const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin, const std::vector<std::string> &pmax , bool docommensurate){
	fmt_="%14.9f";

	// various checks
	plumed_massert(names.size()==gmin.size(),"grid dimensions in input do not match number of arguments");
	plumed_massert(names.size()==nbin.size(),"grid dimensions in input do not match number of arguments");
	plumed_massert(names.size()==gmax.size(),"grid dimensions in input do not match number of arguments");
	dimension_=gmax.size();
	str_min_=gmin; str_max_=gmax;
	argnames.resize( dimension_ );
	min_.resize( dimension_ );
	max_.resize( dimension_ );
	pbc_.resize( dimension_ );
	pmin_.resize( dimension_ );
	pmax_.resize( dimension_ );
	str_pmin_.resize( dimension_ );
	str_pmax_.resize( dimension_ );
	pdelta_.resize( dimension_ );
	pshift_.resize( dimension_ );
	pbin_.resize( dimension_ ); // now many bins are in one periodicity

	for(unsigned int i=0;i<dimension_;++i){

		argnames[i]=names[i];
		if( isperiodic[i] ){
			pbc_[i]=true;
			str_pmin_[i]=pmin[i];
			str_pmax_[i]=pmax[i];
		} else {
			pbc_[i]=false;
			str_pmin_[i]="";
			str_pmax_[i]="";
		}
		Tools::convert(str_min_[i],min_[i]);
		Tools::convert(str_max_[i],max_[i]);
		Tools::convert(str_pmin_[i],pmin_[i]);
		Tools::convert(str_pmax_[i],pmax_[i]);
		nbin_=nbin;
		// preliminary check: the input min max should not exceed the periodicity
		if(isperiodic[i] ){
			plumed_massert(min_[i]<=pmax_[i] && min_[i]>=pmin_[i],"minimum for the grid should be within periodicity of the variable");
			plumed_massert(max_[i]<=pmax_[i] && max_[i]>=pmin_[i],"maximum for the grid should be within periodicity of the variable");
		}
		//Important: adjust the max of the grid so that the bin is commensurate with the periodicity
		if( isperiodic[i] && docommensurate){
			double deltap=pmax_[i]-pmin_[i];
			double delta;
			if(max_[i]<min_[i]){
				delta=deltap-(min_[i]-max_[i]);
			}else{delta=max_[i]-min_[i];}
			//	  cerr<<"VAR "<<i<<" MIN "<<min_[i]<<" MAX "<<max_[i]<<" PMIN "<<pmin_[i]<<" PMAX "<<pmax_[i]<<endl;
			//	  cerr<<" delta "<<delta<<" deltap "<<deltap<<" F1 "<<int(deltap*nbin_[i]/delta+1.e-9)<<" NBIN "<< nbin_[i]
			//	                                                                                                   <<endl;
			double shift=deltap*nbin_[i]/double(int(deltap*nbin_[i]/delta+1.e-9)) -delta;

			if(shift>1.e-9){
				cerr<<"Adjusting periodicity for variable "<<std::setw(3)<<i<<": ";
				// shift in one direction if one of them is pi
				double maxfact=1.;
				double minfact=1.;
				if(max_[i]== pi){minfact=2.;maxfact=0.;};
				if(min_[i]==-pi){minfact=0.;maxfact=2.;};
				// verify that there is no periodic violation
				cerr<<scientific;
				cerr<<" shift is "<<shift;
				cerr<<fixed;
				if(maxfact>0.){
					cerr<<" max shifted from "<<setprecision(9)<<max_[i];
					max_[i]+=shift/2;
					cerr<<" to "<<max_[i]<<" and ";
				}
				if(minfact>0.){
					cerr<<" min shifted from "<<setprecision(9)<<min_[i];
					min_[i]-=shift/2;
					cerr<<" to "<<min_[i]<<endl;
				}

				//		  // safe shifting: in case you shifted too much
				//		  if(max_[i]>min_[i]){
				//			  if(max_[i]>pmax_[i])max_[i]=pmax_[i];
				//			  if(min_[i]<pmin_[i])min_[i]=pmin_[i];
				//		  }else{
				//			  if(max_[i]<pmin_[i])max_[i]=pmin_[i];
				//			  if(min_[i]>pmax_[i])min_[i]=pmax_[i];
				//		  }
				plumed_massert(max_[i]<=pmax_[i] && min_[i]>=pmin_[i] && max_[i]>=pmin_[i] && min_[i]<=pmax_[i] ,"Correction for boundaries is not working");

				ostringstream mm,nn;
				mm<<max_[i];
				str_max_[i]=mm.str();
				nn<<min_[i];
				str_min_[i]=nn.str();

			}
		}

		funcname=funcl;
		if(!isperiodic[i]){// accept the case of min>max when doing across periodic boundary integration
			plumed_massert(max_[i]>min_[i],"maximum in grid must be larger than minimum");
		}else{
			plumed_massert(min_[i]>=pmin_[i] && min_[i]<=pmax_[i],"grid minimum must be enclosed in the periodicity");
			plumed_massert(max_[i]>=pmin_[i] && max_[i]<=pmax_[i],"grid maximum must be enclosed in the periodicity");
			pdelta_[i]=pmax_[i]-pmin_[i];
		}

		plumed_massert(nbin[i]>0,"number of grid points must be greater than zero");
	}

	dospline_=dospline;
	usederiv_=usederiv;
	if(dospline_) plumed_assert(dospline_==usederiv_);
	maxsize_=1;
	for(unsigned int i=0;i<dimension_;++i){
		// find the center of the periodicity, even when the boundaries are reversed
		if(pbc_[i]){
			if(max_[i]<min_[i]){
				pshift_[i]=pdelta_[i];
			}else{
				pshift_[i]=0.;
			}
			dx_.push_back( (max_[i]+pshift_[i]-min_[i])/static_cast<double>( nbin_[i] ) );
			pbin_[i]=floor((pmax_[i]-pmin_[i])/dx_[i]+0.5);
		}else{
			dx_.push_back( (max_[i]-min_[i])/static_cast<double>( nbin_[i] ) );
			pbin_[i]=0;
		}

		if( !pbc_[i] ){ max_[i] += dx_[i]; nbin_[i] += 1; }
		maxsize_*=nbin_[i];
	}
	if(doclear) clear();
}

void Grid2::clear(){
	grid_.resize(maxsize_);
	if(usederiv_) der_.resize(maxsize_);
	for(unsigned int i=0;i<maxsize_;++i){
		grid_[i]=0.0;
		if(usederiv_){
			(der_[i]).resize(dimension_);
			for(unsigned int j=0;j<dimension_;++j) der_[i][j]=0.0;
		}
	}
}

vector<std::string> Grid2::getMin() const {
	return str_min_;
}

vector<std::string> Grid2::getMax() const {
	return str_max_;
}

vector<std::string> Grid2::getPeriodicMin() const {
	return str_pmin_;
}

vector<std::string> Grid2::getPeriodicMax() const {
	return str_pmax_;
}

vector<double> Grid2::getDx() const {
	//	 for(unsigned int i=0;i<dimension_;++i){
	//		 cerr<<" DX "<<dx_[i]<<endl;
	//	 }
	return dx_;
}

double Grid2::getBinVolume() const {
	double vol=1.;
	for(unsigned i=0;i<dx_.size();++i) vol*=dx_[i];
	return vol;
}

vector<bool> Grid2::getIsPeriodic() const {
	return pbc_;
}

vector<unsigned> Grid2::getNbin() const {
	return nbin_;
}

vector<string> Grid2::getArgNames() const {
	return argnames;
}


unsigned Grid2::getSize() const {
	return maxsize_;
}

unsigned Grid2::getDimension() const {
	return dimension_;
}

// if you want to apply pbc in a grid is a bit tricky since
// the grid can be not commensurate to the bin: you need to convert all the time
// into coordinates, then apply pbc, the apply back the index
void Grid2::applyPeriodicityToIndices(vector<int> &indices) const {
	plumed_massert(indices.size()==dimension_,"indices in input must have the same dimension of the grid");
	for(unsigned i=0;i<dimension_;i++){
		if(pbc_[i] && (indices[i]<0 || indices[i]>=int(nbin_[i])) ){
			// convert to coordinate
			//cerr<<"INDINIT "<<indices[i]<<"\n";
			double x=indices[i]*dx_[i]+min_[i];
			//cerr<<"XINIT "<<setprecision(12)<<x<<"\n";
			// bring into -0.5:0.5
			x=(x-pmin_[i])/pdelta_[i]-0.5;
			// back into period
			//cerr<<"XRED "<<setprecision(20)<<x<<endl;
			x=x-int(x+copysign(0.5,x));
			//cerr<<"XREDAFTER "<<setprecision(20)<<x<<endl;
			// back to standard coor
			x=(x+0.5)*pdelta_[i]+pmin_[i];
			//cerr<<"XAFTER "<<setprecision(20)<<x<<"\n";
			//get back to index: use safe rounding
			indices[i]=floor((x-min_[i])/dx_[i]+0.5);
			//cerr<<"INDAFTER "<<indices[i]<<" ratio "<<setprecision(20)<<(x-min_[i])/dx_[i]<<"\n";
		}
	}
}

//// we are flattening arrays using a column-major order
bool Grid2::getIndex(const vector<int> &  indices, unsigned &index) const {
	vector<int> myindices=indices;
	// bring indices to the periodicity to have rapid conversion with positions
	//	 for(unsigned i=0;i<dimension_;i++){
	//		 cerr<<"INDEX BEFORE PBC SHIFT "<<myindices[i]<<" BIN "<<nbin_[i]<<endl;
	//
	//	 }

	// cerr<<"enter pbc calc"<<endl;
	applyPeriodicityToIndices(myindices);
	//cerr<<"enter after calc"<<endl;


	//	 for(unsigned i=0;i<dimension_;i++){
	//		 cerr<<"INDEX AFTER PBC SHIFT "<<myindices[i]<<" BIN "<<nbin_[i]<<endl;
	//	 }

	for(unsigned i=0;i<dimension_;i++){
		// if pbc the index is positive or negative
		// but all the time within the periodicity
		// if negative, first shift up to the next period,
		// then calculate the index
		if(pbc_[i]){
			if(max_[i]>min_[i]){// normal pbc. Minimum is the offset
				if(myindices[i]<0 || myindices[i]>=int(nbin_[i]))return false;
			}else{
				if(myindices[i]<0){
					// shift up of one periodicity
					myindices[i]+=pbin_[i];
					//		 cerr<<"INDEX AFTER_NEXT_PERIOD SHIFT "<<myindices[i]<<" PERBIN "<<(pbin_[i])<<endl;
					if(myindices[i]>int(nbin_[i])) return false;
				}
			}
		}else{
			if(myindices[i]<0 || myindices[i]>=int(nbin_[i]))return false;
		}

	}

	index=unsigned(myindices[dimension_-1]);
	for(unsigned i=dimension_-1;i>0;--i){
		index=index*nbin_[i-1]+unsigned(myindices[i-1]);
	}
	return true;
}

bool Grid2::getIndex(const vector<double> & x, unsigned &i) const {
	plumed_assert(x.size()==dimension_);
	vector<int> indices;
	getIndices(x,indices);
	if(getIndex(indices,i)){return true;}else{return false;}
}

//// we are flattening arrays using a column-major order
// if pbc apply, then return the in-box-index
bool Grid2::getIndices(unsigned index,vector<int> &indices) const {
	unsigned kk=index;
	//verify if this is valid
	if(index>maxsize_)return false;
	indices.resize(dimension_);

	indices[0]=index%nbin_[0];
	for(unsigned int i=1;i<dimension_-1;++i){
		kk=(kk-indices[i-1])/nbin_[i-1];
		indices[i]=kk%nbin_[i];
	}
	if(dimension_>=2){
		indices[dimension_-1]=(kk-indices[dimension_-2])/nbin_[dimension_-2];
	}
	// cerr<<" getIndices\n input_index "<<index<<" output_indices ";
	// for(unsigned int i=0;i<dimension_;++i){
	//	cerr<<indices[i]<<" ";
	// }
	// cerr<<"\n";
	// applyPeriodicityToIndices(indices);
	// cerr<<" after_pbc ";
	// for(unsigned int i=0;i<dimension_;++i){
	//	cerr<<indices[i]<<" ";
	// }
	// cerr<<" getIndices_Out"<<endl;
	return true;
}
//
// indices can be positive or negative:
// if periodic this gives the in-the-box indices
void Grid2::getIndices(const vector<double> & x,vector<int> &indices) const {
	plumed_assert(x.size()==dimension_);
	indices.clear();
	for(unsigned i=0;i<dimension_;++i){
		int ind;
		// special case when the pbc:
		// generates a index that can be positive or negative
		// but that, multiplied by i*dx+min is always within the
		// periodicity
		if (pbc_[i]){
			// normalize bewteen -0.5 and 0.5 within periodicity
			double x1=(x[i]-pmin_[i])/pdelta_[i]-0.5;
			// make a int
			x1=x1-int(x1+copysign(0.5,x1));
			// bring back into original periodicity
			x1=(x1+0.5)*pdelta_[i]+pmin_[i];
			ind=floor((x1-min_[i])/dx_[i]);
		}else{
			ind=floor((x[i]-min_[i])/dx_[i]); //note: this allows also out-of-boundary index (negative or larger than the max)
		}
		// apply pbc to the indices
		indices.push_back(ind);
	}
	applyPeriodicityToIndices(indices);// so that the indices will always look in the cell

}
//
//vector<double> Grid2::getPoint(const vector<int> & indices) const {
// plumed_assert(indices.size()==dimension_);
// vector<double> x;
// for(unsigned int i=0;i<dimension_;++i){
//  x.push_back(min_[i]+(double)(indices[i])*dx_[i]);
// }
// return x;
//}
//
//vector<double> Grid2::getPoint(unsigned index) const {
// plumed_assert(index<maxsize_);
// return getPoint(getIndices(index));
//}
//
void Grid2::getPoint(const vector<double> & x,vector<double> &xfloor) const {
	plumed_assert(x.size()==dimension_);
	vector<int> indices;
	getIndices(x,indices);
	getPoint(indices,xfloor);
}

bool Grid2::getPoint(unsigned index,std::vector<double> & point) const{
	if(index>maxsize_){return false;};// check if the index in input is exceeding dimensions
	//cerr<<"getPoint INDEX "<<index<<" \n";
	vector<int> indices;
	if(getIndices(index,indices)){getPoint(indices,point);return true;}{return false;}
}
//
void Grid2::getPoint(const std::vector<int> & indices,std::vector<double> & point) const{
	plumed_assert(indices.size()==dimension_);
	plumed_assert(point.size()==dimension_);
	vector<int> myindices=indices;

	applyPeriodicityToIndices(myindices);

	for(unsigned int i=0;i<dimension_;++i){
		point[i]=(min_[i]+(double)(myindices[i])*dx_[i]);
	}
}

// note that the indices in input can be positive or negative
// and they can be outside the (periodic) allocated grid
// the routine has to check if they are in the grid, otherwise
// it does not append the index
vector<unsigned> Grid2::getNeighbors
(const vector<int> &indices,const vector<unsigned> &nneigh)const{
	plumed_assert(indices.size()==dimension_ && nneigh.size()==dimension_);
	plumed_assert( nneigh.size()==dimension_);

	//cerr<<endl;
	vector<unsigned> neighbors;
	vector<unsigned> small_bin(dimension_);
	unsigned small_nbin=1;
	for(unsigned j=0;j<dimension_;++j){
		small_bin[j]=(2*nneigh[j]+1);
		small_nbin*=small_bin[j];
	}

	vector<unsigned> small_indices(dimension_);
	vector<int> tmp_indices(dimension_);
	for(unsigned index=0;index<small_nbin;++index){
		unsigned kk=index;
		small_indices[0]=(index%small_bin[0]);
		for(unsigned i=1;i<dimension_-1;++i){
			kk=(kk-small_indices[i-1])/small_bin[i-1];
			small_indices[i]=(kk%small_bin[i]);
		}
		if(dimension_>=2){
			small_indices[dimension_-1]=((kk-small_indices[dimension_-2])/small_bin[dimension_-2]);
		}
		for(unsigned i=0;i<dimension_;++i){
			int i0=small_indices[i]-nneigh[i]+indices[i];
			tmp_indices[i]=(i0);
		}
		unsigned myind;
		// note that getIndex  is true if the indices are in the grid
		// cerr<<"Looking for index"<<endl;
		if(getIndex(tmp_indices,myind)){neighbors.push_back(myind);}
		//else{cerr<<"the point is out of the grid\n";};
		// cerr<<endl;

	}
	return neighbors;
}

vector<unsigned> Grid2::getNeighbors
(const vector<double> & x,const vector<unsigned> & nneigh)const{
	plumed_assert(x.size()==dimension_ && nneigh.size()==dimension_);
	//cerr<<"Neighbor (using x)"<<endl;
	vector<int> indices;
	getIndices(x,indices);
	return getNeighbors(indices,nneigh);
}

vector<unsigned> Grid2::getNeighbors
(unsigned index,const vector<unsigned> & nneigh)const{
	//plumed_assert(index<maxsize_ && nneigh.size()==dimension_);
	vector<int> indices;
	vector<unsigned> neighs;
	if(getIndices(index,indices)){
		neighs=getNeighbors(indices,nneigh);
	}
	return neighs;
}

vector<unsigned> Grid2::getSplineNeighbors(const vector<int> & indices)const{
	plumed_assert(indices.size()==dimension_);
	vector<unsigned> neighbors;
	unsigned nneigh=unsigned(pow(2.0,int(dimension_)));

	for(unsigned int i=0;i<nneigh;++i){
		unsigned tmp=i;
		vector<int> nindices;
		for(unsigned int j=0;j<dimension_;++j){
			unsigned i0=tmp%2+indices[j];
			tmp/=2;
			// just push back all, getIndex will take care of periodicity
			nindices.push_back(i0);
		}
		unsigned index;
		if(getIndex(nindices,index)){neighbors.push_back(index);}
	}
	return neighbors;
}

void Grid2::addKernel( const KernelFunctions& kernel ){
	plumed_assert( kernel.ndim()==dimension_ );
	std::vector<unsigned> nneighb=kernel.getSupport( dx_ );
	std::vector<unsigned> neighbors=getNeighbors( kernel.getCenter(), nneighb );
	std::vector<double> xx( dimension_ ); std::vector<Value*> vv( dimension_ );
	std::string str_min, str_max;
	for(unsigned i=0;i<dimension_;++i){
		vv[i]=new Value();
		if( pbc_[i] ){
			Tools::convert(pmin_[i],str_min);
			Tools::convert(pmax_[i],str_max);
			vv[i]->setDomain( str_min, str_max );
		} else {
			vv[i]->setNotPeriodic();
		}
	}

	double newval; std::vector<double> der( dimension_ );
	for(unsigned i=0;i<neighbors.size();++i){
		unsigned ineigh=neighbors[i];
		if(getPoint( ineigh, xx )){
			for(unsigned j=0;j<dimension_;++j) vv[j]->set(xx[j]);
			newval = kernel.evaluate( vv, der, usederiv_ );
			if( usederiv_ ){
				addValueAndDerivatives( ineigh, newval, der );
			}
			else{
				addValue( ineigh, newval );
			}
		}
	}

	for(unsigned i=0;i<dimension_;++i) delete vv[i];
}

bool Grid2::getValue(unsigned index, double &val) const {
	if(index>maxsize_)return false;
	val=grid_[index];return true;
}

double Grid2::getMinValue() const {
	double minval;
	minval=DBL_MAX;
	for(unsigned i=0;i<grid_.size();++i){
		if(grid_[i]<minval)minval=grid_[i];
	}
	return minval;
}

double Grid2::getMaxValue() const {
	double maxval;
	maxval=DBL_MIN;
	for(unsigned i=0;i<grid_.size();++i){
		if(grid_[i]>maxval)maxval=grid_[i];
	}
	return maxval;
}


bool Grid2::getValue(const vector<int> & indices, double &val) const {
	unsigned ind;
	if(getIndex(indices,ind)){getValue(ind,val);return true;}else{return false;};
}

bool Grid2::getValue(const vector<double> & x,double &val) const {
	if(!dospline_){
		unsigned i;if(!getIndex(x,i)){return false;};
		getValue(i,val); // the check is already done before
		return true;
	} else {
		vector<double> der(dimension_);
		if(getValueAndDerivatives(x,val,der)){ return true; } else{return false;}
	}
}

bool Grid2::getValueAndDerivatives
(unsigned index, double &val, vector<double>& der) const{
	if(index>=maxsize_)return false;
	plumed_assert( usederiv_ && der.size()==dimension_);
	der=der_[index];
	val=grid_[index];
	return true;
}

bool Grid2::getValueAndDerivatives
(const vector<int> & indices, double &val, vector<double>& der) const{
	unsigned index;
	if(getIndex(indices,index)){ getValueAndDerivatives(index,val,der); return true;  }else{return false;}
}

bool Grid2::getValueAndDerivatives
(const vector<double> & x, double &value, vector<double>& der) const {
	plumed_assert(der.size()==dimension_ && usederiv_);

	if(dospline_){
		double X,X2,X3;
		vector<double> fd(dimension_);
		vector<double> C(dimension_);
		vector<double> D(dimension_);
		vector<double> dder(dimension_);
		// reset
		value=0.0;
		for(unsigned int i=0;i<dimension_;++i) der[i]=0.0;

		vector<int> indices;getIndices(x,indices);
		vector<unsigned> neigh=getSplineNeighbors(indices);
		// strict rule: if the number of neighbors is not complete turn there is no
		// way to compute it
		if(neigh.size()<unsigned(pow(2.0,int(dimension_)))){return false;};
		vector<double>   xfloor(dimension_);getPoint(x,xfloor);


		// loop over neighbors
		for(unsigned int ipoint=0;ipoint<neigh.size();++ipoint){
			double grid;
			if(!getValueAndDerivatives(neigh[ipoint],grid,dder)){return false;};
			vector<int> nindices;
			getIndices(neigh[ipoint],nindices);
			double ff=1.0;

			for(unsigned j=0;j<dimension_;++j){
				int x0=1;
				if(nindices[j]==indices[j]) x0=0;
				double dx=getDx()[j];
				X=fabs((x[j]-xfloor[j])/dx-(double)x0);
				X2=X*X;
				X3=X2*X;
				double yy;
				if(fabs(grid)<0.0000001) yy=0.0;
				else yy=-dder[j]/grid;
				C[j]=(1.0-3.0*X2+2.0*X3) - (x0?-1.0:1.0)*yy*(X-2.0*X2+X3)*dx;
				D[j]=( -6.0*X +6.0*X2) - (x0?-1.0:1.0)*yy*(1.0-4.0*X +3.0*X2)*dx;
				D[j]*=(x0?-1.0:1.0)/dx;
				ff*=C[j];
			}
			for(unsigned j=0;j<dimension_;++j){
				fd[j]=D[j];
				for(unsigned i=0;i<dimension_;++i) if(i!=j) fd[j]*=C[i];
			}
			value+=grid*ff;
			for(unsigned j=0;j<dimension_;++j) der[j]+=grid*fd[j];
		}

	}else{

		unsigned index;
		if(!getIndex(x,index))return false;
		getValueAndDerivatives(index,value,der);

	}
	return true;
}

bool Grid2::setValue(const unsigned index, const double value){
	if(index<maxsize_){
		grid_[index]=value;
		plumed_massert(!usederiv_,"you cannot set the value in this grid without setting the derivative");
		return true;
	}else{return false;}
}

bool Grid2::setValue(const vector<int> & indices, double value){
	unsigned ind; if(getIndex(indices,ind)){ setValue(ind,value);return true;}else{return false;}
}

bool Grid2::setValueAndDerivatives
(unsigned index, double value, vector<double>& der){
	if(index>=maxsize_){return false;}else{
		plumed_assert(usederiv_ && der.size()==dimension_);
		grid_[index]=value;
		der_[index]=der;
		return true;
	}
}

bool Grid2::setValueAndDerivatives
(const vector<int> & indices, double value, vector<double>& der){
	unsigned ind;
	if(getIndex(indices,ind)){
		setValueAndDerivatives(ind,value,der);
		return true;
	}else{ return false;}
}

bool Grid2::setValueAndDerivatives
(const vector<double> & x, double value, vector<double>& der){
	unsigned ind;
	if(getIndex(x,ind)){
		setValueAndDerivatives(ind,value,der);
		return true;
	}else{ return false;}
}

bool Grid2::addValue(unsigned index, double value){
	plumed_assert( !usederiv_);
	if(index<maxsize_){
		grid_[index]+=value;return true;
	}else{return false;}
}

bool Grid2::addValue(const vector<int> & indices, double value){
	unsigned ind;if(getIndex(indices,ind)){
		addValue(ind,value);
		return true;
	}else{return false;}
}

// return false if the index is not acceptable
bool Grid2::addValueAndDerivatives
(unsigned index, double value, vector<double>& der){
	if(index<maxsize_){
		plumed_assert(usederiv_ && der.size()==dimension_);
		grid_[index]+=value;
		for(unsigned int i=0;i<dimension_;++i) der_[index][i]+=der[i];
		return true;
	}else{return false;}
}

bool Grid2::addValueAndDerivatives
(const vector<int> & indices, double value, vector<double>& der){
	unsigned ind;
	if(getIndex(indices,ind)){
		addValueAndDerivatives(ind,value,der);
		return true;
	}else{return false;}
}

void Grid2::scaleAllValuesAndDerivatives( const double& scalef ){
	if(usederiv_){
		for(unsigned i=0;i<grid_.size();++i){
			grid_[i]*=scalef;
			for(unsigned j=0;j<dimension_;++j) der_[i][j]*=scalef;
		}
	} else {
		for(unsigned i=0;i<grid_.size();++i) grid_[i]*=scalef;
	}
}

void Grid2::applyFunctionAllValuesAndDerivatives( double (*func)(double val), double (*funcder)(double valder) ){
	if(usederiv_){
		for(unsigned i=0;i<grid_.size();++i){
			grid_[i]=func(grid_[i]);
			for(unsigned j=0;j<dimension_;++j) der_[i][j]=funcder(der_[i][j]);
		}
	} else {
		for(unsigned i=0;i<grid_.size();++i) grid_[i]=func(grid_[i]);
	}
}

void Grid2::writeHeader(OFile& ofile){
	for(unsigned i=0;i<dimension_;++i){
		ofile.addConstantField("min_" + argnames[i]);
		ofile.addConstantField("max_" + argnames[i]);
		ofile.addConstantField("nbins_" + argnames[i]);
		ofile.addConstantField("periodic_" + argnames[i]);
		if(pbc_[i]){
			ofile.addConstantField("pmin_" + argnames[i]);
			ofile.addConstantField("pmax_" + argnames[i]);
		}
	}
}

void Grid2::writeToFile(OFile& ofile){
	vector<double> xx(dimension_);
	vector<double> der(dimension_);
	writeHeader(ofile);
	//cerr<<"FILEWRITING\n";
	for(unsigned i=0;i<getSize();++i){
		plumed_massert(getPoint(i,xx),"Trying to print an out of grid point");
		double  f;
		if(usederiv_){plumed_massert(getValueAndDerivatives(i,f,der),"trying to write out of the grid");}
		else{ plumed_massert(getValue(i,f),"trying to write out of the grid");};
		vector<int> indices; plumed_massert(getIndices(i,indices),"the index is not there");
		if(i>0 && dimension_>1 && indices[dimension_-2]==0) ofile.printf("\n");
		for(unsigned j=0;j<dimension_;++j){
			ofile.printField("min_" + argnames[j], str_min_[j] );
			ofile.printField("max_" + argnames[j], str_max_[j] );
			ofile.printField("nbins_" + argnames[j], static_cast<int>(nbin_[j]) );
			if( pbc_[j] ) {
				ofile.printField("periodic_" + argnames[j],"true");
				ofile.printField("pmin_" + argnames[j], str_pmin_[j] );
				ofile.printField("pmax_" + argnames[j], str_pmax_[j] );
			}
			else          ofile.printField("periodic_" + argnames[j], "false" );
		}
		for(unsigned j=0;j<dimension_;++j){ ofile.fmtField(" "+fmt_); ofile.printField(argnames[j],xx[j]); }
		ofile.fmtField(" "+fmt_); ofile.printField(funcname,f);

		if(usederiv_) for(unsigned j=0;j<dimension_;++j){ ofile.fmtField(" "+fmt_); ofile.printField("der_" + argnames[j] ,der[j]); }
		ofile.printField();
	}
}
// read the grid from input file
Grid2* Grid2::create(const std::string& funcl, std::vector<Value*> args, IFile& ifile,
		const vector<std::string> & gmin,const vector<std::string> & gmax,
		const vector<unsigned> & nbin,bool dosparse, bool dospline, bool doder){
	Grid2* grid=Grid2::create(funcl,args,ifile,dosparse,dospline,doder);
	std::vector<unsigned> cbin( grid->getNbin() );
	std::vector<std::string> cmin( grid->getMin() ), cmax( grid->getMax() );
	for(unsigned i=0;i<args.size();++i){

		double gmax_d,gmin_d,cmax_d,cmin_d;
		Tools::convert(gmax[i],gmax_d);
		Tools::convert(gmin[i],gmin_d);
		Tools::convert(cmin[i],cmin_d);
		Tools::convert(cmax[i],cmax_d);

		plumed_massert( pow(pow(cmax_d-gmax_d,2),0.5)<1.e-6 ||  pow(pow(cmin_d-cmin_d,2),0.5)<1.e-6 , "boundaries of the input grid and output do not coincide" );


		if( args[i]->isPeriodic() ){
			plumed_massert( cbin[i]==nbin[i], "mismatched grid nbins" );
		} else {
			plumed_massert( (cbin[i]-1)==nbin[i], "mismatched grid nbins");
		}

	}
	return grid;
}

Grid2* Grid2::create(const std::string& funcl, std::vector<Value*> args, IFile& ifile, bool dosparse, bool dospline, bool doder)
{
	Grid2* grid=NULL;
	unsigned nvar=args.size(); bool hasder=false; std::string pstring;
	std::vector<int> gbin1(nvar); std::vector<unsigned> gbin(nvar);
	std::vector<std::string> labels(nvar),gmin(nvar),gmax(nvar),permin(nvar),permax(nvar);
	std::vector<std::string> fieldnames;
	ifile.scanFieldList( fieldnames );
	// Retrieve names for fields

	for(unsigned i=0;i<args.size();++i) labels[i]=args[i]->getName();
	// And read the stuff from the header


	for(unsigned i=0;i<args.size();++i){
		ifile.scanField( "min_" + labels[i], gmin[i]);
		ifile.scanField( "max_" + labels[i], gmax[i]);
		ifile.scanField( "periodic_" + labels[i], pstring );
		ifile.scanField( "nbins_" + labels[i], gbin1[i]);
		plumed_assert( gbin1[i]>0 );
		if( args[i]->isPeriodic() ){
			plumed_massert( pstring==string("true"), "input value is periodic but grid is not");
			std::string pmin, pmax;
			args[i]->getDomain( pmin, pmax );
			ifile.scanField( "pmin_" + labels[i], permin[i]);
			ifile.scanField( "pmax_" + labels[i], permax[i]);
			//convert the boundary just read
			double permin_d,permax_d,pmin_d,pmax_d;
			Tools::convert(permin[i],permin_d);
			Tools::convert(permax[i],permax_d);
			Tools::convert(pmax,pmax_d);
			Tools::convert(pmin,pmin_d);
			// now check if the boundaries are the same
			//if( pmin!=gmin[i] || pmax!=permax[i] )
			plumed_massert( pow(pow(pmin_d-permin_d,2),0.5)<1.e-6 ||  pow(pow(pmax_d-permax_d,2),0.5)<1.e-6 , "boundaries and values do not coindcide" );
		} else {
			gbin[i]=gbin1[i]-1;
			//plumed_massert( pstring=="true", "grid value is periodic but value is not");
			plumed_massert(string("true")!=pstring,"grid value is periodic but value is not");

		}
		hasder=ifile.FieldExist( "der_" + args[i]->getName() );
		if( doder && !hasder ) plumed_merror("missing derivatives from grid file");

		for(unsigned j=0;j<fieldnames.size();++j){
			for(unsigned k=i+1;k<args.size();++k){
				if( fieldnames[j]==labels[k] ) plumed_merror("arguments in input are not in same order as in grid file");
			}
			if( fieldnames[j]==labels[i] ) break;
		}

	}
	// disable the rescaling of the limits (last field bool)
	if(!dosparse){grid=new Grid2(funcl,args,gmin,gmax,gbin,dospline,doder,true,false);}
	else{grid=new SparseGrid2(funcl,args,gmin,gmax,gbin,dospline,doder);}

	vector<double> xx(nvar),dder(nvar);
	vector<double> dx=grid->getDx();
	double f,x;
	while( ifile.scanField(funcl,f) ){
		for(unsigned i=0;i<nvar;++i){
			ifile.scanField(labels[i],x); xx[i]=x+dx[i]/2.0;
			ifile.scanField( "min_" + labels[i], gmin[i]);
			ifile.scanField( "max_" + labels[i], gmax[i]);
			ifile.scanField( "nbins_" + labels[i], gbin1[i]);
			ifile.scanField( "periodic_" + labels[i], pstring );
			if(pstring=="true"){
				ifile.scanField( "pmin_" + labels[i], permin[i]);
				ifile.scanField( "pmax_" + labels[i], permax[i]);
			}
		}
		if(hasder){ for(unsigned i=0;i<nvar;++i){ ifile.scanField( "der_" + args[i]->getName(), dder[i] ); } }
		unsigned index;
		plumed_massert(grid->getIndex(xx,index),"the index you asked is out of grid");
		if(doder){grid->setValueAndDerivatives(index,f,dder);}
		else{grid->setValue(index,f);}
		ifile.scanField();
	}
	return grid;
}

// Sparse version of grid with map
void SparseGrid2::clear(){
	map_.clear();
}

unsigned SparseGrid2::getSize() const{
	return map_.size();
}

unsigned SparseGrid2::getMaxSize() const {
	return maxsize_;
}

bool SparseGrid2::getValue(unsigned index, double &val)const{
	if(index<maxsize_){
		val=0.0;
		iterator it=map_.find(index);
		if(it!=map_.end()) val=it->second;
		return true;
	}else{return false;}
}




bool SparseGrid2::getValueAndDerivatives
(unsigned index, double &value, vector<double>& der)const{
	if(index>=maxsize_){return false;};
	plumed_assert(usederiv_ && der.size()==dimension_);
	value=0.0;
	for(unsigned int i=0;i<dimension_;++i) der[i]=0.0;
	iterator it=map_.find(index);
	if(it!=map_.end()) value=it->second;
	iterator_der itder=der_.find(index);
	if(itder!=der_.end()) der=itder->second;
	return true;
}

bool SparseGrid2::setValue(unsigned index, double &value){
	plumed_assert(!usederiv_);
	if(index<maxsize_){
		map_[index]=value;
		return true;
	}else{return false;}
}

bool SparseGrid2::setValueAndDerivatives
(unsigned index, double &value, vector<double>& der){
	plumed_assert(usederiv_ && der.size()==dimension_);
	if(index<maxsize_){
		map_[index]=value;
		der_[index]=der;
		return true;
	}else{return false;}
}

bool SparseGrid2::addValue(unsigned index, double &value){
	plumed_assert(!usederiv_);
	if(index<maxsize_){
		map_[index]+=value;
		return true;
	}else{return false;}
}

bool SparseGrid2::addValueAndDerivatives
(unsigned index,  double value, vector<double>& der ){
	plumed_assert(usederiv_ && der.size()==dimension_);
	if(index<maxsize_){
		map_[index]+=value;
		der_[index].resize(dimension_);
		for(unsigned int i=0;i<dimension_;++i) der_[index][i]+=der[i];
		return true;
	}else{return false;}
}

void SparseGrid2::writeToFile(OFile& ofile){
	vector<double> xx(dimension_);
	vector<double> der(dimension_);
	double f;
	writeHeader(ofile);
	ofile.fmtField(" "+fmt_);
	for(iterator it=map_.begin();it!=map_.end();++it){
		unsigned i=(*it).first;
		if(!getPoint(i,xx)){cerr<<"the point does not exist";};
		if(usederiv_){getValueAndDerivatives(i,f,der);}
		else{getValue(i,f);}
		vector<int> indices(dimension_);getIndices(i,indices);
		if(i>0 && dimension_>1 && indices[dimension_-2]==0) ofile.printf("\n");
		for(unsigned j=0;j<dimension_;++j){
			ofile.printField("min_" + argnames[j], str_min_[j] );
			ofile.printField("max_" + argnames[j], str_max_[j] );
			ofile.printField("nbins_" + argnames[j], static_cast<int>(nbin_[j]) );
			if( pbc_[j] ) {
				ofile.printField("periodic_" + argnames[j],"true");
				ofile.printField("pmin_" + argnames[j], str_pmin_[j] );
				ofile.printField("pmax_" + argnames[j], str_pmax_[j] );
			}
			else          ofile.printField("periodic_" + argnames[j], "false" );
		}
		for(unsigned j=0;j<dimension_;++j) ofile.printField(argnames[j],xx[j]);
		ofile.printField(funcname, f);
		if(usederiv_){ for(unsigned j=0;j<dimension_;++j) ofile.printField("der_" + argnames[j],der[j]); }
		ofile.printField();
	}
}


void Grid2::projectOnLowDimension(double &val, std::vector<int> &vHigh, WeightBase * ptr2obj ){
	unsigned i=0;
	for(i=0;i<vHigh.size();i++){
		if(vHigh[i]<0){// this bin needs to be integrated out
			// parallelize here???
			for(unsigned j=0;j<(getNbin())[i];j++){
				vHigh[i]=int(j);
				projectOnLowDimension(val,vHigh,ptr2obj); // recursive function: this is the core of the mechanism
				vHigh[i]=-1;
			}
			return; //
		}
	}
	// when there are no more bin to dig in then retrieve the value
	if(i==vHigh.size()){
		//std::cerr<<"POINT: ";
		//for(unsigned j=0;j<vHigh.size();j++){
		//   std::cerr<<vHigh[j]<<" ";
		//}
		std::vector<int> vv(vHigh.size());
		for(unsigned j=0;j<vHigh.size();j++)vv[j]=int(vHigh[j]);
		//
		// this is the real assignment !!!!! (hack this to have bias or other stuff)
		//

		// this case: produce fes
		//val+=exp(beta*getValue(vv)) ;

		double myv; plumed_massert(getValue(vv,myv),"this bin is missing");
		val=ptr2obj->projectInnerLoop(val,myv) ;
		// to be added: bias (same as before without negative sign)
		//std::cerr<<" VAL: "<<val <<endl;
	}
}

Grid2 Grid2::project(const std::vector<std::string> & proj , WeightBase *ptr2obj ){
	// find extrema only for the projection
	vector<string>   smallMin,smallMax;
	vector<string>   smallpMin,smallpMax;
	vector<unsigned> smallBin;
	vector<unsigned> dimMapping;
	vector<bool> smallIsPeriodic;
	vector<string> smallName;

	// check if the two key methods are there
	WeightBase* pp = dynamic_cast<WeightBase*>(ptr2obj);
	if (!pp)plumed_merror("This WeightBase is not complete: you need a projectInnerLoop and projectOuterLoop ");

	for(unsigned j=0;j<proj.size();j++){
		for(unsigned i=0;i<getArgNames().size();i++){
			if(proj[j]==getArgNames()[i]){
				unsigned offset;
				// note that at sizetime the non periodic dimension get a bin more
				if(getIsPeriodic()[i]){offset=0;}else{offset=1;}
				smallMax.push_back(getMax()[i]);
				smallMin.push_back(getMin()[i]);
				smallpMax.push_back(getPeriodicMax()[i]);
				smallpMin.push_back(getPeriodicMin()[i]);
				smallBin.push_back(getNbin()[i]-offset);
				smallIsPeriodic.push_back(getIsPeriodic()[i]);
				dimMapping.push_back(i);
				smallName.push_back(getArgNames()[i]);
				break;
			}
		}
	}
	Grid2 smallgrid("projection",smallName,smallMin,smallMax,smallBin,false,false,true,smallIsPeriodic,smallpMin,smallpMax,false);
	// check that the two grids are commensurate
	for(unsigned i=0;i<dimMapping.size();i++){

		plumed_massert(  (smallgrid.getMax())[i] == (getMax())[dimMapping[i]],  "the two input grids are not compatible in max"   );
		plumed_massert(  (smallgrid.getMin())[i] == (getMin())[dimMapping[i]],  "the two input grids are not compatible in min"   );
		plumed_massert(  (smallgrid.getNbin())[i]== (getNbin())[dimMapping[i]], "the two input grids are not compatible in bin"   );
	}
	vector<unsigned> toBeIntegrated;
	for(unsigned i=0;i<getArgNames().size();i++){
		bool doappend=true;
		for(unsigned j=0;j<dimMapping.size();j++){
			if(dimMapping[j]==i){doappend=false; break;}
		}
		if(doappend)toBeIntegrated.push_back(i);
	}
//	for(unsigned i=0;i<dimMapping.size();i++ ){
//	     cerr<<"Dimension to preserve "<<dimMapping[i]<<endl;
//	}
	//for(unsigned i=0;i<toBeIntegrated.size();i++ ){
	//     cerr<<"Dimension to integrate "<<toBeIntegrated[i]<<endl;
	//}

	// loop over all the points in the Grid2, find the corresponding fixed index, rotate over all the other ones
	for(unsigned i=0;i<smallgrid.getSize();i++){
		std::vector<int> v;
		plumed_massert(smallgrid.getIndices(i,v),"this index is not valid in smallgrid");
		std::vector<int> vHigh((getArgNames()).size(),-1);
		for(unsigned j=0;j<dimMapping.size();j++)vHigh[dimMapping[j]]=int(v[j]);
		// the vector vhigh now contains the index of the low dimension and -1 in place of the dimensions that need to be integrated
		double initval=0.;
		projectOnLowDimension(initval,vHigh, ptr2obj);
		plumed_massert(smallgrid.setValue(i,initval),"cannot set the value");
	}
	// reset to zero just for biasing (this option can be evtl enabled in a future...)
	//double vmin;vmin=-smallgrid.getMinValue()+1;
	for(unsigned i=0;i<smallgrid.getSize();i++){
		//         //if(dynamic_cast<BiasWeight*>(ptr2obj)){
		//         //        smallgrid.addValue(i,vmin);// go to 1
		//         //}
		double vv;plumed_massert(smallgrid.getValue(i,vv),"you do not have the grid point");
		smallgrid.setValue(i,ptr2obj->projectOuterLoop(vv));
		//         //if(dynamic_cast<BiasWeight*>(ptr2obj)){
		//         //        smallgrid.addValue(i,-vmin);// bring back to the value
		//         //}
	}

	return smallgrid;
}
/// here test all the things you want to
int Grid2::gridTester(){
	cerr<<"ENTERING GRIDTESTER"<<endl;
	string root=config::getPlumedRoot();
	cerr<<"ROOT IS "<<root<<endl;

	//	cerr<<"create a grid\n";
	//	vector<string> vars;
	//	vars.push_back(string("var1"));
	//	vars.push_back(string("var2"));
	//	vector<string> gmin;
	//	gmin.push_back(string("2.5"));
	//	gmin.push_back(string("-10"));
	//	vector<string> gmax;
	//	gmax.push_back(string("-1."));
	//	gmax.push_back(string("20"));
	//	vector<string> pmin;
	//	pmin.push_back(string("-pi"));
	//	pmin.push_back(string(""));
	//	vector<string> pmax;
	//	pmax.push_back(string("pi"));
	//	pmax.push_back(string(""));
	//	vector<unsigned> nbin;nbin.push_back(40);nbin.push_back(80);
	//	vector<bool> isperiodic;isperiodic.push_back(true);isperiodic.push_back(false);
	//
	//	Grid2 mygrid(string("test"),vars,gmin,gmax,nbin,false,true,true,isperiodic,pmin,pmax,true);
	//
	//	cerr<<" now try to get index for a point "<<endl;
	//
	//	vector<double> x;x.push_back(-1.);x.push_back(3.);
	//	vector<int> indices;
	//	mygrid.getIndices(x,indices);
	//	cerr<<"The point "<<x[0]<<" "<<x[1]<<" corresponds to index "<<indices[0]<<" "<<indices[1]<<endl;
	//
	//	cerr<<" now a grid with reversed boundary in periodic "<<endl;
	//
	//	gmin[0]=string("2.5");gmax[0]=string("-1.0");
	//	Grid2 mygrid2(string("test"),vars,gmin,gmax,nbin,false,true,true,isperiodic,pmin,pmax,true);
	//
	//	// now add a value
	//	for(unsigned i=0;i<100;i++){
	//		vector<double> xx;xx.resize(2);
	//		xx[0]=-3.14+6.28*i/30.;
	//		//xx[1]=double(i);
	//		xx[1]=0.;
	//		double v=sin(x[0])+xx[1];
	//		vector<double> vv;vv.resize(2);
	//		vv[0]=cos(x[0]);
	//		vv[1]=1;
	//		vector<int> indices;
	//		// note that the index is always cleaned up
	//		indices.resize(100);
	//		mygrid2.getIndices(xx,indices);
	//		if(!mygrid2.setValueAndDerivatives(xx,v,vv)){
	//			cerr<<"POINT  "<<setw(4)<<xx[0]<<" "<<setw(4)<<xx[1]<<" IS OUT OF THE GRID. IND "<<setw(4)<<indices[0]<<" "<<setw(4)<<indices[1]<<endl;
	//		}else{
	//			cerr<<"POINT  "<<setw(4)<<xx[0]<<" "<<setw(4)<<xx[1]<<" IS IN THE GRID. IND "<<setw(4)<<indices[0]<<" "<<setw(4)<<indices[1]<<endl;
	//		};
	//	}
	//	OFile gridfile;
	//	gridfile.open(string("testout.dat"));
	//	mygrid2.setOutputFmt(string("%14.9f"));
	//	mygrid2.writeToFile(gridfile);
	//	gridfile.close();
	//
	//	cerr<<"sparsegrids are grids that are made with associative arrays. Might be not completely full at the beginning\n";
	//
	//	SparseGrid2 mygrid3(string("test"),vars,gmin,gmax,nbin,false,true,isperiodic,pmin,pmax,true);

	//
	cerr<<"create a grid\n";
	vector<string> vars;
	vars.push_back(string("var1"));
	vars.push_back(string("var2"));
	vector<string> gmin;
	gmin.push_back(string("-pi"));
	gmin.push_back(string("-pi"));
	vector<string> gmax;
	gmax.push_back(string("pi"));
	gmax.push_back(string("pi"));
	vector<string> pmin;
	pmin.push_back(string("-pi"));
	pmin.push_back(string("-pi"));
	vector<string> pmax;
	pmax.push_back(string("pi"));
	pmax.push_back(string("pi"));
	vector<bool> isperiodic;isperiodic.push_back(true);isperiodic.push_back(true);

	//	// non periodic
	//	pmin[0]="";pmin[1]="";
	//	pmax[0]="";pmax[1]="";
	//	isperiodic[0]=false;
	//	isperiodic[1]=false;


	cerr<<fixed;
	vector<unsigned> nbin;nbin.push_back(500);nbin.push_back(500);
	clock_t start = clock(), diff;
	cerr<<"STANDARD GRID: allocate 500x500 (with clean)\n";
	string ss="grid.dat";
	Grid *mGrid=new Grid(ss,vars,gmin,gmax,nbin,false,true,true,isperiodic,pmin,pmax);
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	cerr<<"NEW GRID: allocate 500x500 (with clean)\n";
	start = clock();
	ss="grid2.dat";
	Grid2 *mGrid2=new Grid2(ss,vars,gmin,gmax,nbin,false,true,true,isperiodic,pmin,pmax);
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	vector<double> xx;xx.push_back(0.);xx.push_back(0.);
	vector<double> der;der.push_back(0.);der.push_back(0.);
	double bias=0.;
	unsigned ineigh;
	cerr<<"add 100000 values on old grid\n";
	start = clock();
	for(unsigned ii=0;ii<100000;ii++){
		mGrid->getPoint(ineigh,xx);
		mGrid->addValueAndDerivatives(ineigh,bias,der);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	cerr<<"add 100000 values on new grid\n";
	start = clock();
	for(unsigned ii=0;ii<100000;ii++){
		mGrid2->getPoint(ineigh,xx);
		mGrid2->addValueAndDerivatives(ineigh,bias,der);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	cerr<<"retrieve 100000 values on old grid\n";
	start = clock();
	for(unsigned ii=0;ii<100000;ii++){
		mGrid->getPoint(ineigh,xx);
		bias=mGrid->getValueAndDerivatives(ineigh,der);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	cerr<<"retrieve 100000 values on new grid\n";
	start = clock();
	for(unsigned ii=0;ii<100000;ii++){
		mGrid2->getPoint(ineigh,xx);
		mGrid2->getValueAndDerivatives(ineigh,bias,der);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	cerr<<"retrieve 100000 indexes on old grid\n";
	start = clock();
	for(unsigned ii=0;ii<100000;ii++){
		mGrid->getPoint(ineigh,xx);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	cerr<<"retrieve 100000 indexes on new grid\n";
	start = clock();
	for(unsigned ii=0;ii<100000;ii++){
		mGrid2->getPoint(ineigh,xx);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	vector<unsigned> indices;
	indices.push_back(0);    indices.push_back(0);
	cerr<<"retrieve 100000 point from indices on old grid\n";
	start = clock();
	for(unsigned ii=0;ii<100000;ii++){
		mGrid->getPoint(indices,xx);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	cerr<<"retrieve 100000 point from indices on new grid\n";
	start = clock();
	vector<int> indices2;
	indices2.push_back(0);    indices2.push_back(0);
	for(unsigned ii=0;ii<100000;ii++){
		mGrid2->getPoint(indices2,xx);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;\


	cerr<<"retrieve 100000 index from indices on old grid\n";
	start = clock();
	unsigned jj;
	for(unsigned ii=0;ii<100000;ii++){
		jj=mGrid->getIndex(indices);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;
	cerr<<"retrieve 100000 index from indices on new grid\n";
	start = clock();
	for(unsigned ii=0;ii<100000;ii++){
		mGrid2->getIndex(indices2,jj);
	};
	diff=clock()-start;
	cerr<<" TIME : "<<diff * 1000 / CLOCKS_PER_SEC<<" MS "<<endl;


	cerr<<"EXITING GRIDTESTER"<<endl;

	plumed_merror("dying because of end of tests");
	return 0;
}

}
