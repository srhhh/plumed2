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
#include "RMSD.h"
#include "PDB.h"
#include "Log.h"
#include "OptimalAlignment.h"
#include "Exception.h"
#include <cmath>
#include <iostream>

using namespace std;
namespace PLMD{

RMSD::RMSD(Log & log ):
  alignmentMethod(SIMPLE),
  myoptimalalignment(NULL),
  log(&log){}

RMSD& RMSD::operator=(const RMSD& v){
  alignmentMethod=v.alignmentMethod;
  reference=v.reference;
  align=v.align;
  displace=v.displace;
// in this manner the new RMSD is built empty and will just allocate its own
// myoptimalalignment when used (in calculate())
  myoptimalalignment=NULL;

  log=v.log;
  return *this ;
}


RMSD::RMSD(const RMSD & oldrmsd):
  alignmentMethod(oldrmsd.alignmentMethod),
  reference(oldrmsd.reference),
  align(oldrmsd.align),
  displace(oldrmsd.align),
// in this manner the new RMSD is built empty and will just allocate its own
// myoptimalalignment when used (in calculate())
  myoptimalalignment(NULL),
  log( oldrmsd.log )
  {  }

void RMSD::set(const PDB&pdb, string mytype ){

	setReference(pdb.getPositions());
	setAlign(pdb.getOccupancy());
	setDisplace(pdb.getBeta());
        setType(mytype);
}

void RMSD::setType(string mytype){
	myoptimalalignment=NULL;

	alignmentMethod=SIMPLE; // initialize with the simplest case: no rotation
	if (mytype=="SIMPLE"){
		alignmentMethod=SIMPLE;
		log->printf("RMSD IS DONE WITH SIMPLE METHOD(NO ROTATION)\n")
	;}
	else if (mytype=="OPTIMAL"){
		alignmentMethod=OPTIMAL;
		log->printf("RMSD IS DONE WITH OPTIMAL ALIGNMENT METHOD\n");
	}
	else plumed_merror("unknown RMSD type" + mytype);

}

void RMSD::clear(){
  reference.clear();
  align.clear();
  displace.clear();
}
RMSD::~RMSD(){
	if(myoptimalalignment!=NULL) delete myoptimalalignment;
}

string RMSD::getMethod(){
	string mystring;
	switch(alignmentMethod){
		case SIMPLE: mystring.assign("SIMPLE");break; 
		case OPTIMAL: mystring.assign("OPTIMAL");break; 
	}	
	return mystring;
}

void RMSD::setReference(const vector<Vector> & reference){
  unsigned n=reference.size();
  this->reference=reference;
  plumed_massert(align.empty(),"you should first clear() an RMSD object, then set a new referece");
  plumed_massert(displace.empty(),"you should first clear() an RMSD object, then set a new referece");
  align.resize(n,1.0);
  displace.resize(n,1.0);
}

void RMSD::setAlign(const vector<double> & align){
  plumed_massert(this->align.size()==align.size(),"mismatch in dimension of align/displace arrays");
  this->align=align;
}

void RMSD::setDisplace(const vector<double> & displace){
  plumed_massert(this->displace.size()==displace.size(),"mismatch in dimension of align/displace arrays");
  this->displace=displace;
}

double RMSD::calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, bool squared){

  double ret=0.;

  switch(alignmentMethod){
	case SIMPLE:
		//	do a simple alignment without rotation 
		ret=simpleAlignment(align,displace,positions,reference,log,derivatives,squared);
		break;	
	case OPTIMAL:
		if (myoptimalalignment==NULL){ // do full initialization	
			//
			// I create the object only here
			// since the alignment object require to know both position and reference
			// and it is possible only at calculate time
			//
			myoptimalalignment=new OptimalAlignment(align,displace,positions,reference,log);
        }
		// this changes the P0 according the running frame
		(*myoptimalalignment).assignP0(positions);

		ret=(*myoptimalalignment).calculate(squared, derivatives);
		//(*myoptimalalignment).weightedFindiffTest(false);

		break;	
  }	

  return ret;

}

double RMSD::simpleAlignment(const  std::vector<double>  & align,
		                     const  std::vector<double>  & displace,
		                     const std::vector<Vector> & positions,
		                     const std::vector<Vector> & reference ,
		                     Log* &log,
		                     std::vector<Vector>  & derivatives, bool squared) {
      double dist(0);
      double anorm(0);
      double dnorm(0);
      unsigned n=reference.size();

      Vector apositions;
      Vector areference;
      Vector dpositions;
      Vector dreference;

      for(unsigned i=0;i<n;i++){
        double aw=align[i];
        double dw=displace[i];
        anorm+=aw;
        dnorm+=dw;
        apositions+=positions[i]*aw;
        areference+=reference[i]*aw;
        dpositions+=positions[i]*dw;
        dreference+=reference[i]*dw;
      }

      apositions/=anorm;
      areference/=anorm;
      dpositions/=dnorm;
      dreference/=dnorm;

      double invdnorm=1.0/dnorm;
      double invanorm=1.0/anorm;
      Vector shift=((apositions-areference)-(dpositions-dreference));
      for(unsigned i=0;i<n;i++){
        Vector d=(positions[i]-apositions)-(reference[i]-areference);
        dist+=displace[i]*d.modulo2();
        derivatives[i]=2*(invdnorm*displace[i]*d+invanorm*align[i]*shift);
      }
     dist*=invdnorm;

     if(!squared){
	// sqrt and normalization
        dist=sqrt(dist);
	///// sqrt and normalization on derivatives
        for(unsigned i=0;i<n;i++){derivatives[i]*=(0.5/dist);}
      }
      return dist;
}
}
