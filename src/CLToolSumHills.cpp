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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "Tools.h"
#include "Plumed.h"
#include "PlumedCommunicator.h"
#include "Random.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include "PlumedFile.h"
#include "Value.h"
#include "Matrix.h"
#include "Gaussian.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC TOOLS sum_hills 
/*
driver is a tool that allows one to to use plumed to post-process an existing trajectory.


*/
//+ENDPLUMEDOC

class CLToolSumHills : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  CLToolSumHills(const CLToolOptions& co );
  int main(FILE* in,FILE*out,PlumedCommunicator& pc);
  string description()const;
};

void CLToolSumHills::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys );
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  //keys.add("compulsory","--plumed","plumed.dat","specify the name of the plumed input file");
  keys.add("compulsory","--hills","HILLS","specify the name of the hills file");
}

CLToolSumHills::CLToolSumHills(const CLToolOptions& co ):
CLTool(co)
{
 inputdata=commandline;
}

string CLToolSumHills::description()const{ return "sum the hills with  plumed"; }


int CLToolSumHills::main(FILE* in,FILE*out,PlumedCommunicator& pc){
  cerr<<"sum_hills utility  "<<endl;
  
// Read the hills input file name  
  string hillsFile; parse("--hills",hillsFile);
  cerr<<"FILENAME "<<hillsFile<<endl;  
// create the plumed object first
   Plumed p;
// parse it as it was a restart
   PlumedIFile *ifile=new PlumedIFile();
   vector<string> myfields;
   vector<Gaussian> hills;
   if(ifile->FileExist(hillsFile)){
      ifile->open(hillsFile);
      double dummy;
      ifile->scanFieldList(myfields);
      //for(int i=0;i<myfields.size();i++)cerr<<myfields[i]<<endl; 
      // now find first sigma 
      size_t found;
      vector<bool> issigma;issigma.resize(myfields.size());
      for(int i=0;i<myfields.size();i++){
         size_t pos = 0;
	found=myfields[i].find("sigma_", pos);
        if (found!=string::npos){
            issigma[i]=true;
        }else{
            issigma[i]=false;
        }
      }
      int ncv=0;
      std::vector<Value> tmpvalues; tmpvalues.resize(ncv);
      for(int i=1;i<myfields.size();i++){
        if (issigma[i])break;
        cerr<<"found variables "<<myfields[i]<<endl; 
        tmpvalues.push_back( Value() ); 
        tmpvalues.back().setName(myfields[i]);
        ncv++;
      }
      vector<bool> isconstant;
      bool multivariate=false;
      vector<double> center(ncv);
      vector<double> sigma(ncv);
      double height;
      int nhills=0; 
      vector<Gaussian> hills_;
      while(ifile->scanField("time",dummy)){
	for(unsigned i=0;i<ncv;++i){
   	 ifile->scanField( &tmpvalues[i] );
         center[i]=tmpvalues[i].get();
      	}
        // sigmas 
        std::string sss;
        ifile->scanField("multivariate",sss);
        if(sss=="true") multivariate=true;
        else if(sss=="false") multivariate=false;
        else plumed_merror("cannot parse multivariate = "+ sss);
	 if(multivariate){
	        sigma.resize(ncv*(ncv+1)/2);
	        Matrix<double> upper(ncv,ncv);
	        Matrix<double> lower(ncv,ncv);
		for (unsigned i=0;i<ncv;i++){
	              for (unsigned j=0;j<ncv-i;j++){
	                      ifile->scanField("sigma_"+tmpvalues[j+i].getName()+"_"+tmpvalues[j].getName(),lower(j+i,j));
	                      upper(j,j+i)=lower(j+i,j);
	              }
	         }
	        Matrix<double> mymult(ncv,ncv);       
	        Matrix<double> invmatrix(ncv,ncv);       
	        mult(lower,upper,mymult);          
	        // now invert and get the sigmas
	        Invert(mymult,invmatrix);
	        // put the sigmas in the usual order 
	        unsigned k=0;
		for (unsigned i=0;i<ncv;i++){
			for (unsigned j=i;j<ncv;j++){
				sigma[k]=invmatrix(i,j);
				k++;
			}
		}
	  }else{
	  	for(unsigned i=0;i<ncv;++i)ifile->scanField("sigma_"+tmpvalues[i].getName(),sigma[i]);
	  }
          // height and biasfactor  
          double height;
          ifile->scanField("height",height);
	  ifile->scanField("biasf",dummy);
	  if(ifile->FieldExist("clock")) ifile->scanField("clock",dummy);
          // check that all the fields are read
          ifile->scanField();
          // now put all in a structure or vector  
          // make the simple case, not with grids now 
          hills.push_back(Gaussian(center,sigma,height,multivariate));
          nhills++;
      }
   }
   else { 
     cerr<<"FILE NOT EXISTING"<<endl;
   }
   cerr<<"FOUND "<<hills.size()<<" HILLS "<<endl;
   // find min and max
  cerr<<"end of sum_hills"<<endl;
  return 0;
}

PLUMED_REGISTER_CLTOOL(CLToolSumHills,"sum_hills")

}
