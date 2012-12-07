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
#include "Action.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
//#include "Plumed.h"
#include "PlumedCommunicator.h"
#include "Random.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include "PlumedFile.h"
#include "Value.h"
#include "Matrix.h"

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
// parse it as it was a restart
   PlumedIFile *ifile=new PlumedIFile();
   vector<string> myfields;
   vector<string> mycvs;
   vector<string> myperiod_min;
   vector<string> myperiod_max;

   if(ifile->FileExist(hillsFile)){
      ifile->open(hillsFile);
      ifile->scanFieldList(myfields);
      //for(int i=0;i<myfields.size();i++)cerr<<myfields[i]<<endl; 
      // now find first sigma 
      size_t found;
      bool before_sigma=true;
      for(int i=0;i<myfields.size();i++){
        size_t pos = 0;
	found=myfields[i].find("sigma_", pos);
        if (found!=string::npos)before_sigma=false;
        // cvs are after time and before sigmas 
        found=myfields[i].find("time", pos); 
        if( found==string::npos && before_sigma){
             mycvs.push_back(myfields[i]);
	     cerr<<"found variable "<<mycvs.back()<<endl;
             // get periodicity
             myperiod_min.push_back("none");
             myperiod_max.push_back("none");
             if(ifile->FieldExist("min_"+mycvs.back())){
			string val;
			ifile->scanField("min_"+mycvs.back(),val);
                        myperiod_min[myperiod_min.size()-1]=val; 
             }
 	     if(ifile->FieldExist("max_"+mycvs.back())){
			string val;
			ifile->scanField("max_"+mycvs.back(),val);
                        myperiod_max[myperiod_max.size()-1]=val; 
             }
             if(  myperiod_min.back()!=string("none") &&  myperiod_max.back()!=string("none")  ){
                   cerr<<"variable "<<mycvs.back()<<" is periodic in range "<<myperiod_min.back()<<" and "<<myperiod_max.back()<<endl;
             } 

        }
      }
      // is multivariate ???
      std::string sss;
      bool multivariate=false;
      ifile->scanField("multivariate",sss);
      if(sss=="true"){ cerr<<"IS A MULTIVARIATE CALC "<<endl; multivariate=true;}
      else if(sss=="false"){cerr<<"DIAGONAL HILLS"<<endl; multivariate=false;}
      else plumed_merror("cannot parse multivariate = "+ sss);
      ifile->close();
   }
   else { 
     cerr<<"FILE NOT EXISTING"<<endl;
   }

   PlumedMain plumed;
   std::string ss;
   unsigned nn=1;
   ss="setNatoms";
   plumed.cmd(ss,&nn);  
   ss="init";
   plumed.cmd("init",&nn);  
   for(int i=0;i<mycvs.size();i++){
        std::string actioninput; 
        //actioninput=std::string("DISTANCE ATOMS=1,2 LABEL=")+myfields[i];           //the CV 
        actioninput=std::string("FAKE  ATOMS=1 LABEL=")+mycvs[i];           //the CV 
        // periodicity
        if (myperiod_max[i]==string("none")){
		actioninput+=string(" PERIODIC=NO "); 
        }else{
		actioninput+=string(" PERIODIC=")+myperiod_min[i]+string(",")+myperiod_max[i]; 
        } 
   //     cerr<<"FAKELINE: "<<actioninput<<endl;
        plumed.readInputString(actioninput);
   }
   // define the metadynamics
   int ncv=mycvs.size();
   std::string actioninput=std::string("METAD ARG=");
   for(unsigned i=0;i<ncv-1;i++)actioninput+=std::string(mycvs[i])+",";
   actioninput+=myfields[ncv-1];
   actioninput+=std::string(" SIGMA=");
   for(unsigned i=1;i<ncv;i++)actioninput+=std::string("0.1,");
    actioninput+=std::string("0.1 HEIGHT=1.0 PACE=1 FILE=");
    actioninput+=hillsFile;
    // multivariate? welltemp? grids? restart from grid? automatically generate it? which projection? stride?    
    cerr<<"METASTRING:  "<<actioninput<<endl;
    plumed.readInputString(actioninput);
//      //plumed.cmd("calc");
  cerr<<"end of sum_hills"<<endl;
  return 0;
}

PLUMED_REGISTER_CLTOOL(CLToolSumHills,"sum_hills")

}
