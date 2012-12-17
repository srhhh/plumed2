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
#include "tools/Tools.h"
#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/PlumedCommunicator.h"
#include "tools/Random.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include "tools/PlumedFile.h"
#include "core/Value.h"
#include "tools/Matrix.h"

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
  keys.add("optional","--stride","specify the stride for integrating hills file (default 0=never)");
  keys.add("compulsory","--GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","--GRID_MAX","the upper bounds for the grid");
  keys.add("compulsory","--GRID_BIN","the number of bins for the grid");
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
	     cerr<<"found variable number  "<<mycvs.size()<<" :  "<<mycvs.back()<<endl;
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

   // setup grids
   vector<std::string> gmin(mycvs.size());
   parseVector("--GRID_MIN",gmin);
   //std::string myerr;myerr="not enough values for GRID_MIN: needed "+std::string(mycvs.size())
   if(gmin.size()!=mycvs.size() && gmin.size()!=0) plumed_merror("not enough values for --GRID_MIN");
   vector<std::string> gmax(mycvs.size() );
   parseVector("--GRID_MAX",gmax);
   if(gmax.size()!=mycvs.size() && gmax.size()!=0) plumed_merror("not enough values for --GRID_MAX");
   vector<std::string> gbin(mycvs.size());
   parseVector("--GRID_BIN",gbin);
   if(gbin.size()!=mycvs.size() && gbin.size()!=0) plumed_merror("not enough values for --GRID_BIN");
   plumed_assert(gmin.size()==gmax.size() && gmin.size()==gbin.size());
   // TODO set automatic grid

   PlumedMain plumed;
   std::string ss;
   unsigned nn=1;
   ss="setNatoms";
   plumed.cmd(ss,&nn);  
   ss="init";
   plumed.cmd("init",&nn);  
   // it is a restart with HILLS  
   plumed.readInputString(string("RESTART"));
   for(int i=0;i<mycvs.size();i++){
        std::string actioninput; 
        //actioninput=std::string("DISTANCE ATOMS=1,2 LABEL=")+myfields[i];           //the CV 
        actioninput=std::string("FAKE  ATOMS=1 LABEL=")+mycvs[i];           //the CV 
        // periodicity
        if (myperiod_max[i]==string("none")){
		actioninput+=string(" PERIODIC=NO "); 
        }else{
		actioninput+=string(" PERIODIC=")+myperiod_min[i]+string(",")+myperiod_max[i]; 
                // check if min and max values are ok
                double gm; Tools::convert(gmin[i],gm);              
                double pm; Tools::convert(myperiod_min[i],pm);              
                if(  gm<pm ){
                     plumed_merror("Periodicity issue : GRID_MIN value ( "+gmin[i]+" ) is less than periodicity in HILLS file in "+mycvs[i]+ " ( "+myperiod_min[i]+" ) ");
                } 
                Tools::convert(gmax[i],gm);              
                Tools::convert(myperiod_max[i],pm);              
		if(  gm>pm ){
                     plumed_merror("Periodicity issue : GRID_MAX value ( "+gmax[i]+" ) is more than periodicity in HILLS file in "+mycvs[i]+ " ( "+myperiod_max[i]+" ) ");
                } 
        } 
   //     cerr<<"FAKELINE: "<<actioninput<<endl;
        plumed.readInputString(actioninput);
   }
   // define the metadynamics
   unsigned ncv=mycvs.size();
   std::string actioninput=std::string("METAD ARG=");
   for(unsigned i=0;i<(ncv-1);i++)actioninput+=std::string(mycvs[i])+",";
   actioninput+=mycvs[ncv-1];
   actioninput+=std::string(" SIGMA=");
   for(unsigned i=1;i<ncv;i++)actioninput+=std::string("0.1,");
    actioninput+=std::string("0.1 HEIGHT=1.0 PACE=1 FILE=");
    // this sets the restart 
    actioninput+=hillsFile;

    // set the grid 
    actioninput+=std::string(" GRID_MAX=");
    for(unsigned i=0;i<(ncv-1);i++)actioninput+=gmax[i]+",";
    actioninput+=gmax[ncv-1];
    actioninput+=std::string(" GRID_MIN=");
    for(unsigned i=0;i<(ncv-1);i++)actioninput+=gmin[i]+",";
    actioninput+=gmin[ncv-1];
    actioninput+=std::string(" GRID_BIN=");
    for(unsigned i=0;i<(ncv-1);i++)actioninput+=gbin[i]+",";
    actioninput+=gbin[ncv-1];
    actioninput+=std::string(" GRID_WFILE=fes.dat  ");
    // take the stride 
    std::string  stride; 
    if(parse("--stride",stride)){
   	 actioninput+=std::string("  GRID_WSTRIDE=")+stride;
    }
    // the input keyword
    actioninput+=std::string(" SUMHILLS ");
    // multivariate? welltemp? grids? restart from grid? automatically generate it? which projection? stride?    
    cerr<<"METASTRING:  "<<actioninput<<endl;
    plumed.readInputString(actioninput);

    // if not a grid, then set it up automatically
  cerr<<"end of sum_hills"<<endl;
  return 0;
}

PLUMED_REGISTER_CLTOOL(CLToolSumHills,"sum_hills")

}
