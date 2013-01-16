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
#include "tools/Communicator.h"
#include "tools/Random.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include "tools/File.h"
#include "core/Value.h"
#include "tools/Matrix.h"

using namespace std;

namespace PLMD {

//+PLUMEDOC TOOLS sum_hills 
/*
sumhills is a tool that allows one to to use plumed to post-process an existing hills/colvar file 

a typical case is about the integration of a hills file: 

\verbatim
plumed sum_hills  --hills PATHTOMYHILLSFILE 
\endverbatim

The default name for the output file will be fes.dat 
Note that starting from this version plumed will automatically detect the 
number of the variables you have and their periodicity. 
Additionally, if you use flexible hills (multivariate gaussians), plumed will understand it from the HILLS file.

if you want to integrate out some variable you do

\verbatim
plumed sum_hills  --hills PATHTOMYHILLSFILE   --idw t1 --kt 0.6 
\endverbatim

where with --idw you define the variables that you want
all the others will be integrated out. --kt defines the temperature of the system in energy units.
(be consistent with the units you have in your hills: plumed will not check this for you)  
If you need more variables then you may use a comma separated syntax

\verbatim
plumed sum_hills  --hills PATHTOMYHILLSFILE   --idw t1,t2 --kt 0.6  
\endverbatim

You can define the output grid only with the number of bins you want 
while min/max will be detected for you

\verbatim
plumed sum_hills --bin 99,99 --hills PATHTOMYHILLSFILE   
\endverbatim

or full grid specification

\verbatim
plumed sum_hills --bin 99,99 --min -pi,-pi --max pi,pi --hills PATHTOMYHILLSFILE   
\endverbatim

You can of course use numbers instead of -pi/pi.

You can use a --stride keyword to have a dump each bunch of hills you read 
\verbatim
plumed sum_hills --stride 300 --hills PATHTOMYHILLSFILE   
\endverbatim
 
You can also have, in case of welltempered metadynamics, only the negative 
bias instead of the free energy through the keyword --negbias

\verbatim
plumed sum_hills --negbias --hills PATHTOMYHILLSFILE   
\endverbatim

Here the default name will be negativebias.dat 

From time to time you might need to use HILLS or a COLVAR file 
as it was just a simple set  of points from which you want to build 
a free energy by using -(1/beta)log(P)
then you use --histo

\verbatim
plumed sum_hills --histo PATHTOMYCOLVARORHILLSFILE  --sigma 0.2,0.2 --kt 0.6  
\endverbatim

in this case you need a --kt to do the reweighting and then you
need also some width (with the --sigma keyword) for the histogram calculation (actually will be done with 
gaussians, so it will be a continuous histogram)
Here the default output will be correction.dat.

*/
//+ENDPLUMEDOC

class CLToolSumHills : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  CLToolSumHills(const CLToolOptions& co );
  int main(FILE* in,FILE*out,Communicator& pc);
  string description()const;
};

void CLToolSumHills::registerKeywords( Keywords& keys ){
  CLTool::registerKeywords( keys );
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  keys.add("optional","--hills","specify the name of the hills file");
  keys.add("optional","--histo","specify the name of the file for histogram a colvar/hills file is good");
  keys.add("optional","--stride","specify the stride for integrating hills file (default 0=never)");
  keys.add("optional","--min","the lower bounds for the grid");
  keys.add("optional","--max","the upper bounds for the grid");
  keys.add("optional","--bin","the number of bins for the grid");
  keys.add("optional","--idw","specify the variables to be integrated (default is all)");
  keys.add("optional","--outfile","specify the outputfile for sumhills");
  keys.add("optional","--outhisto","specify the outputfile for the histogram");
  keys.add("optional","--kt","specify temperature for integrating out variables");
  keys.add("optional","--sigma"," a vector that specify the sigma for binning (only needed when doing histogram ");
  keys.addFlag("--negbias",false," print the negative bias instead of the free energy (only needed with welltempered runs and flexible hills) ");
}

CLToolSumHills::CLToolSumHills(const CLToolOptions& co ):
CLTool(co)
{
 inputdata=commandline;
}

string CLToolSumHills::description()const{ return "sum the hills with  plumed"; }

int CLToolSumHills::main(FILE* in,FILE*out,Communicator& pc){
  
// Read the hills input file name  
  vector<string> hillsFiles; 
  bool dohills;
  dohills=parseVector("--hills",hillsFiles);
// Read the histogram file
  vector<string> histoFiles; 
  bool dohisto;
  dohisto=parseVector("--histo",histoFiles);

  plumed_massert(dohisto || dohills,"you should use --histo or/and --hills command");

  vector<string> vcvs;
  vector<string> vpmin;
  vector<string> vpmax;
  bool vmultivariate;
  if(dohills){
       // parse it as it was a restart
       IFile *ifile=new IFile();
       ifile->findCvsAndPeriodic(hillsFiles[0], vcvs, vpmin, vpmax, vmultivariate);
       free(ifile);
  }

  vector<string> hcvs;
  vector<string> hpmin;
  vector<string> hpmax;
  bool hmultivariate;
 
  vector<std::string> sigma; 
  if(dohisto){
       IFile *ifile=new IFile();
       ifile->findCvsAndPeriodic(histoFiles[0], hcvs, hpmin, hpmax, hmultivariate);
       free(ifile);
       // here need also the vector of sigmas
       parseVector("--sigma",sigma);
       if(sigma.size()!=hcvs.size())plumed_merror("you should define --sigma vector when using histogram");
  }

  if(dohisto && dohills){
    plumed_massert(vcvs==hcvs,"variables for histogram and bias should have the same labels"); 
    plumed_massert(hpmin==vpmin,"variables for histogram and bias should have the same min for periodicity"); 
    plumed_massert(hpmax==vpmax,"variables for histogram and bias should have the same max for periodicity"); 
  }

  // now put into a neutral vector  
  
  vector<string> cvs;
  vector<string> pmin;
  vector<string> pmax;

  if(dohills){
    cvs=vcvs;
    pmin=vpmin;
    pmax=vpmax;
  }
  if(dohisto){
    cvs=hcvs;
    pmin=hpmin;
    pmax=hpmax;
  }

  // setup grids
  unsigned grid_check=0; 
  vector<std::string> gmin(cvs.size());
  if(parseVector("--min",gmin)){
  	if(gmin.size()!=cvs.size() && gmin.size()!=0) plumed_merror("not enough values for --min");
       grid_check++;
  }
  vector<std::string> gmax(cvs.size() );
  if(parseVector("--max",gmax)){
  	if(gmax.size()!=cvs.size() && gmax.size()!=0) plumed_merror("not enough values for --max");
       grid_check++;
  }
  vector<std::string> gbin(cvs.size());
  bool grid_has_bin; grid_has_bin=false;	
  if(parseVector("--bin",gbin)){
       if(gbin.size()!=cvs.size() && gbin.size()!=0) plumed_merror("not enough values for --bin");
       grid_has_bin=true;
  }
  // allowed: no grids only bin
  // not allowed: partial grid definition 
  plumed_massert( gmin.size()==gmax.size() && (gmin.size()==0 ||  gmin.size()==cvs.size() ) ,"you should specify --min and --max together with same number of components");
  //plumed_massert(( (grid_check==0 && grid_has_bin==false ) || (grid_check==2 && grid_has_bin==true) ),"you should define all the --min --max --bin keys");

  PlumedMain plumed;
  std::string ss;
  unsigned nn=1;
  ss="setNatoms";
  plumed.cmd(ss,&nn);  
  ss="init";
  plumed.cmd("init",&nn);  
  for(int i=0;i<cvs.size();i++){
       std::string actioninput; 
       actioninput=std::string("FAKE  ATOMS=1 LABEL=")+cvs[i];           //the CV 
       // periodicity
       if (pmax[i]==string("none")){
       	actioninput+=string(" PERIODIC=NO "); 
       }else{
       	actioninput+=string(" PERIODIC=")+pmin[i]+string(",")+pmax[i]; 
               // check if min and max values are ok with grids
               if(grid_check==2){  
                   double gm; Tools::convert(gmin[i],gm);              
                   double pm; Tools::convert(pmin[i],pm);              
                   if(  gm<pm ){
                        plumed_merror("Periodicity issue : GRID_MIN value ( "+gmin[i]+" ) is less than periodicity in HILLS file in "+cvs[i]+ " ( "+pmin[i]+" ) ");
                   } 
                   Tools::convert(gmax[i],gm);              
                   Tools::convert(pmax[i],pm);              
                   if(  gm>pm ){
                        plumed_merror("Periodicity issue : GRID_MAX value ( "+gmax[i]+" ) is more than periodicity in HILLS file in "+cvs[i]+ " ( "+pmax[i]+" ) ");
                   }
               } 
       } 
       plumed.readInputString(actioninput);
  }
  // define the metadynamics
  unsigned ncv=cvs.size();
  std::string actioninput;
  vector<std::string> idw;
  // check if the variables to be used are correct 
  if(parseVector("--idw",idw)){
      for(unsigned i=0;i<idw.size();i++){
          bool found=false;
          for(unsigned j=0;j<cvs.size();j++){
                if(idw[i]==cvs[j])found=true;
          }
          if(!found)plumed_merror("variable "+idw[i]+" is not found in the bunch of cvs: revise your --idw option" ); 
      } 
      plumed_massert( idw.size()<=cvs.size() ,"the number of variables to be integrated should be at most equal to the total number of cvs  "); 
      // in this case you neeed a beta factor!
  } 

  std::string kt; kt=std::string("1.");// assign an arbitrary value just in case that idw.size()==cvs.size() 
  if ( dohisto || idw.size()!=0  ) {
  		plumed_massert(parse("--kt",kt),"if you make a dimensionality reduction (--idw) or a histogram (--histo) then you need to define --kt ");
  }

  /*

	different implementation through function

  */

  actioninput="FUNCSUMHILLS ISCLTOOL ARG=";
  for(unsigned i=0;i<(ncv-1);i++)actioninput+=std::string(cvs[i])+",";
  actioninput+=cvs[ncv-1];
  if(dohills){ actioninput+=" HILLSFILES=";
  for(unsigned i=0;i<hillsFiles.size()-1;i++)actioninput+=hillsFiles[i]+",";
     actioninput+=hillsFiles[hillsFiles.size()-1];
  // set the grid 
  }
  if(grid_check==2){
     actioninput+=std::string(" GRID_MAX=");
     for(unsigned i=0;i<(ncv-1);i++)actioninput+=gmax[i]+",";
     actioninput+=gmax[ncv-1];
     actioninput+=std::string(" GRID_MIN=");
     for(unsigned i=0;i<(ncv-1);i++)actioninput+=gmin[i]+",";
     actioninput+=gmin[ncv-1];
  }
  if(grid_has_bin){
     actioninput+=std::string(" GRID_BIN=");
     for(unsigned i=0;i<(ncv-1);i++)actioninput+=gbin[i]+",";
     actioninput+=gbin[ncv-1];
  }
  std::string  stride; stride="";
  if(parse("--stride",stride)){
    actioninput+=std::string(" INITSTRIDE=")+stride;
  }
  if(idw.size()!=0){ 
     actioninput+=std::string(" PROJ=");
     for(unsigned i=0;i<idw.size()-1;i++){actioninput+=idw[i]+",";}
     actioninput+=idw.back();  
  }
  if(idw.size()!=0 || dohisto){
     actioninput+=std::string(" KT=")+kt ;
  } 
  if(dohisto){
	actioninput+=" HISTOFILES=";
        for(unsigned i=0;i<histoFiles.size()-1;i++){actioninput+=histoFiles[i]+",";}
        actioninput+=histoFiles[histoFiles.size()-1];
 
        actioninput+=std::string(" HISTOSIGMA=");
        for(unsigned i=0;i<sigma.size()-1;i++){actioninput+=sigma[i]+",";}
        actioninput+=sigma.back();  
  } 
 
  bool negbias;
  parseFlag("--negbias",negbias);
  if(negbias){
 	actioninput+=" NEGBIAS ";
  }

  //cerr<<"FUNCSTRING:  "<<actioninput<<endl;
  plumed.readInputString(actioninput);
  // if not a grid, then set it up automatically
  return 0;
}

PLUMED_REGISTER_CLTOOL(CLToolSumHills,"sum_hills")



}
