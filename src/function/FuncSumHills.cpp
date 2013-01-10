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
#include "ActionRegister.h"
#include "Function.h"
#include "tools/Exception.h"

using namespace std;

namespace PLMD{
namespace function{


//+PLUMEDOC FUNCTION FUNCSUMHILLS 
/*

*/
//+ENDPLUMEDOC


class FuncSumHills :
  public Function
{
  vector<string> hillsFiles,histoFiles; 
  vector<string> proj; 
  unsigned initstride,stride;
  bool iscltool;
  double beta;
public:
  FuncSumHills(const ActionOptions&);
  ~FuncSumHills();
  void calculate(); // this probably is not needed
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FuncSumHills,"FUNCSUMHILLS")

void FuncSumHills::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); 
  keys.add("optional","HILLSFILES"," source file for hills creation(may be the same as HILLS)"); // this can be a vector! 
  keys.add("optional","HISTOFILES"," source file for histogram creation(may be the same as HILLS)"); // also this can be a vector!
  keys.add("optional","PROJ"," only with sumhills: the projection on the cvs");
  keys.add("optional","KT"," only with sumhills: the kt factor when projection on cvs");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid"); 
  keys.add("optional","OUTHILLS"," output file for hills ");
  keys.add("optional","INITSTRIDE"," stride if you want an initial dump ");
  keys.add("optional","STRIDE"," stride when you do it on the fly ");
  keys.addFlag("ISCLTOOL",true,"use via plumed commandline: calculate at read phase and then go");
}

FuncSumHills::FuncSumHills(const ActionOptions&ao):
Action(ao),
Function(ao),
initstride(-1),
stride(-1),
iscltool(false),
beta(-1.)
{
  histoFiles.push_back("COLVAR");
  hillsFiles.push_back("HILLS");
  // here read 
  // Grid Stuff
  vector<std::string> gmin(getNumberOfArguments());
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=getNumberOfArguments() && gmin.size()!=0) error("not enough values for GRID_MIN");
  vector<std::string> gmax(getNumberOfArguments());
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=getNumberOfArguments() && gmax.size()!=0) error("not enough values for GRID_MAX");
  vector<unsigned> gbin(getNumberOfArguments());
  parseVector("GRID_BIN",gbin);
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) error("not enough values for GRID_BIN");
  plumed_assert(gmin.size()==gmax.size() && gmin.size()==gbin.size());
  // add some automatic hills width: not in case stride is defined  
  // since when you start from zero the automatic size will be zero!
 
  // hills file: 
  parseVector("HILLSFILES",hillsFiles);
  // histo file: 
  parseVector("HISTOFILES",histoFiles);
  // needs a projection? 
  proj.clear();
  parseVector("PROJ",proj);
  plumed_massert(proj.size()<getNumberOfArguments()," The number of projection must be less than the full list of arguments ");
  string kt="";
  if(proj.size()>0 ) {
    parse("KT",kt);
    plumed_massert(kt!="","if you make a projection then you need KT flag!"); 
    Tools::convert(beta,kt); beta=1./beta; 
  }
  // is a cltool: then you start and then die
  parseFlag("ISCLTOOL",iscltool);
  //what might it be this? 
  checkRead();
  // here start 
  // want something right now?? do it and return
  // your argument is a set of cvs 
  // then you need: a hills / a colvar-like file (to do a histogram) 
  // create a bias representation for this
  if(iscltool){

    // check if the file exists 

    //  

    return;
  } 
  // just an initialization but you need to do something on the fly?: need to connect with a metad run and its grid representation 
  // your argument is a metad run
  // if the grid does not exist crash and say that you need some data 
  // otherwise just link with it

}

void FuncSumHills::calculate(){
  // this should be connected only with a grid representation to metadynamics 
  // at regular time just dump it
}

FuncSumHills::~FuncSumHills(){
}

}
}


