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
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC COLVAR FAKE 
/*
This is a fake colvar container used by cltools or various other actions
and just support input and period definition

\par Examples

FAKE ATOMS=1 PERIODIC=-3.14,3.14   LABEL=d2
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC
   
class ColvarFake : public Colvar {
  bool pbc;

public:
  static void registerKeywords( Keywords& keys );
  ColvarFake(const ActionOptions&);
// active methods:
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(ColvarFake,"FAKE")

void ColvarFake::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the pair of atom that we are calculating the distance between");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  keys.use("PERIODIC");
}

ColvarFake::ColvarFake(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  //assert(atoms.size()==2);
  addValueWithDerivatives(); 
  if( keywords.exists("PERIODIC") ){
     std::vector<std::string> period;
     parseVector("PERIODIC",period);
     if(period.size()==1 && period[0]=="NO"){
        setNotPeriodic();
        log.printf("  this variable is not periodic \n");
     } else if(period.size()==2){
        setPeriodic(period[0],period[1]);
        log.printf("  this variable is periodic with period %s : %s \n",period[0].c_str(),period[1].c_str());
     } else error("missing PERIODIC keyword");
  }
  checkRead();

  //log.printf("  between atoms %d %d\n",atoms[0].serial(),atoms[1].serial());

  //setNotPeriodic();
  //void Value::setDomain(const std::string& pmin,const std::string& pmax)

  requestAtoms(atoms);
}


// calculator
void ColvarFake::calculate(){

    setValue  (0.);

}

}



