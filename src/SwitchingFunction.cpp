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
#include "SwitchingFunction.h"
#include "Tools.h"
#include "Keywords.h"
#include <vector>
#include <limits>

using namespace std;
using namespace PLMD;

std::string SwitchingFunction::documentation(){
  std::ostringstream ostr;
  ostr<<"A switching function is a function that is one for values less than \\f$d_0\\f$ and which smoothly decays to zero for values greater than \\f$d_0\\f$. "; 
  ostr<<"Within plumed you can use the following functions to control the decay of the switching function ";
  ostr<<NonLinearFunction::printLatexFunctions( "\\frac{ r - d_0 }{ r_0 }", true )<<" ";
  ostr<<NonLinearFunction::printInput( "R_0=\\f$r_0\\f$ D_0=\\f$d_0\\f$", true )<<" "; 
  ostr<<"For all the switching functions if the D_0 is missing the \\f$d_0\\f$ parameter is assumed equal to zero.";
  ostr<<"You can also add the D_MAX flag to all switching function definitions to specify that if the input \\f$r\\f$ ";
  ostr<<"is greater than this value the value of the switching function is very small and can thus be assumed to equal zero.";

//  ostr<<"\\f$s(r)=\\frac{ 1 - \\left(\\frac{ r - d_0 }{ r_0 }\\right)^{n} }{ 1 - \\left(\\frac{ r - d_0 }{ r_0 }\\right)^{m} } \\f$, ";
//  ostr<<"\\f$s(r)=\\exp\\left(-\\frac{ r - d_0 }{ r_0 }\\right)\\f$ or using \\f$s(r)=\\exp\\left(-\\frac{ (r - d_0)^2 }{ 2r_0^2 }\\right)\\f$. ";
//  ostr<<"The first of these options is specified using the syntax {RATIONAL R_0=\\f$r_0\\f$ D_0=\\f$d_0\\f$ NN=\\f$n\\f$ MM=\\f$m\\f$} and if ";
//  ostr<<"the D_0, NN and MM keywords are missing they are assumed equal to 0, 6 and 12 respectively.  The second form is specified using ";
//  ostr<<"{EXP R_0=\\f$r_0\\f$ D_0=\\f$d_0\\f$} and if the D_0 is missing it is assumed equal to 0.  The third form is specified using ";
//  ostr<<"{GAUSSIAN R_0=\\f$r_0\\f$ D_0=\\f$d_0\\f$} and if the D_0 is missing it is assumed equal to 0. You can add the D_MAX flag to ";
//  ostr<<"all switching function definitions to specify that if the input \\f$r\\f$ is greater than this value the value of the switching ";
//  ostr<<"function is very small and can thus be assumed to equal zero.";
  return ostr.str();
}

void SwitchingFunction::set(const std::string & definition){
  // Get the name of the non linear function we are using
  vector<string> data=Tools::getWords(definition);
  plumed_assert(data.size()>=1);
  string name=data[0];
  data.erase(data.begin());
  // Setup the non linear function
  nlfunc.set(name,data,true);

  init=true;
  dmax=std::numeric_limits<double>::max();

  // Now read the actual switching funciton parameters
  double r0;
  bool found_r0=Tools::parse(data,"R_0",r0);
  plumed_massert(found_r0,"R_0 is needed");
  invr0=1.0/r0;
  d0=0; Tools::parse(data,"D_0",d0);
  Tools::parse(data,"D_MAX",dmax);

  if( !data.empty() ){
      std::string errormsg="found the following rogue keywords in switching function input : ";
      for(unsigned i=0;i<data.size();++i) errormsg = errormsg + data[i] + " "; 
      plumed_merror(errormsg);
  }
}

std::string SwitchingFunction::description() const {
  std::ostringstream ostr;
  ostr<<1./invr0<<".  Using "<<nlfunc.getName()<<" swiching function with parameters d0="<<d0<<" "<<nlfunc.writeParameters();
  return ostr.str(); 
}

double SwitchingFunction::calculate(double distance,double&dfunc)const{
  plumed_massert(init,"you are trying to use an unset SwitchingFunction");
  const double rdist = (distance-d0)*invr0;
  double result;
  if(rdist<=0.){
     result=1.;
     dfunc=0.0;
  }else if(rdist>dmax){
     result=0.;
     dfunc=0.0;
  }else{
     if( nlfunc.inputXSquared() ){
         result=nlfunc.calculate( rdist*rdist, dfunc );
         dfunc*=rdist;
     } else { 
         result=nlfunc.calculate( rdist, dfunc );
     }
// this is for the chain rule:
    dfunc*=invr0;
// this is because calculate() sets dfunc to the derivative divided times the distance.
// (I think this is misleading and I would like to modify it - GB)
    dfunc/=distance;
  }
  return result;
}

SwitchingFunction::SwitchingFunction():
  init(false){
}

void SwitchingFunction::set(int nn,int mm,double r0,double d0){
  init=true;
  std::vector<std::string> data;
  std::string nn_str; Tools::convert(nn,nn_str); data.push_back( "NN="+nn_str );
  std::string mm_str; Tools::convert(mm,mm_str); data.push_back( "MM="+mm_str );
  nlfunc.set("RATIONAL",data,true);
  this->invr0=1.0/r0;
  this->d0=d0;
  this->dmax=pow(0.00001,1./(nn-mm));
}

double SwitchingFunction::get_r0() const {
  return 1./invr0;
}
